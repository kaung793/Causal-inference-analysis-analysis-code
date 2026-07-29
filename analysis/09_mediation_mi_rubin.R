file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "mi_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_mediation) || cli$force, "MI/Rubin mediation")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "09_mediation_mi_rubin")
require_packages(c("mice", "mediation"))
d <- read_analysis_file(cfg$mediation_patient)
aliases <- list(subject_id = c("subject_id", "id"), fluid_auc = c("fluid_auc", "early_fluid_auc"),
                iap_day5 = c("iap_day5"), iap_day1 = c("iap_day1", "iap_baseline"), apache = c("apache", "apache_ii"),
                creatinine = c("creatinine", "cr"), map = c("map"), death_hospital = c("death_hospital", "death_in_hosp"),
                age = "age", sex = "sex", etiology = "etiology")
d <- rename_aliases(d, aliases)
required <- c("subject_id", "fluid_auc", "iap_day5", "iap_day1", "age", "sex", "etiology", "apache", "creatinine", "map", "death_hospital")
assert_columns(d, required, "mediation patient data")
d$subject_id <- as.character(d$subject_id); d$sex <- factor(d$sex); d$etiology <- factor(d$etiology); d$death_hospital <- as.integer(d$death_hospital)
if (isTRUE(cfg$strict_cohort_checks)) {
  if (nrow(d) != 801L || anyDuplicated(d$subject_id)) stop("Mediation data must have one row for each of 801 patients")
  if (sum(d$death_hospital, na.rm = TRUE) != 145L) stop("Expected 145 in-hospital deaths")
  if (anyNA(d[c("fluid_auc", "iap_day5", "death_hospital")])) stop("Exposure, mediator, and outcome must not be imputed")
}

ini <- mice::mice(d, m = 1, maxit = 0, printFlag = FALSE)
method <- ini$method; method[] <- ""; targets <- c("apache", "creatinine", "map"); method[targets] <- "pmm"
pred <- ini$predictorMatrix; pred[,] <- 0
predictors <- setdiff(required, "subject_id")
for (target in targets) pred[target, setdiff(predictors, target)] <- 1
pred[, "subject_id"] <- 0; pred["subject_id", ] <- 0
imp <- mice::mice(d, m = cfg$mi_m, maxit = cfg$mi_maxit, method = method, predictorMatrix = pred,
                  donors = cfg$mi_donors, seed = cfg$seed, printFlag = FALSE,
                  remove.constant = FALSE, remove.collinear = FALSE)
saveRDS(imp, file.path(out, "mice_mediation.rds"), compress = "xz")
auc_median <- median(d$fluid_auc)

run_one <- function(z, binary, imputation) {
  z$exposure <- if (binary) as.integer(z$fluid_auc >= auc_median) else z$fluid_auc / 1000
  mm <- lm(iap_day5 ~ exposure + iap_day1 + age + sex + etiology + apache + creatinine + map, data = z)
  om <- glm(death_hospital ~ exposure + iap_day5 + iap_day1 + age + sex + etiology + apache + creatinine + map,
            data = z, family = binomial())
  set.seed(cfg$seed + imputation + ifelse(binary, 10000L, 0L))
  med <- mediation::mediate(mm, om, treat = "exposure", mediator = "iap_day5", boot = TRUE,
                            sims = cfg$mediation_sims, boot.ci.type = "perc", long = TRUE)
  draws <- cbind(ACME = as.numeric(med$d.avg.sims), ADE = as.numeric(med$z.avg.sims), Total = as.numeric(med$tau.sims))
  estimates <- c(ACME = med$d.avg, ADE = med$z.avg, Total = med$tau.coef)
  list(
    row = data.frame(imputation = imputation,
                     exposure = ifelse(binary, "High vs low at cohort median", "Per 1,000 mL-day"),
                     effect = names(estimates), estimate = estimates,
                     within_variance = diag(var(draws))),
    covariance = var(draws), models = list(mediator = mm, outcome = om)
  )
}

runs <- list(); rows <- list(); idx <- 0L
for (i in seq_len(cfg$mi_m)) {
  z <- mice::complete(imp, i); z$sex <- factor(z$sex, levels = levels(d$sex)); z$etiology <- factor(z$etiology, levels = levels(d$etiology))
  for (binary in c(FALSE, TRUE)) { idx <- idx + 1L; runs[[idx]] <- run_one(z, binary, i); rows[[idx]] <- runs[[idx]]$row }
}
per_imp <- do.call(rbind, rows); pooled_scalar <- rubin_pool_grouped(per_imp, c("exposure", "effect"), "estimate", "within_variance")

prop_rows <- list()
for (exposure_label in unique(per_imp$exposure)) {
  selected <- which(vapply(runs, function(x) x$row$exposure[1] == exposure_label, logical(1)))
  Q <- t(vapply(runs[selected], function(x) setNames(x$row$estimate, x$row$effect)[c("ACME", "ADE", "Total")], numeric(3)))
  U <- array(NA_real_, c(3, 3, length(selected)))
  for (j in seq_along(selected)) U[, , j] <- runs[[selected[j]]]$covariance[c("ACME", "ADE", "Total"), c("ACME", "ADE", "Total")]
  pv <- pool_vector_rubin(Q, U); A <- pv$estimate["ACME"]; T <- pv$estimate["Total"]
  gradient <- c(1 / T, 0, -A / T^2); variance <- as.numeric(t(gradient) %*% pv$total %*% gradient)
  estimate <- A / T; se <- sqrt(variance); bounded <- bounded_proportion_ci(estimate, se)
  prop_rows[[exposure_label]] <- data.frame(exposure = exposure_label, effect = "Proportion mediated", m = length(selected),
      estimate = estimate, within_variance = NA_real_, between_variance = NA_real_, total_variance = variance,
      standard_error = se, df = Inf, ci_lower = bounded[1], ci_upper = bounded[2], p_value = NA_real_,
      fraction_missing_information = NA_real_, ci_method = "delta method on logit scale using pooled ACME-total covariance")
}
pooled_scalar$ci_method <- "Rubin total variance"
pooled <- rbind(pooled_scalar, do.call(rbind, prop_rows)); row.names(pooled) <- NULL
write_csv_atomic(per_imp, file.path(out, "mediation_estimates_by_imputation.csv"))
write_csv_atomic(pooled, file.path(out, "mediation_mi_rubin_pooled_results.csv"))
write_csv_atomic(data.frame(variable = names(d), missing_n = vapply(d, function(x) sum(is.na(x)), integer(1))), file.path(out, "missingness.csv"))
write_run_metadata(cfg, cfg$mediation_patient, out, list(m = cfg$mi_m, maxit = cfg$mi_maxit, simulations_per_imputation = cfg$mediation_sims, auc_median_ml_day = auc_median))
