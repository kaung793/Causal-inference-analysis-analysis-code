file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R")); source(file.path(root, "R", "mi_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "08_equal_lag_sensitivity")
raw <- read_analysis_file(cfg$equal_lag_long)
d <- prepare_primary(raw)
# Normalize labels from either D1-D3 or D1_to_D3 source conventions by start day.
d$window <- c(`1` = "D1_to_D3", `3` = "D3_to_D5", `5` = "D5_to_D7")[as.character(d$day)]
d$window <- factor(d$window, levels = c("D1_to_D3", "D3_to_D5", "D5_to_D7"))
if (isTRUE(cfg$strict_cohort_checks) && (nrow(d) != 2403L || length(unique(d$subject_id)) != 801L)) stop("Equal-lag source must contain 2,403 windows from 801 patients")

fit_equal_model <- function(data, outcome) {
  z <- scale_covariates(data, c("iap_current", "apache_ii", "age", "creatinine", "map"))
  needed <- c(outcome, "estimated_fluid_balance_l", "window", "iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z", "subject_id")
  z <- z[complete.cases(z[needed]), , drop = FALSE]
  rhs_full <- "estimated_fluid_balance_l * window + iap_current_z + apache_ii_z + age_z + sex + etiology + creatinine_z + map_z + (1 | subject_id)"
  rhs_reduced <- "estimated_fluid_balance_l + window + iap_current_z + apache_ii_z + age_z + sex + etiology + creatinine_z + map_z + (1 | subject_id)"
  binary <- outcome != "iap_next"
  full <- if (binary) fit_glmm_strict(as.formula(paste(outcome, "~", rhs_full)), z) else fit_lmm_strict(as.formula(paste(outcome, "~", rhs_full)), z)
  reduced <- if (binary) fit_glmm_strict(as.formula(paste(outcome, "~", rhs_reduced)), z) else fit_lmm_strict(as.formula(paste(outcome, "~", rhs_reduced)), z)
  rows <- list(); main <- "estimated_fluid_balance_l"
  for (w in levels(z$window)) {
    it <- paste0(main, ":window", w); terms <- setNames(1, main)
    if (w != levels(z$window)[1L]) terms <- c(terms, setNames(1, it))
    e <- linear_contrast(full, terms, exponentiate = binary)
    rows[[w]] <- data.frame(outcome = outcome, window = w, result_type = "window_specific",
                            estimate = e["estimate"], ci_lower = e["ci_lower"], ci_upper = e["ci_upper"], p_value = e["p_value"],
                            linear_estimate = e["linear_estimate"], variance = e["variance"],
                            effect_measure = ifelse(binary, "OR per 1,000 mL", "beta mmHg per 1,000 mL"),
                            windows = nrow(z), patients = length(unique(z$subject_id)))
    if (w != levels(z$window)[1L]) {
      c0 <- linear_contrast(full, setNames(1, it), exponentiate = binary)
      rows[[paste0(w, "_contrast")]] <- data.frame(outcome = outcome, window = w, result_type = "contrast_vs_D1_to_D3",
                            estimate = c0["estimate"], ci_lower = c0["ci_lower"], ci_upper = c0["ci_upper"], p_value = c0["p_value"],
                            linear_estimate = c0["linear_estimate"], variance = c0["variance"],
                            effect_measure = ifelse(binary, "ratio of ORs", "delta beta"), windows = nrow(z), patients = length(unique(z$subject_id)))
    }
  }
  lrt <- suppressMessages(anova(reduced, full))
  global <- data.frame(outcome = outcome, statistic = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2],
                       windows = nrow(z), patients = length(unique(z$subject_id)))
  list(fit = full, reduced = reduced, rows = do.call(rbind, rows), global = global,
       diagnostics = model_diagnostics(full, paste0("equal_lag_", outcome)), data = z)
}

outcomes <- c("iap_next", "iap15_next", "iap20_next")
cc <- lapply(outcomes, function(o) fit_equal_model(d, o)); names(cc) <- outcomes
cc_rows <- do.call(rbind, lapply(cc, `[[`, "rows")); cc_rows$analysis <- "complete_case_Model3"
cc_global <- do.call(rbind, lapply(cc, `[[`, "global")); cc_global$analysis <- "complete_case_Model3"
cc_diag <- do.call(rbind, lapply(cc, `[[`, "diagnostics"))
write_csv_atomic(cc_rows, file.path(out, "complete_case_results.csv")); write_csv_atomic(cc_global, file.path(out, "complete_case_global_interaction.csv"))
write_csv_atomic(cc_diag, file.path(out, "complete_case_diagnostics.csv")); saveRDS(lapply(cc, `[[`, "fit"), file.path(out, "complete_case_fits.rds"), compress = "xz")

if (isTRUE(cfg$run_equal_lag_mi)) {
  base <- d[!duplicated(d$subject_id), c("subject_id", "age", "sex", "etiology")]
  tv <- c("fluid_intake_ml", "estimated_fluid_balance_ml", "iap_current", "iap_next", "apache_ii", "creatinine", "map")
  wide <- base
  for (w in levels(d$window)) {
    z <- d[d$window == w, c("subject_id", tv), drop = FALSE]
    names(z)[-1L] <- paste0(tv, "__", w)
    wide <- merge(wide, z, by = "subject_id", all.x = TRUE, sort = FALSE)
  }
  targets <- grep("^(apache_ii|creatinine|map)__", names(wide), value = TRUE)
  imp <- make_patient_wide_imputation(wide, targets, cfg$mi_m, cfg$equal_lag_mi_maxit, cfg$seed, cfg$mi_donors)
  saveRDS(imp, file.path(out, "equal_lag_mice.rds"), compress = "xz")
  to_long <- function(wide_complete) {
    parts <- lapply(levels(d$window), function(w) {
      z <- data.frame(subject_id = wide_complete$subject_id, window = w, age = wide_complete$age,
                      sex = wide_complete$sex, etiology = wide_complete$etiology, stringsAsFactors = FALSE)
      for (v in tv) z[[v]] <- wide_complete[[paste0(v, "__", w)]]
      z
    })
    z <- do.call(rbind, parts); z$window <- factor(z$window, levels = levels(d$window)); z$sex <- factor(z$sex); z$etiology <- factor(z$etiology)
    z$estimated_fluid_balance_l <- z$estimated_fluid_balance_ml / 1000
    z$iap15_next <- as.integer(z$iap_next >= 15); z$iap20_next <- as.integer(z$iap_next >= 20)
    z
  }
  mi_fits <- setNames(lapply(outcomes, function(x) vector("list", cfg$mi_m)), outcomes)
  mi_rows <- list(); qidx <- 0L
  for (i in seq_len(cfg$mi_m)) {
    li <- to_long(mice::complete(imp, i))
    for (o in outcomes) {
      ans <- fit_equal_model(li, o); mi_fits[[o]][[i]] <- ans$fit
      qidx <- qidx + 1L; mi_rows[[qidx]] <- transform(ans$rows, imputation = i)
    }
  }
  per_imp <- do.call(rbind, mi_rows)
  pooled <- rubin_pool_grouped(per_imp, c("outcome", "window", "result_type"), "linear_estimate", "variance")
  pooled$binary <- pooled$outcome != "iap_next"
  pooled$effect_measure <- ifelse(pooled$binary & pooled$result_type == "window_specific", "MI-pooled OR per 1,000 mL",
                                  ifelse(pooled$binary, "MI-pooled ratio of ORs", ifelse(pooled$result_type == "window_specific", "MI-pooled beta", "MI-pooled delta beta")))
  pooled$estimate_display <- ifelse(pooled$binary, exp(pooled$estimate), pooled$estimate)
  pooled$ci_lower_display <- ifelse(pooled$binary, exp(pooled$ci_lower), pooled$ci_lower)
  pooled$ci_upper_display <- ifelse(pooled$binary, exp(pooled$ci_upper), pooled$ci_upper)
  globals <- list()
  for (o in outcomes) {
    fits <- mi_fits[[o]]; ints <- grep("estimated_fluid_balance_l:window", names(fixed_coef(fits[[1L]])), value = TRUE, fixed = TRUE)
    Q <- t(vapply(fits, function(f) fixed_coef(f)[ints], numeric(length(ints))))
    U <- array(NA_real_, c(length(ints), length(ints), length(fits)))
    for (i in seq_along(fits)) U[, , i] <- vcov(fits[[i]])[ints, ints, drop = FALSE]
    pv <- pool_vector_rubin(Q, U); stat <- as.numeric(t(pv$estimate) %*% qr.solve(pv$total, pv$estimate))
    globals[[o]] <- data.frame(outcome = o, statistic = stat, df = length(ints), p_value = pchisq(stat, length(ints), lower.tail = FALSE), m = length(fits))
  }
  write_csv_atomic(per_imp, file.path(out, "mi_estimates_by_imputation.csv")); write_csv_atomic(pooled, file.path(out, "mi_rubin_pooled_results.csv"))
  write_csv_atomic(do.call(rbind, globals), file.path(out, "mi_global_interaction.csv")); saveRDS(mi_fits, file.path(out, "mi_fits.rds"), compress = "xz")
}
write_run_metadata(cfg, cfg$equal_lag_long, out)
