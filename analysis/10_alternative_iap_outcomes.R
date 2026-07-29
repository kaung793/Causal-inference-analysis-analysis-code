file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_revision_sensitivities) || cli$force, "alternative IAP outcomes")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "10_alternative_iap_outcomes")
require_packages(c("geepack", "ordinal", "sandwich"))
d <- prepare_primary(read_analysis_file(cfg$primary_long))
et <- suppressWarnings(as.numeric(as.character(d$etiology)))
d$etiology <- factor(ifelse(et %in% c(4, 5), 4, et))
d$wsacs_grade <- cut(d$iap_next, breaks = c(-Inf, 12, 16, 21, 26, Inf), right = FALSE,
                     labels = c("0", "I", "II", "III", "IV"), ordered_result = TRUE)
d <- scale_covariates(d, c("iap_current", "apache_ii", "age", "creatinine", "map"))
needed <- c("subject_id", "estimated_fluid_balance_l", "iap12_next", "wsacs_grade", "iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z")
z <- droplevels(d[complete.cases(d[needed]), , drop = FALSE])
fixed_rhs <- "estimated_fluid_balance_l + iap_current_z + apache_ii_z + age_z + sex + etiology + creatinine_z + map_z"

gee <- geepack::geeglm(as.formula(paste("iap12_next ~", fixed_rhs)), id = subject_id, data = z,
                       family = binomial(), corstr = "independence", std.err = "san.se")
gc <- summary(gee)$coefficients["estimated_fluid_balance_l", ]
rows <- data.frame(
  analysis = "IAP >=12 mmHg, Model 3 GEE", effect_measure = "OR per 1,000 mL",
  estimate = unname(exp(gc["Estimate"])), ci_lower = unname(exp(gc["Estimate"] - 1.96 * gc["Std.err"])),
  ci_upper = unname(exp(gc["Estimate"] + 1.96 * gc["Std.err"])), p_value = unname(gc["Pr(>|W|)"]),
  windows = nrow(z), patients = length(unique(z$subject_id)), events = sum(z$iap12_next),
  ci_method = "patient-cluster sandwich"
)

ord_formula <- as.formula(paste("wsacs_grade ~", fixed_rhs))
po <- ordinal::clm(ord_formula, data = z, link = "logit", Hess = TRUE)
po_v <- sandwich::vcovCL(po, cluster = z$subject_id, type = "HC0")
b <- unname(coef(po)["estimated_fluid_balance_l"]); se <- unname(sqrt(diag(po_v))["estimated_fluid_balance_l"])
rows <- rbind(rows, data.frame(
  analysis = "WSACS grade 0-IV proportional odds", effect_measure = "common OR per 1,000 mL",
  estimate = exp(b), ci_lower = exp(b - 1.96 * se), ci_upper = exp(b + 1.96 * se),
  p_value = 2 * pnorm(-abs(b / se)), windows = nrow(z), patients = length(unique(z$subject_id)), events = NA_integer_,
  ci_method = "patient-cluster sandwich"
))

partial <- ordinal::clm(ord_formula, nominal = ~ apache_ii_z + map_z, data = z, link = "logit", Hess = TRUE)
pb <- unname(coef(partial)["estimated_fluid_balance_l"])
set.seed(cfg$seed); ids <- unique(z$subject_id); partial_boot <- rep(NA_real_, cfg$bootstrap_reps)
for (i in seq_len(cfg$bootstrap_reps)) {
  sampled <- sample(ids, length(ids), replace = TRUE)
  boot_data <- do.call(rbind, unname(lapply(sampled, function(id) z[z$subject_id == id, , drop = FALSE])))
  partial_boot[i] <- tryCatch({
    fit <- ordinal::clm(ord_formula, nominal = ~ apache_ii_z + map_z, data = boot_data, link = "logit", Hess = FALSE)
    unname(coef(fit)["estimated_fluid_balance_l"])
  }, error = function(e) NA_real_)
}
successful <- sum(is.finite(partial_boot))
if (successful < .8 * cfg$bootstrap_reps) stop("Too many failed partial proportional-odds bootstrap replicates")
pse <- sd(partial_boot, na.rm = TRUE); pci <- quantile(exp(partial_boot), c(.025, .975), na.rm = TRUE, names = FALSE)
rows <- rbind(rows, data.frame(
  analysis = "Partial proportional odds (APACHE II and MAP non-proportional)", effect_measure = "common OR per 1,000 mL",
  estimate = exp(pb), ci_lower = pci[1], ci_upper = pci[2], p_value = 2 * pnorm(-abs(pb / pse)),
  windows = nrow(z), patients = length(unique(z$subject_id)), events = NA_integer_,
  ci_method = sprintf("patient-cluster percentile bootstrap (%s/%s successful)", successful, cfg$bootstrap_reps)
))

nominal <- as.data.frame(ordinal::nominal_test(po)); nominal$term <- row.names(nominal); row.names(nominal) <- NULL
grade <- as.data.frame(table(z$wsacs_grade)); names(grade) <- c("grade", "n"); grade$percent <- 100 * grade$n / sum(grade$n)
write_csv_atomic(rows, file.path(out, "alternative_outcome_results.csv")); write_csv_atomic(nominal, file.path(out, "proportional_odds_nominal_test.csv"))
write_csv_atomic(grade, file.path(out, "wsacs_grade_distribution.csv")); saveRDS(list(gee = gee, proportional_odds = po, partial = partial), file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out, list(partial_proportional_odds_bootstrap_reps = cfg$bootstrap_reps, bootstrap_unit = "patient"))
