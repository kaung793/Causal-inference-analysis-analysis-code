file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "07_time_window_heterogeneity")
d <- prepare_primary(read_analysis_file(cfg$primary_long))
window_codes <- c(`1` = "D1_D2", `2` = "D2_D3", `3` = "D3_D5", `5` = "D5_D7")
d$window <- factor(unname(window_codes[as.character(d$day)]), levels = unname(window_codes))
d <- scale_covariates(d, c("iap_current", "apache_ii", "age", "creatinine", "map"))
required <- c("estimated_fluid_balance_l", "iap_next", "iap15_next", "iap_current_z", "apache_ii_z",
              "age_z", "sex", "etiology", "creatinine_z", "map_z", "subject_id", "window")
d <- d[complete.cases(d[required]), , drop = FALSE]
if (isTRUE(cfg$strict_cohort_checks) && (nrow(d) != 2188L || length(unique(d$subject_id)) != 771L)) {
  stop("Time-window Model 3 complete-case cohort must contain 2,188 windows from 771 patients")
}

rhs <- "estimated_fluid_balance_l * window + iap_current_z + apache_ii_z + age_z + sex + etiology + creatinine_z + map_z + (1 | subject_id)"
fits <- list(
  continuous = fit_lmm_strict(as.formula(paste("iap_next ~", rhs)), d),
  iap15 = fit_glmm_strict(as.formula(paste("iap15_next ~", rhs)), d)
)

display <- c(D1_D2 = "Day 1 to 2", D2_D3 = "Day 2 to 3", D3_D5 = "Day 3 to 5", D5_D7 = "Day 5 to 7")
extract_windows <- function(fit, outcome, binary) {
  bnames <- names(fixed_coef(fit)); main <- "estimated_fluid_balance_l"
  V <- as.matrix(vcov(fit))
  rows <- list(); audits <- list()
  for (w in levels(d$window)) {
    interaction <- paste0(main, ":window", w)
    terms <- setNames(1, main)
    if (w != levels(d$window)[1L]) terms <- c(terms, setNames(1, interaction))
    e <- linear_contrast(fit, terms, exponentiate = binary)
    interaction_result <- if (w == levels(d$window)[1L]) c(estimate = ifelse(binary, 1, 0), p_value = NA) else {
      x <- extract_effect(fit, interaction, exponentiate = binary); x[c("estimate", "p_value")]
    }
    z <- d[d$window == w, , drop = FALSE]
    rows[[w]] <- data.frame(
      outcome = outcome, window = w, window_label = display[w],
      effect_measure = ifelse(binary, "OR per 1,000 mL", "beta mmHg per 1,000 mL"),
      estimate = e["estimate"], ci_lower = e["ci_lower"], ci_upper = e["ci_upper"], p_value = e["p_value"],
      contrast_vs_day1 = interaction_result["estimate"], p_interaction_vs_day1 = interaction_result["p_value"],
      window_n = nrow(z), events = if (binary) sum(z$iap15_next) else NA_integer_
    )
    var_legacy <- V[main, main]
    if (w != levels(d$window)[1L]) var_legacy <- var_legacy + V[interaction, interaction]
    se_legacy <- sqrt(var_legacy)
    linear <- e["linear_estimate"]
    audits[[w]] <- data.frame(
      outcome = outcome, window = w, linear_estimate = linear,
      corrected_se_full_vcov = e["standard_error"], legacy_se_ignoring_covariance = se_legacy,
      corrected_lower_linear = linear - 1.96 * e["standard_error"], corrected_upper_linear = linear + 1.96 * e["standard_error"],
      legacy_lower_linear = linear - 1.96 * se_legacy, legacy_upper_linear = linear + 1.96 * se_legacy,
      covariance_main_interaction = if (w == levels(d$window)[1L]) NA_real_ else V[main, interaction]
    )
  }
  interactions <- grep("estimated_fluid_balance_l:window", bnames, value = TRUE, fixed = TRUE)
  list(results = do.call(rbind, rows), audit = do.call(rbind, audits), global = cbind(outcome = outcome, joint_wald_test(fit, interactions)))
}

cont <- extract_windows(fits$continuous, "continuous_iap", FALSE)
bin <- extract_windows(fits$iap15, "iap15", TRUE)
results <- rbind(cont$results, bin$results); audit <- rbind(cont$audit, bin$audit); global <- rbind(cont$global, bin$global)
diag <- rbind(model_diagnostics(fits$continuous, "time_window_continuous"), model_diagnostics(fits$iap15, "time_window_iap15"))
write_csv_atomic(results, file.path(out, "preferred_converged_full_covariance_results.csv"))
write_csv_atomic(audit, file.path(out, "legacy_ci_formula_audit.csv"))
write_csv_atomic(global, file.path(out, "global_interaction_tests.csv"))
write_csv_atomic(diag, file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out, list(analysis_sample = "Model 3 complete case", windows = nrow(d), patients = length(unique(d$subject_id))))
