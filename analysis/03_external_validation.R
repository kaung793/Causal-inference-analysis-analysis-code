file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_external_validation) || cli$force, "external validation")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "03_external_validation")
d <- prepare_external(read_analysis_file(cfg$external_long))
validate_panel(d, cfg$expected_external_patients, cfg$expected_external_windows,
               expected_days = c(1, 2, 3, 5), strict = isTRUE(cfg$strict_cohort_checks))
d <- scale_covariates(d, c("iap_current", "age", "map"))

forms <- list(
  iap_next = iap_next ~ fluid_intake_l + iap_current_z + age_z + sex + etiology + map_z + (1 | subject_id),
  iap15_next = iap15_next ~ fluid_intake_l + iap_current_z + age_z + sex + etiology + map_z + (1 | subject_id)
)
fits <- list(
  iap_next = fit_lmm_strict(forms$iap_next, d),
  iap15_next = fit_glmm_strict(forms$iap15_next, d)
)
rows <- do.call(rbind, lapply(names(fits), function(outcome) {
  binary <- outcome != "iap_next"; e <- extract_effect(fits[[outcome]], "fluid_intake_l", exponentiate = binary)
  mf <- model.frame(fits[[outcome]])
  data.frame(outcome = outcome, effect_measure = ifelse(binary, "OR per 1,000 mL", "beta mmHg per 1,000 mL"),
             estimate = e["estimate"], ci_lower = e["ci_lower"], ci_upper = e["ci_upper"], p_value = e["p_value"],
             windows = nrow(mf), patients = length(unique(mf$subject_id)))
}))
diag <- do.call(rbind, Map(model_diagnostics, fits, paste0("external_", names(fits))))
write_csv_atomic(rows, file.path(out, "external_validation_results.csv"))
write_csv_atomic(diag, file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$external_long, out)
