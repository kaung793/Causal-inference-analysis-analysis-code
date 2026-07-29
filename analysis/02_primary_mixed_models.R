file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "02_primary_mixed_models")
d <- prepare_primary(read_analysis_file(cfg$primary_long))
d <- scale_covariates(d, c("iap_current", "apache_ii", "age", "creatinine", "map"))

specs <- expand.grid(
  exposure = c("fluid_intake_l", "fluid_balance_l"),
  outcome = c("iap_next", "iap15_next", "iap20_next"),
  model = 1:3, stringsAsFactors = FALSE
)
rows <- list(); diagnostics <- list(); fits <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  covars <- if (s$model == 1L) character() else c("iap_current_z", "apache_ii_z", "age_z", "sex", "etiology")
  if (s$model == 3L) covars <- c(covars, "creatinine_z", "map_z")
  rhs <- paste(c(s$exposure, covars, "(1 | subject_id)"), collapse = " + ")
  form <- as.formula(paste(s$outcome, "~", rhs))
  binary <- s$outcome != "iap_next"
  fit <- if (binary) fit_glmm_strict(form, d) else fit_lmm_strict(form, d)
  label <- paste(s$exposure, s$outcome, paste0("Model", s$model), sep = "__")
  fits[[label]] <- fit
  effect <- extract_effect(fit, s$exposure, exponentiate = binary)
  used <- model.frame(fit)
  rows[[i]] <- data.frame(
    exposure = ifelse(s$exposure == "fluid_intake_l", "Fluid intake", "Fluid balance"),
    outcome = s$outcome, model = paste0("Model ", s$model), effect_measure = ifelse(binary, "OR per 1,000 mL", "beta mmHg per 1,000 mL"),
    estimate = effect["estimate"], ci_lower = effect["ci_lower"], ci_upper = effect["ci_upper"], p_value = effect["p_value"],
    windows = nrow(used), patients = length(unique(as.character(used$subject_id))),
    stringsAsFactors = FALSE
  )
  diagnostics[[i]] <- model_diagnostics(fit, label)
}
results <- do.call(rbind, rows); row.names(results) <- NULL
diag <- do.call(rbind, diagnostics); row.names(diag) <- NULL
write_csv_atomic(results, file.path(out, "primary_mixed_model_results.csv"))
write_csv_atomic(diag, file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out)
message("Primary models complete: ", nrow(results), " estimates")
