file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_revision_sensitivities) || cli$force, "additional-adjustment analyses")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "11_additional_adjustments")
revision_source <- cfg$revision_primary_long %||% cfg$primary_long
if (!file.exists(revision_source)) revision_source <- cfg$primary_long

comorbidity_fields <- c("smoking", "drinking", "hypertension", "diabetes", "hyperlipidemia_history", "copd")
base_data <- prepare_primary(read_analysis_file(revision_source))
base_data <- merge_optional_fields(base_data, cfg$patient_flags, comorbidity_fields, "patient flags")
et <- suppressWarnings(as.numeric(as.character(base_data$etiology)))
base_data$etiology <- factor(ifelse(et %in% c(4, 5), 4, et))

revision_fields <- c("fluid_balance_ml", "iap_current", "iap_next", "apache_ii", "creatinine", "map",
                     "age", "weight_kg", "heart_rate", "albumin", "ph", "spo2")
corrected_data <- merge_optional_fields(base_data, cfg$revision_covariates, revision_fields, "revision covariates")
corrected_data$fluid_balance_l <- corrected_data$fluid_balance_ml / 1000
corrected_data$iap12_next <- as.integer(corrected_data$iap_next >= 12)
corrected_data$iap15_next <- as.integer(corrected_data$iap_next >= 15)
corrected_data$iap20_next <- as.integer(corrected_data$iap_next >= 20)

base_vars <- c("iap_current", "apache_ii", "age", "creatinine", "map")
optional_numeric <- c("weight_kg", "heart_rate", "albumin", "ph", "spo2")
base_data <- scale_covariates(base_data, c(base_vars, intersect(optional_numeric, names(base_data))))
corrected_data <- scale_covariates(corrected_data, c(base_vars, intersect(optional_numeric, names(corrected_data))))
datasets <- list(base = base_data, corrected = corrected_data)
model3_covars <- c("iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z")
scenarios <- list(Model3 = list(covars = model3_covars, estimator = "mixed", dataset = "base"))
if (all(c("weight_kg_z", "heart_rate_z", "albumin_z") %in% names(corrected_data))) {
  scenarios$Model4 <- list(covars = c(model3_covars, "weight_kg_z", "heart_rate_z", "albumin_z"), estimator = "mixed", dataset = "corrected")
  if ("ph_z" %in% names(corrected_data)) scenarios$Model4_plus_pH <- list(covars = c(scenarios$Model4$covars, "ph_z"), estimator = "mixed", dataset = "corrected")
}
if ("spo2_z" %in% names(base_data)) {
  scenarios$Model3_plus_SpO2 <- list(covars = c(model3_covars, "spo2_z"), estimator = "mixed", dataset = "base")
  scenarios$Model3_plus_SpO2_GEE <- list(covars = c(model3_covars, "spo2_z"), estimator = "GEE with patient-clustered sandwich variance", dataset = "base")
}
if (all(comorbidity_fields %in% names(base_data))) {
  for (v in comorbidity_fields) datasets$base[[v]] <- factor(datasets$base[[v]])
  scenarios$Model3_plus_comorbidities <- list(covars = c(model3_covars, comorbidity_fields), estimator = "mixed", dataset = "base", allow_warnings = TRUE)
}

rows <- list(); diagnostics <- list(); fits <- list(); idx <- 0L
for (scenario in names(scenarios)) for (outcome in c("iap_next", "iap15_next", "iap20_next")) {
  spec <- scenarios[[scenario]]; data <- datasets[[spec$dataset]]; binary <- outcome != "iap_next"
  need <- c(outcome, "fluid_balance_l", "subject_id", spec$covars)
  z <- droplevels(data[complete.cases(data[need]), , drop = FALSE])
  fixed_formula <- as.formula(paste(outcome, "~ fluid_balance_l +", paste(spec$covars, collapse = " + ")))
  if (spec$estimator == "mixed") {
    mixed_formula <- update(fixed_formula, . ~ . + (1 | subject_id))
    fit <- if (binary && isTRUE(spec$allow_warnings)) fit_glmm_audited(mixed_formula, z) else if (binary) fit_glmm_strict(mixed_formula, z) else fit_lmm_strict(mixed_formula, z)
    e <- extract_effect(fit, "fluid_balance_l", exponentiate = binary)
    estimate <- unname(e["estimate"]); lower <- unname(e["ci_lower"]); upper <- unname(e["ci_upper"]); p <- unname(e["p_value"])
  } else {
    require_packages("geepack")
    fit <- geepack::geeglm(fixed_formula, id = subject_id, data = z, family = if (binary) binomial() else gaussian(),
                           corstr = "independence", std.err = "san.se")
    co <- summary(fit)$coefficients
    b <- as.numeric(co["fluid_balance_l", "Estimate"]); se <- as.numeric(co["fluid_balance_l", "Std.err"])
    p <- as.numeric(co["fluid_balance_l", "Pr(>|W|)"])
    estimate <- if (binary) exp(b) else b; lower <- if (binary) exp(b - 1.96 * se) else b - 1.96 * se
    upper <- if (binary) exp(b + 1.96 * se) else b + 1.96 * se
  }
  estimator_label <- if (isTRUE(spec$allow_warnings)) "mixed; audited convergence warnings retained" else spec$estimator
  idx <- idx + 1L
  rows[[idx]] <- data.frame(scenario = scenario, estimator = estimator_label, data_version = spec$dataset, outcome = outcome,
                            effect_measure = ifelse(binary, "OR per 1,000 mL", "beta per 1,000 mL"),
                            estimate = estimate, ci_lower = lower, ci_upper = upper, p_value = p,
                            windows = nrow(z), patients = length(unique(z$subject_id)))
  label <- paste(scenario, outcome, sep = "__"); fits[[label]] <- fit
  diagnostics[[idx]] <- transform(model_diagnostics(fit, label), estimator = estimator_label, data_version = spec$dataset)
}
write_csv_atomic(do.call(rbind, rows), file.path(out, "additional_adjustment_results.csv"))
write_csv_atomic(do.call(rbind, diagnostics), file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, c(revision_source, cfg$revision_covariates, cfg$patient_flags), out)
