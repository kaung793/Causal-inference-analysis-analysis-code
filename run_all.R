file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else normalizePath(getwd(), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
scripts <- c(
  "00_validate_inputs.R", "01_descriptive_missingness.R", "02_primary_mixed_models.R",
  "03_external_validation.R", "04_rcs_dose_response.R", "05_msm_iptw.R", "06_aipw.R",
  "07_time_window_heterogeneity.R", "08_equal_lag_sensitivity.R", "09_mediation_mi_rubin.R",
  "10_alternative_iap_outcomes.R", "11_additional_adjustments.R", "12_exclusions_subgroups.R",
  "13_standardized_absolute_risks.R", "14_model3_multiple_imputation.R", "15_repeated_threshold.R", "16_figures.R"
)
if (!is.null(cli$only)) {
  requested <- trimws(strsplit(cli$only, ",", fixed = TRUE)[[1L]])
  scripts <- scripts[basename(tools::file_path_sans_ext(scripts)) %in% requested | scripts %in% requested]
  if (!length(scripts)) stop("--only did not match an analysis script")
}
rscript <- file.path(R.home("bin"), "Rscript")
config_arg <- paste0("--config=", cfg$config_file)
ensure_dir(file.path(root, "logs"))
status <- data.frame(script = scripts, exit_status = NA_integer_)
for (i in seq_along(scripts)) {
  path <- file.path(root, "analysis", scripts[i])
  message("RUN: ", scripts[i])
  code <- system2(rscript, args = c(shQuote(path), shQuote(config_arg)))
  status$exit_status[i] <- code
  if (!identical(code, 0L)) {
    write.csv(status, file.path(root, "logs", "run_all_status.csv"), row.names = FALSE)
    stop("Analysis failed: ", scripts[i], " (exit status ", code, ")")
  }
}
ensure_dir(file.path(root, "logs")); write.csv(status, file.path(root, "logs", "run_all_status.csv"), row.names = FALSE)
message("All selected analyses completed")
