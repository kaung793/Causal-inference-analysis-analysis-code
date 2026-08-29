file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else
  normalizePath(getwd(), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (!isTRUE(cfg$run_manuscript_release) && !cli$force) {
  stop("Set run_manuscript_release=TRUE in the local configuration, or use --force")
}

scripts <- c(
  "analysis/00_validate_inputs.R",
  "analysis/02_primary_mixed_models.R",
  "analysis/03_external_validation.R",
  "analysis/04_rcs_dose_response.R",
  "analysis/05_msm_iptw.R",
  "analysis/07_time_window_heterogeneity.R",
  "analysis/08_equal_lag_sensitivity.R",
  "analysis/09_mediation_mi_rubin.R",
  "analysis/14_model3_multiple_imputation.R",
  "analysis/15_repeated_threshold.R",
  "release/01_manuscript_aipw.R",
  "release/02_manuscript_time_window.R",
  "release/03_check_manuscript_release.R"
)
if (!is.null(cli$only)) {
  requested <- trimws(strsplit(cli$only, ",", fixed = TRUE)[[1L]])
  id <- tools::file_path_sans_ext(basename(scripts))
  scripts <- scripts[id %in% requested | basename(scripts) %in% requested | scripts %in% requested]
  if (!length(scripts)) stop("--only did not match a manuscript-release script")
}

rscript <- file.path(R.home("bin"), "Rscript")
config_arg <- paste0("--config=", cfg$config_file)
ensure_dir(file.path(root, "logs"))
status <- data.frame(script = scripts, exit_status = NA_integer_)
for (i in seq_along(scripts)) {
  path <- file.path(root, scripts[i])
  message("RUN: ", scripts[i])
  code <- system2(rscript, args = c(shQuote(path), shQuote(config_arg), "--force"))
  status$exit_status[i] <- code
  if (!identical(code, 0L)) {
    write.csv(status, file.path(root, "logs", "run_manuscript_release_status.csv"), row.names = FALSE)
    stop("Manuscript-release analysis failed: ", scripts[i], " (exit status ", code, ")")
  }
}
write.csv(status, file.path(root, "logs", "run_manuscript_release_status.csv"), row.names = FALSE)
message("All enabled manuscript-release analyses completed")
