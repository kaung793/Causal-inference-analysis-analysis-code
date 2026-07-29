file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "01_descriptive_missingness")
d <- prepare_primary(read_analysis_file(cfg$primary_long))

vars <- c("fluid_intake_ml", "fluid_balance_ml", "iap_current", "iap_next", "apache_ii", "creatinine", "map")
by_day <- do.call(rbind, lapply(sort(unique(d$day)), function(day_value) {
  z <- d[d$day == day_value, , drop = FALSE]
  do.call(rbind, lapply(vars, function(v) data.frame(
    day = day_value, variable = v, n = nrow(z), observed_n = sum(!is.na(z[[v]])),
    missing_n = sum(is.na(z[[v]])), missing_percent = 100 * mean(is.na(z[[v]])),
    mean = mean(z[[v]], na.rm = TRUE), sd = sd(z[[v]], na.rm = TRUE),
    median = median(z[[v]], na.rm = TRUE), q1 = unname(quantile(z[[v]], .25, na.rm = TRUE)),
    q3 = unname(quantile(z[[v]], .75, na.rm = TRUE))
  )))
}))

patient_level <- d[!duplicated(d$subject_id), c("subject_id", "age", "sex", "etiology")]
cohort <- data.frame(
  patients = nrow(patient_level), windows = nrow(d),
  age_mean = mean(patient_level$age), age_sd = sd(patient_level$age),
  female_n = sum(as.character(patient_level$sex) %in% c("0", "F", "Female")),
  iap15_windows = sum(d$iap15_next), iap20_windows = sum(d$iap20_next)
)
write_csv_atomic(by_day, file.path(out, "descriptive_by_day.csv"))
write_csv_atomic(cohort, file.path(out, "cohort_summary.csv"))
write_run_metadata(cfg, cfg$primary_long, out)
