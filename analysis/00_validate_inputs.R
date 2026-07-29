file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript")
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args()
cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "00_input_validation")

raw <- read_analysis_file(cfg$primary_long)
dat <- prepare_primary(raw)
strict <- isTRUE(cfg$strict_cohort_checks)
panel <- validate_panel(
  dat, cfg$expected_primary_patients, cfg$expected_primary_windows,
  expected_days = c(1, 2, 3, 5), strict = strict
)

never_impute <- c("fluid_intake_ml", "fluid_balance_ml", "iap_current", "iap_next", "age", "sex", "etiology")
missing_never <- vapply(dat[never_impute], function(x) sum(is.na(x)), integer(1))
if (strict && any(missing_never > 0L)) {
  stop("Observed exposure, outcome, or baseline fields contain missing values: ",
       paste(names(missing_never)[missing_never > 0L], collapse = ", "))
}

lag_checks <- data.frame(comparison = c("D1 next vs D2 current", "D2 next vs D3 current", "D3 next vs D5 current"),
                         compared = 0L, max_absolute_difference = NA_real_)
wide <- reshape(dat[c("subject_id", "day", "iap_current", "iap_next")],
                idvar = "subject_id", timevar = "day", direction = "wide")
pairs <- list(c("iap_next.1", "iap_current.2"), c("iap_next.2", "iap_current.3"), c("iap_next.3", "iap_current.5"))
for (i in seq_along(pairs)) {
  a <- wide[[pairs[[i]][1L]]]; b <- wide[[pairs[[i]][2L]]]
  keep <- is.finite(a) & is.finite(b)
  lag_checks$compared[i] <- sum(keep)
  lag_checks$max_absolute_difference[i] <- if (any(keep)) max(abs(a[keep] - b[keep])) else NA_real_
}
if (strict && any(lag_checks$max_absolute_difference > 1e-8, na.rm = TRUE)) stop("Lagged IAP sequence is inconsistent")

balance_check <- data.frame(compared = 0L, median_absolute_difference_ml = NA_real_, max_absolute_difference_ml = NA_real_)
if ("urine_output_ml" %in% names(dat)) {
  delta <- dat$fluid_balance_ml - (dat$fluid_intake_ml - dat$urine_output_ml)
  keep <- is.finite(delta)
  balance_check <- data.frame(compared = sum(keep), median_absolute_difference_ml = median(abs(delta[keep])),
                              max_absolute_difference_ml = max(abs(delta[keep])))
}

missingness <- data.frame(
  variable = names(dat), missing_n = vapply(dat, function(x) sum(is.na(x)), integer(1)),
  missing_percent = vapply(dat, function(x) 100 * mean(is.na(x)), numeric(1))
)
write_csv_atomic(panel, file.path(out, "panel_dimensions.csv"))
write_csv_atomic(lag_checks, file.path(out, "lag_sequence_checks.csv"))
write_csv_atomic(balance_check, file.path(out, "fluid_balance_identity_check.csv"))
write_csv_atomic(missingness, file.path(out, "missingness.csv"))
write_run_metadata(cfg, cfg$primary_long, out)
message("Input validation passed: ", panel$patients, " patients / ", panel$windows, " windows")
