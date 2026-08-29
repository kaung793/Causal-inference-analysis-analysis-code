file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_manuscript_release) || cli$force,
                        "manuscript-release verification")) quit(save = "no", status = 0)

contract_dir <- file.path(root, "results", "manuscript_release_v2")
out <- analysis_output_dir(cfg, "manuscript_release_verification")
checks <- list(); check_index <- 0L

append_checks <- function(component, actual, expected, keys, metrics, tolerances,
                          verification_class) {
  key_actual <- do.call(paste, c(actual[keys], sep = "::"))
  key_expected <- do.call(paste, c(expected[keys], sep = "::"))
  if (anyDuplicated(key_actual) || anyDuplicated(key_expected)) stop(component, " has duplicate verification keys")
  m <- match(key_expected, key_actual)
  if (anyNA(m)) stop(component, " is missing expected result rows")
  actual <- actual[m, , drop = FALSE]
  for (metric in metrics) {
    tol <- tolerances[[metric]] %||% tolerances[[1L]]
    delta <- abs(as.numeric(actual[[metric]]) - as.numeric(expected[[metric]]))
    ok_na <- is.na(actual[[metric]]) & is.na(expected[[metric]])
    pass <- ok_na | (!is.na(delta) & delta <= tol)
    for (i in seq_len(nrow(expected))) {
      check_index <<- check_index + 1L
      checks[[check_index]] <<- data.frame(
        component = component, result_key = key_expected[i], metric = metric,
        expected = expected[[metric]][i], actual = actual[[metric]][i],
        absolute_difference = if (ok_na[i]) 0 else delta[i], tolerance = tol,
        verification_class = verification_class,
        status = if (pass[i]) "PASS" else "FAIL", stringsAsFactors = FALSE
      )
    }
  }
}

aipw_actual_path <- file.path(cfg$output_dir, "manuscript_release_aipw",
                              "manuscript_aipw_patient_bootstrap.csv")
rcs_actual_path <- file.path(cfg$output_dir, "04_rcs_dose_response",
                             "rcs_tests_and_exploratory_inflections.csv")
time_actual_path <- file.path(cfg$output_dir, "manuscript_release_time_window",
                              "document_display_formula_current_runtime.csv")
for (path in c(aipw_actual_path, rcs_actual_path, time_actual_path)) {
  if (!file.exists(path)) stop("Required manuscript-release output is missing: ", path)
}

aipw_actual <- read.csv(aipw_actual_path, check.names = FALSE)
aipw_expected <- read.csv(file.path(contract_dir, "aipw.csv"), check.names = FALSE)
append_checks(
  "AIPW", aipw_actual, aipw_expected, c("scenario", "estimand"),
  c("estimate", "ci_lower", "ci_upper"),
  list(estimate = 1e-8, ci_lower = 1e-8, ci_upper = 1e-8),
  "exact deterministic rerun from dedicated freeze"
)

rcs_actual <- read.csv(rcs_actual_path, check.names = FALSE)
names(rcs_actual)[names(rcs_actual) == "exploratory_inflection_ml"] <- "inflection_ml"
names(rcs_actual)[names(rcs_actual) == "exploratory_plateau_ml"] <- "plateau_ml"
rcs_expected <- read.csv(file.path(contract_dir, "rcs.csv"), check.names = FALSE)
append_checks(
  "RCS", rcs_actual, rcs_expected, c("exposure", "outcome"),
  c("overall_p", "nonlinearity_p", "inflection_ml", "plateau_ml"),
  list(overall_p = 1e-8, nonlinearity_p = 1e-8,
       inflection_ml = 1e-6, plateau_ml = 1e-6),
  "exact rerun from dedicated RCS freeze"
)

time_actual <- read.csv(time_actual_path, check.names = FALSE)
time_expected <- read.csv(file.path(contract_dir, "time_window.csv"), check.names = FALSE)
append_checks(
  "time-window", time_actual, time_expected, c("outcome", "window"),
  c("estimate", "ci_lower", "ci_upper", "p_interaction_vs_day1"),
  list(estimate = .01, ci_lower = .02, ci_upper = .02, p_interaction_vs_day1 = .015),
  "runtime concordance to archived document output; package-sensitive binary fit"
)

summary <- do.call(rbind, checks)
write_csv_atomic(summary, file.path(out, "verification_summary.csv"))
write_run_metadata(cfg, c(aipw_actual_path, rcs_actual_path, time_actual_path), out,
                   list(checks = nrow(summary), failures = sum(summary$status == "FAIL")))
if (any(summary$status == "FAIL")) {
  failed <- summary[summary$status == "FAIL", ]
  stop("Manuscript-release verification failed for ", nrow(failed), " metric(s)")
}
message("PASS: manuscript-release AIPW and RCS reproduced exactly; time-window runtime output remained within the documented package-sensitive tolerance")
