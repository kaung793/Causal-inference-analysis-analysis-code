file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R")); source(file.path(root, "R", "mi_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_model3_mi) || cli$force, "Model 3 multiple imputation")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "14_model3_multiple_imputation")
input_override <- Sys.getenv("MODEL3_MI_LONG", unset = "")
d_input <- if (nzchar(input_override)) normalizePath(input_override, winslash = "/", mustWork = FALSE) else
  cfg$model3_mi_long %||% cfg$primary_long
d <- prepare_primary(read_analysis_file(d_input))
validate_panel(d, cfg$expected_primary_patients, cfg$expected_primary_windows, c(1, 2, 3, 5), isTRUE(cfg$strict_cohort_checks))

days <- c(1L, 2L, 3L, 5L)
if (!all(table(d$subject_id) == length(days))) stop("Model 3 MI requires four rows per patient.")
assert_columns(d, c("sbp", "dbp"), "Model 3 MI source")
non_imputed <- c("fluid_intake_ml", "estimated_fluid_balance_ml", "iap_current", "iap_next",
                 "age", "sex", "etiology")
if (anyNA(d[non_imputed])) stop("Exposure, IAP, or complete baseline predictors unexpectedly contain missing values.")

# The adjudicated source freeze contains one pre-identified implausible MAP
# measurement at Day 3 that is deterministically reconstructed from observed
# SBP/DBP. Selecting it by data state and day avoids publishing a subject ID.
reconstruction_day <- as.integer(cfg$model3_mi_map_reconstruction_day %||% 3L)
map_reconstruction <- is.na(d$map) & d$day == reconstruction_day &
  stats::complete.cases(d[c("sbp", "dbp")])
expected_reconstructions <- as.integer(cfg$model3_mi_map_reconstruction_n %||% 1L)
if (isTRUE(cfg$strict_cohort_checks) && sum(map_reconstruction) != expected_reconstructions) {
  stop("Unexpected number of pre-identified MAP reconstruction rows: ", sum(map_reconstruction))
}
d$map[map_reconstruction] <- (d$sbp[map_reconstruction] + 2 * d$dbp[map_reconstruction]) / 3

id_order <- unique(d$subject_id)
key <- paste(d$subject_id, d$day, sep = "__")
row_for <- function(day_value) {
  idx <- match(paste(id_order, day_value, sep = "__"), key)
  if (anyNA(idx)) stop("Incomplete patient-day panel while constructing Model 3 MI data.")
  idx
}
at_day <- function(variable, day_value) d[[variable]][row_for(day_value)]

if (!isTRUE(all.equal(at_day("iap_next", 1L), at_day("iap_current", 2L), tolerance = 0)) ||
    !isTRUE(all.equal(at_day("iap_next", 2L), at_day("iap_current", 3L), tolerance = 0)) ||
    !isTRUE(all.equal(at_day("iap_next", 3L), at_day("iap_current", 5L), tolerance = 0))) {
  stop("Adjacent lagged IAP values are inconsistent.")
}

baseline_idx <- match(id_order, d$subject_id)
etiology_value <- as.character(d$etiology[baseline_idx])
wide <- data.frame(
  analysis_id = seq_along(id_order),
  age = d$age[baseline_idx],
  sex = factor(as.character(d$sex[baseline_idx])),
  etiology4 = factor(ifelse(etiology_value %in% c("4", "5"), "4", etiology_value)),
  stringsAsFactors = FALSE
)
wide_specs <- c(
  fluid_t = "fluid_intake_ml",
  fluid_balance_t = "estimated_fluid_balance_ml",
  apache_t = "apache_ii",
  cr_t = "creatinine",
  map_t = "map"
)
for (prefix in names(wide_specs)) {
  for (day_value in days) {
    wide[[paste0(prefix, "_d", day_value)]] <- at_day(wide_specs[[prefix]], day_value)
  }
}
wide$iap_d1 <- at_day("iap_current", 1L)
wide$iap_d2 <- at_day("iap_next", 1L)
wide$iap_d3 <- at_day("iap_next", 2L)
wide$iap_d5 <- at_day("iap_next", 3L)
wide$iap_d7 <- at_day("iap_next", 5L)

targets <- unlist(lapply(c("apache_t", "cr_t", "map_t"), function(v) paste0(v, "_d", days)))
if (any(vapply(wide[setdiff(names(wide), targets)], anyNA, logical(1)))) {
  stop("A non-target imputation predictor contains missing values.")
}
seed_override <- Sys.getenv("MODEL3_MI_SEED", unset = "")
mi_seed <- if (nzchar(seed_override)) as.integer(seed_override) else
  as.integer(cfg$model3_mi_seed %||% cfg$seed)
imp <- make_patient_wide_imputation(
  wide, targets, cfg$mi_m, cfg$mi_maxit, mi_seed, cfg$mi_donors,
  id_column = "analysis_id"
)
saveRDS(imp, file.path(out, "model3_patient_wide_mice.rds"), compress = "xz")
write_csv_atomic(wide, file.path(out, "analysis_data_wide_preimputation.csv"))

missingness <- do.call(rbind, lapply(days, function(day_value) {
  z <- d[d$day == day_value, , drop = FALSE]
  data.frame(
    day = day_value, windows = nrow(z),
    apache_missing_n = sum(is.na(z$apache_ii)),
    creatinine_missing_n = sum(is.na(z$creatinine)),
    map_missing_n = sum(is.na(z$map)),
    model3_complete_n = sum(stats::complete.cases(z[c("apache_ii", "creatinine", "map")]))
  )
}))
write_csv_atomic(missingness, file.path(out, "missingness_by_day.csv"))
logged <- imp$loggedEvents
if (is.null(logged) || !nrow(logged)) logged <- data.frame(message = "No logged MICE events")
write_csv_atomic(logged, file.path(out, "mice_logged_events.csv"))

to_long <- function(w) {
  next_day <- c(`1` = 2L, `2` = 3L, `3` = 5L, `5` = 7L)
  parts <- lapply(days, function(day_value) {
    data.frame(
      analysis_id = factor(w$analysis_id),
      day = day_value,
      fluid_intake_l = w[[paste0("fluid_t_d", day_value)]] / 1000,
      estimated_fluid_balance_l = w[[paste0("fluid_balance_t_d", day_value)]] / 1000,
      iap_current = w[[paste0("iap_d", day_value)]],
      iap_next = w[[paste0("iap_d", next_day[as.character(day_value)])]],
      apache_ii = w[[paste0("apache_t_d", day_value)]],
      creatinine = w[[paste0("cr_t_d", day_value)]],
      map = w[[paste0("map_t_d", day_value)]],
      age = w$age,
      sex = factor(w$sex),
      etiology4 = factor(w$etiology4)
    )
  })
  z <- do.call(rbind, parts)
  z <- z[order(as.integer(z$analysis_id), z$day), , drop = FALSE]
  z$iap15_next <- as.integer(z$iap_next >= 15); z$iap20_next <- as.integer(z$iap_next >= 20)
  z$iap5 <- z$iap_current / 5; z$apache10 <- z$apache_ii / 10
  z$age10 <- z$age / 10; z$creatinine100 <- z$creatinine / 100; z$map10 <- z$map / 10
  if (nrow(z) != cfg$expected_primary_windows || anyNA(z)) {
    stop("Completed Model 3 MI long data have unexpected dimensions or missing values.")
  }
  z
}

specs <- expand.grid(exposure = c("fluid_intake_l", "estimated_fluid_balance_l"),
                     outcome = c("iap_next", "iap15_next", "iap20_next"), stringsAsFactors = FALSE)
per <- list(); diag <- list(); idx <- 0L
for (i in seq_len(cfg$mi_m)) {
  z <- to_long(mice::complete(imp, i))
  for (j in seq_len(nrow(specs))) {
    s <- specs[j, ]; binary <- s$outcome != "iap_next"
    rhs <- paste(c(s$exposure, "iap5", "apache10", "age10", "sex", "etiology4",
                   "creatinine100", "map10", "(1 | analysis_id)"), collapse = " + ")
    fit <- if (binary) fit_glmm_strict(as.formula(paste(s$outcome, "~", rhs)), z) else fit_lmm_strict(as.formula(paste(s$outcome, "~", rhs)), z)
    e <- extract_effect(fit, s$exposure, exponentiate = FALSE); idx <- idx + 1L
    per[[idx]] <- data.frame(imputation = i, exposure = s$exposure, outcome = s$outcome,
                             estimate = e["linear_estimate"], variance = e["standard_error"]^2)
    diag[[idx]] <- transform(model_diagnostics(fit, paste(s$exposure, s$outcome, sep = "__")), imputation = i)
  }
}
per <- do.call(rbind, per); pooled <- rubin_pool_grouped(per, c("exposure", "outcome"), "estimate", "variance")
pooled$binary <- pooled$outcome != "iap_next"; pooled$effect_measure <- ifelse(pooled$binary, "OR per 1,000 mL", "beta per 1,000 mL")
pooled$estimate_display <- ifelse(pooled$binary, exp(pooled$estimate), pooled$estimate)
pooled$ci_lower_display <- ifelse(pooled$binary, exp(pooled$ci_lower), pooled$ci_lower)
pooled$ci_upper_display <- ifelse(pooled$binary, exp(pooled$ci_upper), pooled$ci_upper)
pooled$patients <- cfg$expected_primary_patients; pooled$windows <- cfg$expected_primary_windows
write_csv_atomic(per, file.path(out, "estimates_by_imputation.csv")); write_csv_atomic(pooled, file.path(out, "model3_mi_rubin_pooled.csv"))
write_csv_atomic(do.call(rbind, diag), file.path(out, "model_diagnostics.csv"))
write_run_metadata(
  cfg, d_input, out,
  list(m = cfg$mi_m, maxit = cfg$mi_maxit, seed = mi_seed,
       map_reconstructed_n = sum(map_reconstruction))
)