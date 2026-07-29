file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_revision_sensitivities) || cli$force, "repeated-threshold analysis")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "15_repeated_threshold")
d <- prepare_primary(read_analysis_file(cfg$primary_long))
next_rows <- d[c("subject_id", "day", "iap_next")]; names(next_rows) <- c("subject_id", "next_start_day", "second_iap")
z <- d[d$day %in% c(1, 2, 3), , drop = FALSE]
z$next_start_day <- c(`1` = 2, `2` = 3, `3` = 5)[as.character(z$day)]
z$first_iap <- z$iap_next
z <- merge(z, next_rows, by = c("subject_id", "next_start_day"), all.x = TRUE, sort = FALSE)
z$repeated_iap12 <- as.integer(z$first_iap >= 12 & z$second_iap >= 12)
z$repeated_iap15 <- as.integer(z$first_iap >= 15 & z$second_iap >= 15)
z <- scale_covariates(z, c("iap_current", "apache_ii", "age", "creatinine", "map"))

rows <- list(); fits <- list(); diagnostics <- list()
for (outcome in c("repeated_iap12", "repeated_iap15")) {
  need <- c(outcome, "fluid_balance_l", "iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z", "subject_id")
  a <- droplevels(z[complete.cases(z[need]), , drop = FALSE])
  f <- as.formula(paste(outcome, "~ fluid_balance_l + iap_current_z + apache_ii_z + age_z + sex + etiology + creatinine_z + map_z + (1 | subject_id)"))
  fit <- fit_glmm_strict(f, a); e <- extract_effect(fit, "fluid_balance_l", exponentiate = TRUE)
  rows[[outcome]] <- data.frame(outcome = outcome, effect_measure = "OR per 1,000 mL", estimate = e["estimate"],
    ci_lower = e["ci_lower"], ci_upper = e["ci_upper"], p_value = e["p_value"], windows = nrow(a),
    patients = length(unique(a$subject_id)), events = sum(a[[outcome]]), event_percent = 100 * mean(a[[outcome]]))
  fits[[outcome]] <- fit; diagnostics[[outcome]] <- model_diagnostics(fit, outcome)
}
write_csv_atomic(do.call(rbind, rows), file.path(out, "repeated_threshold_results.csv")); write_csv_atomic(do.call(rbind, diagnostics), file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz"); write_run_metadata(cfg, cfg$primary_long, out)
