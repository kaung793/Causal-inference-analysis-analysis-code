file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_revision_sensitivities) || cli$force, "exclusion/subgroup analyses")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "12_exclusions_subgroups")
d <- prepare_primary(read_analysis_file(cfg$primary_long))
d <- merge_optional_fields(d, cfg$patient_flags, c("pcd", "crrt", "shock"), "patient flags")

et <- suppressWarnings(as.numeric(as.character(d$etiology))); d$etiology <- factor(ifelse(et %in% c(4, 5), 4, et))
d <- scale_covariates(d, c("iap_current", "apache_ii", "age", "creatinine", "map"))
base_covars <- c("iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z")

fit_set <- function(data, scenario, exposure = "estimated_fluid_balance_l", multiplier = 1) {
  rows <- list()
  for (outcome in c("iap_next", "iap15_next", "iap20_next")) {
    binary <- outcome != "iap_next"; need <- c(outcome, exposure, "subject_id", base_covars)
    z <- droplevels(data[complete.cases(data[need]), , drop = FALSE])
    f <- as.formula(paste(outcome, "~", exposure, "+", paste(base_covars, collapse = " + "), "+ (1 | subject_id)"))
    fit <- if (binary) fit_glmm_strict(f, z) else fit_lmm_strict(f, z)
    e <- extract_effect(fit, exposure, multiplier = multiplier, exponentiate = binary)
    rows[[outcome]] <- data.frame(scenario = scenario, outcome = outcome, effect_measure = ifelse(binary, "OR", "beta"),
      estimate = e["estimate"], ci_lower = e["ci_lower"], ci_upper = e["ci_upper"], p_value = e["p_value"],
      windows = nrow(z), patients = length(unique(z$subject_id)))
  }
  do.call(rbind, rows)
}

results <- list(); idx <- 0L
if ("weight_kg" %in% names(d)) {
  d$balance_per_kg_10 <- d$estimated_fluid_balance_ml / d$weight_kg / 10
  idx <- idx + 1L; results[[idx]] <- fit_set(d[is.finite(d$balance_per_kg_10) & d$weight_kg > 0, ], "weight-normalized, per 10 mL/kg", "balance_per_kg_10")
  idx <- idx + 1L; results[[idx]] <- fit_set(d[is.finite(d$balance_per_kg_10) & d$weight_kg > 0, ], "same weight-complete subset, per 1,000 mL")
}
for (flag in c("pcd", "crrt")) if (flag %in% names(d)) {
  keep <- is.na(d[[flag]]) | as.numeric(as.character(d[[flag]])) == 0
  idx <- idx + 1L; results[[idx]] <- fit_set(d[keep, ], paste0("exclude_", flag))
}
if ("shock" %in% names(d)) {
  shock_num <- as.integer(as.numeric(as.character(d$shock)) != 0)
  for (g in 0:1) { idx <- idx + 1L; results[[idx]] <- fit_set(d[shock_num == g, ], paste0("shock_", g)) }
  z <- d[complete.cases(d[c("shock", "iap_next", "iap15_next", "estimated_fluid_balance_l", "subject_id", base_covars)]), , drop = FALSE]
  z$shock_group <- factor(as.integer(as.numeric(as.character(z$shock)) != 0))
  for (outcome in c("iap_next", "iap15_next", "iap20_next")) {
    binary <- outcome != "iap_next"
    f <- as.formula(paste(outcome, "~ estimated_fluid_balance_l * shock_group +", paste(base_covars, collapse = " + "), "+ (1 | subject_id)"))
    fit <- if (binary) fit_glmm_strict(f, z) else fit_lmm_strict(f, z)
    term <- "estimated_fluid_balance_l:shock_group1"; e <- extract_effect(fit, term, exponentiate = binary)
    idx <- idx + 1L; results[[idx]] <- data.frame(scenario = "shock_interaction", outcome = outcome,
      effect_measure = ifelse(binary, "ratio of ORs", "delta beta"), estimate = e["estimate"], ci_lower = e["ci_lower"],
      ci_upper = e["ci_upper"], p_value = e["p_value"], windows = nrow(z), patients = length(unique(z$subject_id)))
  }
}
if (!length(results)) stop("No weight, exclusion, or subgroup fields were available; see data/README.md")
write_csv_atomic(do.call(rbind, results), file.path(out, "exclusion_subgroup_results.csv"))
write_run_metadata(cfg, c(cfg$primary_long, cfg$patient_flags), out)
