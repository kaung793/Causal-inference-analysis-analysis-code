file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."),
                      winslash = "/")
source(file.path(root, "R", "utils.R"))
source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args()
cfg <- load_config(root, cli$config)

require_packages(c("lme4", "lmerTest"))
input_path <- cfg$time_window_long %||% cfg$primary_long
out <- analysis_output_dir(cfg, "07_time_window_heterogeneity")
d <- prepare_primary(read_analysis_file(input_path))
d <- scale_covariates(d, c("iap_current", "apache_ii", "age", "creatinine", "map"))
d$day_factor <- factor(d$day, levels = c(1, 2, 3, 5))

needed <- c(
  "estimated_fluid_balance_l", "iap_next", "iap15_next", "iap_current_z",
  "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z",
  "subject_id", "day_factor"
)
d <- d[complete.cases(d[needed]), , drop = FALSE]
if (isTRUE(cfg$strict_cohort_checks) &&
    (nrow(d) != 2188L || length(unique(d$subject_id)) != 771L)) {
  stop("Time-window analysis requires 2,188 windows from 771 patients")
}

rhs <- paste(
  "estimated_fluid_balance_l * day_factor + iap_current_z + apache_ii_z + age_z +",
  "sex + etiology + creatinine_z + map_z + (1 | subject_id)"
)
fits <- list(
  continuous = lmerTest::lmer(
    as.formula(paste("iap_next ~", rhs)), data = d, REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
    )
  ),
  iap15 = lme4::glmer(
    as.formula(paste("iap15_next ~", rhs)), data = d,
    family = binomial(), nAGQ = 1,
    control = lme4::glmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
    )
  )
)

window_labels <- c(
  `1` = "Day 1 to 2", `2` = "Day 2 to 3",
  `3` = "Day 3 to 5", `5` = "Day 5 to 7"
)

extract_windows <- function(fit, outcome, binary) {
  coefficients <- fixed_coef(fit)
  covariance <- as.matrix(vcov(fit))
  main <- "estimated_fluid_balance_l"
  rows <- vector("list", length(window_labels))

  for (i in seq_along(window_labels)) {
    day <- names(window_labels)[i]
    interaction <- paste0(main, ":day_factor", day)
    linear <- unname(coefficients[main])
    interval_variance <- covariance[main, main]
    contrast_terms <- setNames(1, main)

    if (i > 1L) {
      linear <- linear + unname(coefficients[interaction])
      interval_variance <- interval_variance + covariance[interaction, interaction]
      contrast_terms <- c(contrast_terms, setNames(1, interaction))
    }

    interval_se <- sqrt(interval_variance)
    full <- linear_contrast(fit, contrast_terms, exponentiate = binary)
    interaction_result <- if (i == 1L) {
      c(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        p_value = NA_real_)
    } else {
      extract_effect(fit, interaction, exponentiate = binary)
    }
    transform_effect <- function(x) if (binary) exp(x) else x
    subset <- d[d$day_factor == day, , drop = FALSE]

    # The main CI columns follow the calculation reported in Table 4. CIs based
    # on the full coefficient covariance matrix are included alongside them.
    rows[[i]] <- data.frame(
      outcome = outcome,
      window = paste0("D", day, "_D", c(2, 3, 5, 7)[i]),
      window_label = unname(window_labels[i]),
      effect_measure = if (binary) "OR per 1,000 mL" else
        "beta mmHg per 1,000 mL",
      estimate = transform_effect(linear),
      ci_lower = transform_effect(linear - 1.96 * interval_se),
      ci_upper = transform_effect(linear + 1.96 * interval_se),
      full_covariance_ci_lower = full["ci_lower"],
      full_covariance_ci_upper = full["ci_upper"],
      p_value = full["p_value"],
      contrast_vs_day1 = interaction_result["estimate"],
      contrast_ci_lower = interaction_result["ci_lower"],
      contrast_ci_upper = interaction_result["ci_upper"],
      p_interaction_vs_day1 = interaction_result["p_value"],
      window_n = nrow(subset),
      events = if (binary) sum(subset$iap15_next) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  interaction_terms <- grep(
    "estimated_fluid_balance_l:day_factor",
    names(coefficients), value = TRUE, fixed = TRUE
  )
  list(
    results = do.call(rbind, rows),
    global = cbind(outcome = outcome, joint_wald_test(fit, interaction_terms))
  )
}

continuous <- extract_windows(fits$continuous, "continuous_iap", FALSE)
binary <- extract_windows(fits$iap15, "iap15", TRUE)
results <- rbind(continuous$results, binary$results)
global <- rbind(continuous$global, binary$global)
diagnostics <- rbind(
  model_diagnostics(fits$continuous, "time_window_continuous"),
  model_diagnostics(fits$iap15, "time_window_iap15")
)

write_csv_atomic(results, file.path(out, "time_window_results.csv"))
write_csv_atomic(global, file.path(out, "global_interaction_tests.csv"))
write_csv_atomic(diagnostics, file.path(out, "model_diagnostics.csv"))
saveRDS(fits, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(
  cfg, input_path, out,
  list(windows = nrow(d), patients = length(unique(d$subject_id)))
)
