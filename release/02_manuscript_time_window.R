file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R")); source(file.path(root, "R", "model_helpers.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_manuscript_release) || cli$force,
                        "manuscript-release time-window analysis")) quit(save = "no", status = 0)

require_packages(c("lme4", "lmerTest"))
input_path <- cfg$time_window_manuscript_long %||% cfg$primary_long
out <- analysis_output_dir(cfg, "manuscript_release_time_window")
d <- prepare_primary(read_analysis_file(input_path))
d$day_factor <- factor(d$day, levels = c(1, 2, 3, 5))
needed <- c("estimated_fluid_balance_l", "iap_next", "iap15_next", "iap_current", "apache_ii",
            "age", "sex", "etiology", "creatinine", "map", "subject_id", "day_factor")
d <- d[complete.cases(d[needed]), , drop = FALSE]
if (isTRUE(cfg$strict_cohort_checks) &&
    (nrow(d) != 2188L || length(unique(d$subject_id)) != 771L)) {
  stop("Manuscript time-window cohort must contain 2,188 windows from 771 patients")
}

rhs <- paste(
  "estimated_fluid_balance_l * day_factor + iap_current + apache_ii + age +",
  "sex + etiology + creatinine + map + (1 | subject_id)"
)
continuous_fit <- lmerTest::lmer(
  as.formula(paste("iap_next ~", rhs)), data = d, REML = FALSE,
  control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)
binary_fit <- lme4::glmer(
  as.formula(paste("iap15_next ~", rhs)), data = d, family = binomial(), nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

window_labels <- c(`1` = "D1_D2", `2` = "D2_D3", `3` = "D3_D5", `5` = "D5_D7")
extract_display <- function(fit, outcome, binary) {
  b <- fixed_coef(fit); V <- as.matrix(vcov(fit)); main <- "estimated_fluid_balance_l"
  rows <- vector("list", length(window_labels)); full <- vector("list", length(window_labels))
  for (i in seq_along(window_labels)) {
    day <- names(window_labels)[i]; term <- paste0(main, ":day_factor", day)
    linear <- unname(b[main])
    legacy_var <- V[main, main]
    terms <- setNames(1, main)
    if (i > 1L) {
      linear <- linear + unname(b[term])
      legacy_var <- legacy_var + V[term, term]
      terms <- c(terms, setNames(1, term))
    }
    legacy_se <- sqrt(legacy_var)
    full_effect <- linear_contrast(fit, terms, exponentiate = binary)
    int <- if (i == 1L) c(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_) else
      extract_effect(fit, term, exponentiate = binary)
    convert <- function(x) if (binary) exp(x) else x
    rows[[i]] <- data.frame(
      outcome = outcome, window = unname(window_labels[i]),
      effect_measure = if (binary) "OR per 1,000 mL" else "beta mmHg per 1,000 mL",
      estimate = convert(linear),
      ci_lower = convert(linear - 1.96 * legacy_se),
      ci_upper = convert(linear + 1.96 * legacy_se),
      contrast_vs_day1 = int["estimate"], contrast_ci_lower = int["ci_lower"],
      contrast_ci_upper = int["ci_upper"], p_interaction_vs_day1 = int["p_value"],
      display_interval_method = if (i == 1L) "single coefficient" else
        "historical variance-sum display interval; covariance omitted",
      stringsAsFactors = FALSE
    )
    full[[i]] <- data.frame(
      outcome = outcome, window = unname(window_labels[i]),
      estimate = full_effect["estimate"], ci_lower = full_effect["ci_lower"],
      ci_upper = full_effect["ci_upper"], p_value = full_effect["p_value"],
      interval_method = "full coefficient covariance matrix", stringsAsFactors = FALSE
    )
  }
  list(display = do.call(rbind, rows), full = do.call(rbind, full))
}

cont <- extract_display(continuous_fit, "continuous_iap", FALSE)
binary <- extract_display(binary_fit, "iap15", TRUE)
write_csv_atomic(rbind(cont$display, binary$display),
                 file.path(out, "document_display_formula_current_runtime.csv"))
write_csv_atomic(rbind(cont$full, binary$full),
                 file.path(out, "full_covariance_current_runtime.csv"))
diagnostics <- rbind(model_diagnostics(continuous_fit, "time_window_continuous_unscaled"),
                     model_diagnostics(binary_fit, "time_window_iap15_unscaled"))
write_csv_atomic(diagnostics, file.path(out, "model_diagnostics.csv"))
saveRDS(list(continuous = continuous_fit, iap15 = binary_fit),
        file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(
  cfg, input_path, out,
  list(analysis_sample = "Model 3 complete case", windows = nrow(d),
       patients = length(unique(d$subject_id)),
       document_interval_method = "historical variance-sum display; full-covariance companion also emitted")
)
