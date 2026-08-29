file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_manuscript_release) || cli$force,
                        "manuscript-release AIPW")) quit(save = "no", status = 0)

weights_path <- cfg$aipw_manuscript_weights %||% ""
long_path <- cfg$aipw_manuscript_long %||% ""
if (!nzchar(weights_path) || !nzchar(long_path)) {
  stop("Set aipw_manuscript_weights and aipw_manuscript_long in the local configuration")
}
out <- analysis_output_dir(cfg, "manuscript_release_aipw")

weights <- read_analysis_file(weights_path)
long <- read_analysis_file(long_path)
assert_columns(weights, c("subject_id", "day", "high_fluid_balance", "fluid_balance_t",
                          "age", "sex", "etiology", "iap_t", "apache_t", "IAH15_next"),
               "manuscript AIPW weight freeze")
assert_columns(long, c("subject_id", "day", "fluid_balance_t", "age", "sex", "etiology",
                       "iap_t", "apache_t", "cr_t", "map_t", "iap_next"),
               "manuscript AIPW longitudinal freeze")

weights$subject_id <- as.character(weights$subject_id)
long$subject_id <- as.character(long$subject_id)
key_weights <- paste(weights$subject_id, weights$day, sep = "::")
key_long <- paste(long$subject_id, long$day, sep = "::")
if (anyDuplicated(key_long)) stop("AIPW longitudinal freeze has duplicate subject-day rows")
match_long <- match(key_weights, key_long)
if (anyNA(match_long)) stop("AIPW weight rows could not all be matched to the longitudinal freeze")

primary <- weights
primary$iap_next <- long$iap_next[match_long]
primary <- primary[complete.cases(primary[c("high_fluid_balance", "IAH15_next", "iap_next",
                                             "age", "sex", "etiology", "iap_t", "apache_t")]), , drop = FALSE]

extended <- long[long$day %in% c(1, 2, 3, 5), , drop = FALSE]
extended$high_fluid_balance <- as.integer(extended$fluid_balance_t >= 1010)
extended$IAH15_next <- as.integer(extended$iap_next >= 15)
extended <- extended[complete.cases(extended[c("high_fluid_balance", "IAH15_next", "iap_next",
                                               "age", "sex", "etiology", "iap_t", "apache_t",
                                               "cr_t", "map_t")]), , drop = FALSE]

if (isTRUE(cfg$strict_cohort_checks)) {
  dims <- c(primary_patients = length(unique(primary$subject_id)), primary_windows = nrow(primary),
            extended_patients = length(unique(extended$subject_id)), extended_windows = nrow(extended))
  expected <- c(primary_patients = 773L, primary_windows = 2863L,
                extended_patients = 771L, extended_windows = 2188L)
  if (!identical(as.integer(dims), as.integer(expected))) {
    stop("Manuscript AIPW dimensions differ from the frozen contract: ",
         paste(names(dims), dims, sep = "=", collapse = "; "))
  }
}

estimate_once <- function(data, covars, frequency_weights = rep(1, nrow(data))) {
  rhs <- paste(covars, collapse = " + ")
  ps_model <- suppressWarnings(glm(
    as.formula(paste("high_fluid_balance ~", rhs)), data = data,
    weights = frequency_weights, family = binomial(link = "logit")
  ))
  ps <- fitted(ps_model)
  if (any(!is.finite(ps)) || any(ps <= 0) || any(ps >= 1)) stop("Invalid propensity scores")

  binary_model <- suppressWarnings(glm(
    as.formula(paste("IAH15_next ~ high_fluid_balance +", rhs)), data = data,
    weights = frequency_weights, family = binomial(link = "logit")
  ))
  continuous_model <- lm(
    as.formula(paste("iap_next ~ high_fluid_balance +", rhs)), data = data,
    weights = frequency_weights
  )

  treated <- control <- data
  treated$high_fluid_balance <- 1L; control$high_fluid_balance <- 0L
  A <- data$high_fluid_balance

  mu1_binary <- predict(binary_model, newdata = treated, type = "response")
  mu0_binary <- predict(binary_model, newdata = control, type = "response")
  phi1_binary <- (A / ps) * data$IAH15_next + (1 - A / ps) * mu1_binary
  phi0_binary <- ((1 - A) / (1 - ps)) * data$IAH15_next +
    (1 - (1 - A) / (1 - ps)) * mu0_binary
  risk_high <- weighted.mean(phi1_binary, frequency_weights)
  risk_low <- weighted.mean(phi0_binary, frequency_weights)

  mu1_continuous <- predict(continuous_model, newdata = treated)
  mu0_continuous <- predict(continuous_model, newdata = control)
  phi1_continuous <- (A / ps) * data$iap_next + (1 - A / ps) * mu1_continuous
  phi0_continuous <- ((1 - A) / (1 - ps)) * data$iap_next +
    (1 - (1 - A) / (1 - ps)) * mu0_continuous

  c(
    risk_difference = risk_high - risk_low,
    odds_ratio = (risk_high / (1 - risk_high)) / (risk_low / (1 - risk_low)),
    ate = weighted.mean(phi1_continuous, frequency_weights) -
      weighted.mean(phi0_continuous, frequency_weights)
  )
}

cluster_frequency_weights <- function(data) {
  ids <- unique(data$subject_id)
  sampled <- sample(ids, length(ids), replace = TRUE)
  counts <- table(factor(sampled, levels = ids))
  as.numeric(counts[match(data$subject_id, ids)])
}

run_scenario <- function(name, data, covars, n_boot) {
  point <- estimate_once(data, covars)
  bootstrap <- matrix(NA_real_, nrow = n_boot, ncol = length(point),
                      dimnames = list(NULL, names(point)))
  for (i in seq_len(n_boot)) {
    bootstrap[i, ] <- tryCatch(
      estimate_once(data, covars, cluster_frequency_weights(data)),
      error = function(e) rep(NA_real_, length(point))
    )
  }
  bootstrap <- bootstrap[complete.cases(bootstrap), , drop = FALSE]
  if (nrow(bootstrap) < .95 * n_boot) stop(name, " had too many failed bootstrap replicates")
  p_from_zero <- function(x) 2 * min(mean(x >= 0), mean(x <= 0))
  data.frame(
    scenario = name,
    estimand = names(point),
    estimate = unname(point),
    ci_lower = apply(bootstrap, 2, quantile, probs = .025, na.rm = TRUE),
    ci_upper = apply(bootstrap, 2, quantile, probs = .975, na.rm = TRUE),
    p_value = c(p_from_zero(bootstrap[, "risk_difference"]),
                p_from_zero(bootstrap[, "risk_difference"]),
                p_from_zero(bootstrap[, "ate"])),
    patients = length(unique(data$subject_id)), windows = nrow(data),
    successful_replicates = nrow(bootstrap), bootstrap_unit = "patient",
    stringsAsFactors = FALSE
  )
}

n_boot <- cfg$manuscript_aipw_bootstrap_reps %||% 1000L
set.seed(cfg$manuscript_aipw_seed %||% 123L)
results <- rbind(
  run_scenario("primary", primary, c("age", "sex", "etiology", "iap_t", "apache_t"), n_boot),
  run_scenario("extended", extended,
               c("age", "sex", "etiology", "iap_t", "apache_t", "cr_t", "map_t"), n_boot)
)
write_csv_atomic(results, file.path(out, "manuscript_aipw_patient_bootstrap.csv"))
write_run_metadata(
  cfg, c(weights_path, long_path), out,
  list(exposure_cutoff_ml = 1010, bootstrap_seed = cfg$manuscript_aipw_seed %||% 123L,
       bootstrap_reps = n_boot, categorical_code_handling = "stored numeric codes, matching frozen workflow")
)
