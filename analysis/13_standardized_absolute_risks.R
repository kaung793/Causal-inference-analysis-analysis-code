file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_revision_sensitivities) || cli$force, "standardized absolute risks")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "13_standardized_absolute_risks")
d <- prepare_primary(read_analysis_file(cfg$primary_long))
cutoff <- median(d$estimated_fluid_balance_ml, na.rm = TRUE)
d$above_median <- as.integer(d$estimated_fluid_balance_ml >= cutoff)
d$day_factor <- factor(d$day)
need <- c("subject_id", "above_median", "iap15_next", "age", "sex", "etiology", "iap_current", "apache_ii", "day_factor")
z <- droplevels(d[complete.cases(d[need]), , drop = FALSE])

estimate_weighted_risks <- function(data) {
  numerator <- glm(above_median ~ age + sex + etiology, data = data, family = binomial())
  denominator <- glm(above_median ~ age + sex + etiology + iap_current + apache_ii, data = data, family = binomial())
  pn <- pmin(pmax(predict(numerator, type = "response"), 1e-6), 1 - 1e-6)
  pd <- pmin(pmax(predict(denominator, type = "response"), 1e-6), 1 - 1e-6)
  data$stabilized_weight <- ifelse(data$above_median == 1, pn / pd, (1 - pn) / (1 - pd))
  outcome <- suppressWarnings(glm(iap15_next ~ above_median + day_factor, data = data,
                                  family = binomial(), weights = stabilized_weight))
  below <- data; below$above_median <- 0L
  above <- data; above$above_median <- 1L
  risk_below <- mean(predict(outcome, newdata = below, type = "response"))
  risk_above <- mean(predict(outcome, newdata = above, type = "response"))
  values <- c(risk_below = risk_below, risk_above = risk_above,
              risk_difference = risk_above - risk_below, risk_ratio = risk_above / risk_below,
              odds_ratio = exp(unname(coef(outcome)["above_median"])))
  list(values = values, weights = data$stabilized_weight, models = list(numerator = numerator, denominator = denominator, outcome = outcome))
}

point_object <- estimate_weighted_risks(z); point <- point_object$values
set.seed(cfg$seed); ids <- unique(z$subject_id)
boot <- matrix(NA_real_, cfg$bootstrap_reps, length(point), dimnames = list(NULL, names(point)))
for (i in seq_len(cfg$bootstrap_reps)) {
  sampled <- sample(ids, length(ids), replace = TRUE)
  pieces <- lapply(seq_along(sampled), function(j) {
    x <- z[z$subject_id == sampled[j], , drop = FALSE]
    x$bootstrap_subject <- j
    x
  })
  b <- do.call(rbind, pieces)
  boot[i, ] <- tryCatch(estimate_weighted_risks(b)$values, error = function(e) rep(NA_real_, length(point)))
}
if (sum(complete.cases(boot)) < .8 * cfg$bootstrap_reps) stop("Too many failed standardized-risk bootstrap replicates")
rows <- data.frame(
  estimand = names(point), estimate = as.numeric(point),
  ci_lower = apply(boot, 2, quantile, .025, na.rm = TRUE),
  ci_upper = apply(boot, 2, quantile, .975, na.rm = TRUE),
  patients = length(ids), windows = nrow(z), cutoff_ml = cutoff,
  successful_replicates = colSums(is.finite(boot)),
  method = "stabilized-IPTW weighted logistic model; patient-cluster bootstrap"
)
weight_summary <- data.frame(mean = mean(point_object$weights), sd = sd(point_object$weights),
                             min = min(point_object$weights), median = median(point_object$weights), max = max(point_object$weights))
write_csv_atomic(rows, file.path(out, "standardized_absolute_risks.csv"))
write_csv_atomic(weight_summary, file.path(out, "weight_diagnostics.csv"))
saveRDS(point_object$models, file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out, list(bootstrap_unit = "patient", bootstrap_reps = cfg$bootstrap_reps, exposure_cutoff_ml = cutoff))
