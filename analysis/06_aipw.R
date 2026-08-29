file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_aipw) || cli$force, "AIPW")) quit(save = "no", status = 0)
# This is the portable current-input implementation. The second-round
# manuscript's frozen AIPW analysis is reproduced separately by
# release/01_manuscript_aipw.R because it used a dedicated input freeze and
# patient-cluster bootstrap resampling.
out <- analysis_output_dir(cfg, "06_aipw")
d0 <- prepare_primary(read_analysis_file(cfg$primary_long))
cutoff <- median(d0$estimated_fluid_balance_ml, na.rm = TRUE)
d0$high_balance <- as.integer(d0$estimated_fluid_balance_ml >= cutoff)

estimate_aipw <- function(data, covars) {
  needed <- c("high_balance", "iap_next", "iap15_next", covars)
  d <- data[complete.cases(data[needed]), , drop = FALSE]
  rhs <- paste(covars, collapse = " + ")
  ps_fit <- glm(as.formula(paste("high_balance ~", rhs)), data = d, family = binomial())
  ps <- pmin(pmax(predict(ps_fit, type = "response"), .001), .999)
  A <- d$high_balance
  rows <- list(); models <- list(ps = ps_fit)
  for (outcome in c("iap_next", "iap15_next")) {
    binary <- outcome == "iap15_next"
    f <- as.formula(paste(outcome, "~ high_balance +", rhs))
    om <- if (binary) glm(f, data = d, family = binomial()) else lm(f, data = d)
    d1 <- d; d1$high_balance <- 1L; d0x <- d; d0x$high_balance <- 0L
    mu1 <- predict(om, newdata = d1, type = if (binary) "response" else "response")
    mu0 <- predict(om, newdata = d0x, type = if (binary) "response" else "response")
    Y <- d[[outcome]]
    phi1 <- mu1 + A * (Y - mu1) / ps
    phi0 <- mu0 + (1 - A) * (Y - mu0) / (1 - ps)
    ey1 <- mean(phi1); ey0 <- mean(phi0)
    if (binary) {
      vals <- c(risk_high = ey1, risk_low = ey0, risk_difference = ey1 - ey0,
                risk_ratio = ey1 / ey0, odds_ratio = (ey1 / (1 - ey1)) / (ey0 / (1 - ey0)))
    } else vals <- c(mean_high = ey1, mean_low = ey0, ate = ey1 - ey0)
    rows[[outcome]] <- data.frame(outcome = outcome, estimand = names(vals), estimate = as.numeric(vals),
                                  windows = nrow(d), patients = length(unique(d$subject_id)))
    models[[outcome]] <- om
  }
  list(results = do.call(rbind, rows), models = models, n = nrow(d))
}

scenarios <- list(
  primary = c("age", "sex", "etiology", "iap_current", "apache_ii"),
  extended = c("age", "sex", "etiology", "iap_current", "apache_ii", "creatinine", "map")
)
point <- lapply(scenarios, function(v) estimate_aipw(d0, v))
point_rows <- do.call(rbind, Map(function(x, nm) transform(x$results, scenario = nm), point, names(point)))

bootstrap_once <- function(data, covars, unit) {
  if (unit == "row") {
    b <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
  } else {
    ids <- unique(data$subject_id); sampled <- sample(ids, length(ids), replace = TRUE)
    b <- do.call(rbind, lapply(seq_along(sampled), function(i) {
      z <- data[data$subject_id == sampled[i], , drop = FALSE]; z$subject_id <- paste0("boot_", i); z
    }))
  }
  estimate_aipw(b, covars)$results[c("outcome", "estimand", "estimate")]
}

set.seed(cfg$seed)
boot_rows <- list(); idx <- 0L
for (scenario in names(scenarios)) for (unit in cfg$aipw_bootstrap_units) {
  reps <- vector("list", cfg$bootstrap_reps)
  for (b in seq_len(cfg$bootstrap_reps)) reps[[b]] <- tryCatch(bootstrap_once(d0, scenarios[[scenario]], unit), error = function(e) NULL)
  good <- Filter(Negate(is.null), reps)
  if (length(good) < .8 * cfg$bootstrap_reps) stop("Too many failed AIPW bootstrap replicates")
  long <- do.call(rbind, Map(function(z, i) transform(z, replicate = i), good, seq_along(good)))
  key <- interaction(long$outcome, long$estimand, drop = TRUE)
  ci <- do.call(rbind, lapply(split(long, key), function(z) data.frame(
    outcome = z$outcome[1], estimand = z$estimand[1], ci_lower = quantile(z$estimate, .025),
    ci_upper = quantile(z$estimate, .975), successful_replicates = nrow(z)
  )))
  idx <- idx + 1L; boot_rows[[idx]] <- transform(ci, scenario = scenario, bootstrap_unit = unit)
}
ci <- do.call(rbind, boot_rows)
results <- merge(point_rows, ci, by = c("scenario", "outcome", "estimand"), all.x = TRUE)
write_csv_atomic(results, file.path(out, "aipw_results.csv"))
saveRDS(lapply(point, `[[`, "models"), file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out, list(exposure_cutoff_ml = cutoff, bootstrap_reps = cfg$bootstrap_reps))
