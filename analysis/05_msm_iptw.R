file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_msm) || cli$force, "MSM/IPTW")) quit(save = "no", status = 0)
out <- analysis_output_dir(cfg, "05_msm_iptw")
require_packages(c("sandwich"))
d <- prepare_primary(read_analysis_file(cfg$primary_long))
cutoff <- median(d$fluid_balance_ml, na.rm = TRUE)
d$high_balance <- as.integer(d$fluid_balance_ml >= cutoff)
vars <- c("high_balance", "age", "sex", "etiology", "iap_current", "apache_ii", "iap_next", "iap15_next", "day", "subject_id")
d <- d[complete.cases(d[vars]), , drop = FALSE]

num <- glm(high_balance ~ age + sex + etiology, data = d, family = binomial())
den <- glm(high_balance ~ age + sex + etiology + iap_current + apache_ii, data = d, family = binomial())
pn <- pmin(pmax(predict(num, type = "response"), .001), .999)
pd <- pmin(pmax(predict(den, type = "response"), .001), .999)
d$stabilized_weight <- ifelse(d$high_balance == 1, pn / pd, (1 - pn) / (1 - pd))
limits <- quantile(d$stabilized_weight, c(.01, .99))
d$truncated_weight <- pmin(pmax(d$stabilized_weight, limits[1]), limits[2])

weighted_mean <- function(x, w) sum(w * x) / sum(w)
weighted_var <- function(x, w) sum(w * (x - weighted_mean(x, w))^2) / sum(w)
smd <- function(x, a, w) {
  m1 <- weighted_mean(x[a == 1], w[a == 1]); m0 <- weighted_mean(x[a == 0], w[a == 0])
  v <- (weighted_var(x[a == 1], w[a == 1]) + weighted_var(x[a == 0], w[a == 0])) / 2
  (m1 - m0) / sqrt(v)
}
X <- model.matrix(~ age + sex + etiology + iap_current + apache_ii - 1, data = d)
balance <- do.call(rbind, lapply(seq_len(ncol(X)), function(j) data.frame(
  covariate = colnames(X)[j], smd_unweighted = smd(X[, j], d$high_balance, rep(1, nrow(d))),
  smd_weighted = smd(X[, j], d$high_balance, d$stabilized_weight),
  smd_weighted_truncated = smd(X[, j], d$high_balance, d$truncated_weight)
)))

extract_weighted <- function(outcome, weight_name, binary) {
  f <- as.formula(paste(outcome, "~ high_balance + factor(day)"))
  fit <- if (binary) glm(f, data = d, family = binomial(), weights = d[[weight_name]]) else lm(f, data = d, weights = d[[weight_name]])
  b <- coef(fit)["high_balance"]
  naive_se <- sqrt(vcov(fit)["high_balance", "high_balance"])
  robust_v <- sandwich::vcovCL(fit, cluster = d$subject_id, type = "HC0")
  robust_se <- sqrt(robust_v["high_balance", "high_balance"])
  make <- function(se, inference) {
    lo <- b - 1.96 * se; hi <- b + 1.96 * se
    data.frame(outcome = outcome, weights = weight_name, inference = inference,
               effect_measure = ifelse(binary, "OR high vs low balance", "mean difference, mmHg"),
               estimate = ifelse(binary, exp(b), b), ci_lower = ifelse(binary, exp(lo), lo),
               ci_upper = ifelse(binary, exp(hi), hi), p_value = 2 * pnorm(-abs(b / se)),
               windows = nrow(d), patients = length(unique(d$subject_id)), cutoff_ml = cutoff)
  }
  list(rows = rbind(make(naive_se, "model-based"), make(robust_se, "patient-cluster robust")), fit = fit)
}

items <- list(); fits <- list(); k <- 0L
for (w in c("stabilized_weight", "truncated_weight")) for (o in c("iap_next", "iap15_next")) {
  z <- extract_weighted(o, w, o != "iap_next"); k <- k + 1L; items[[k]] <- z$rows; fits[[paste(w, o, sep = "__")]] <- z$fit
}
results <- do.call(rbind, items)
weight_summary <- data.frame(
  weight = c("stabilized", "truncated"), mean = c(mean(d$stabilized_weight), mean(d$truncated_weight)),
  sd = c(sd(d$stabilized_weight), sd(d$truncated_weight)), min = c(min(d$stabilized_weight), min(d$truncated_weight)),
  p01 = c(quantile(d$stabilized_weight, .01), quantile(d$truncated_weight, .01)),
  median = c(median(d$stabilized_weight), median(d$truncated_weight)),
  p99 = c(quantile(d$stabilized_weight, .99), quantile(d$truncated_weight, .99)),
  max = c(max(d$stabilized_weight), max(d$truncated_weight))
)
write_csv_atomic(results, file.path(out, "msm_results.csv")); write_csv_atomic(balance, file.path(out, "covariate_balance.csv"))
write_csv_atomic(weight_summary, file.path(out, "weight_diagnostics.csv")); saveRDS(list(num = num, den = den, outcomes = fits), file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, cfg$primary_long, out, list(exposure_cutoff_ml = cutoff))
