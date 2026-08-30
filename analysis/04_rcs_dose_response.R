file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
out <- analysis_output_dir(cfg, "04_rcs_dose_response")
require_packages(c("splines", "zoo"))
input_path <- cfg$rcs_long %||% cfg$primary_long
d <- prepare_primary(read_analysis_file(input_path))
needed <- c("fluid_intake_ml", "estimated_fluid_balance_ml", "iap_next", "iap15_next", "iap20_next",
            "iap_current", "apache_ii", "age", "sex", "etiology", "creatinine", "map")
d <- d[complete.cases(d[needed]), , drop = FALSE]

mode_value <- function(x) names(which.max(table(x)))[1L]
reference_data <- function(data, exposure, grid) {
  z <- data.frame(matrix(nrow = length(grid), ncol = 0))
  z[[exposure]] <- grid
  for (v in c("iap_current", "apache_ii", "age", "creatinine", "map")) z[[v]] <- median(data[[v]])
  z$sex <- factor(mode_value(data$sex), levels = levels(data$sex))
  z$etiology <- factor(mode_value(data$etiology), levels = levels(data$etiology))
  z
}

fit_rcs <- function(data, exposure, outcome, binary) {
  knots <- unname(quantile(data[[exposure]], c(.05, .35, .65, .95)))
  spline_term <- sprintf("splines::ns(%s, knots=c(%.12g,%.12g), Boundary.knots=c(%.12g,%.12g))",
                         exposure, knots[2], knots[3], knots[1], knots[4])
  covars <- "iap_current + apache_ii + age + sex + etiology + creatinine + map"
  f_spline <- as.formula(paste(outcome, "~", spline_term, "+", covars))
  f_linear <- as.formula(paste(outcome, "~", exposure, "+", covars))
  f_null <- as.formula(paste(outcome, "~", covars))
  if (binary) {
    fit <- glm(f_spline, data = data, family = binomial())
    linear <- glm(f_linear, data = data, family = binomial())
  } else {
    fit <- lm(f_spline, data = data); linear <- lm(f_linear, data = data)
  }
  p_nonlin <- anova(linear, fit, test = "Chisq")$`Pr(>Chi)`[2]
  spline_rows <- grep("^splines::ns\\(", row.names(summary(fit)$coefficients))
  p_overall <- if (length(spline_rows)) summary(fit)$coefficients[spline_rows[1], 4] else NA_real_
  grid <- seq(knots[1], knots[4], length.out = 200L)
  nd <- reference_data(data, exposure, grid)
  ref <- reference_data(data, exposure, median(data[[exposure]]))
  tt <- delete.response(terms(fit))
  X <- model.matrix(tt, nd); X0 <- model.matrix(tt, ref)
  L <- X - matrix(X0[1, ], nrow(X), ncol(X), byrow = TRUE)
  b <- coef(fit); V <- vcov(fit)
  lp <- as.numeric(L %*% b)
  se <- sqrt(pmax(0, rowSums((L %*% V) * L)))
  effect <- if (binary) exp(lp) else lp
  lower <- if (binary) exp(lp - 1.96 * se) else lp - 1.96 * se
  upper <- if (binary) exp(lp + 1.96 * se) else lp + 1.96 * se
  h <- grid[2] - grid[1]
  first <- c(NA, diff(lp) / h)
  second <- c(NA, diff(first) / h)
  first_smooth <- zoo::rollmean(first, k = 5, fill = NA, align = "center")
  second_smooth <- zoo::rollmean(second, k = 5, fill = NA, align = "center")
  valid <- which(is.finite(second_smooth))
  inflection <- if (p_nonlin < .05 && length(valid)) grid[valid[which.max(abs(second_smooth[valid]))]] else NA_real_
  slope_threshold <- quantile(abs(first_smooth), .25, na.rm = TRUE)
  plateau_candidates <- which(grid > median(grid) & is.finite(first_smooth) & abs(first_smooth) <= slope_threshold)
  plateau <- if (p_nonlin < .05 && length(plateau_candidates)) grid[plateau_candidates[1]] else NA_real_
  curve <- data.frame(exposure = exposure, outcome = outcome, exposure_ml = grid,
                      effect = effect, ci_lower = lower, ci_upper = upper,
                      first_derivative = first, second_derivative = second,
                      first_derivative_smooth = first_smooth, second_derivative_smooth = second_smooth)
  summary <- data.frame(exposure = exposure, outcome = outcome, n = nrow(data),
                        knot_05 = knots[1], knot_35 = knots[2], knot_65 = knots[3], knot_95 = knots[4],
                        overall_p = p_overall, nonlinearity_p = p_nonlin,
                        exploratory_inflection_ml = inflection, exploratory_plateau_ml = plateau,
                        interpretation = ifelse(p_nonlin < .05, "nonlinear; curve-derived values are descriptive", "no evidence of nonlinearity"))
  list(summary = summary, curve = curve, fit = fit, linear = linear)
}

specs <- expand.grid(exposure = c("fluid_intake_ml", "estimated_fluid_balance_ml"),
                     outcome = c("iap_next", "iap15_next", "iap20_next"), stringsAsFactors = FALSE)
ans <- lapply(seq_len(nrow(specs)), function(i) fit_rcs(d, specs$exposure[i], specs$outcome[i], specs$outcome[i] != "iap_next"))
summary <- do.call(rbind, lapply(ans, `[[`, "summary")); curves <- do.call(rbind, lapply(ans, `[[`, "curve"))
write_csv_atomic(summary, file.path(out, "rcs_tests_and_exploratory_inflections.csv"))
write_csv_atomic(curves, file.path(out, "rcs_prediction_curves.csv"))
saveRDS(lapply(ans, function(x) x[c("fit", "linear")]), file.path(out, "fitted_models.rds"), compress = "xz")
write_run_metadata(cfg, input_path, out, list(complete_case_windows = nrow(d), complete_case_patients = length(unique(d$subject_id))))
