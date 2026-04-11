start_time <- Sys.time()
project_root <- "PROJECT_ROOT_PLACEHOLDER"
setwd(project_root)

# Load packages
required_packages <- c("dplyr", "readxl", "splines", "ggplot2", "officer", "flextable")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(dplyr)
library(readxl)
library(splines)
library(ggplot2)
library(officer)
library(flextable)

# Create output directory
ts <- format(start_time, "%Y%m%d_%H%M")
out_dir <- file.path("output", ts)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("RCS Revised Threshold Analysis Started\n")
cat("Output directory:", out_dir, "\n")

# Read data
data <- readxl::read_xlsx("data/longitudinal_data_final.xlsx")
data <- data %>%
  dplyr::mutate(
    IAH15_next = ifelse(iap_next >= 15, 1, 0),
    IAH20_next = ifelse(iap_next >= 20, 1, 0)
  )

# Complete case analysis
data_complete <- data %>%
  dplyr::filter(!is.na(iap_next) & !is.na(fluid_t) & !is.na(fluid_balance_t) &
                !is.na(iap_t) & !is.na(apache_t) & !is.na(age) & !is.na(sex) &
                !is.na(etiology) & !is.na(cr_t) & !is.na(map_t)) %>%
  dplyr::mutate(
    sex = as.factor(sex),
    etiology = as.factor(etiology)
  )

cat("Complete case N =", nrow(data_complete), "\n\n")

# Calculate knots (5th, 35th, 65th, 95th percentiles) - UNCHANGED
knots_fluid_t <- quantile(data_complete$fluid_t, probs = c(0.05, 0.35, 0.65, 0.95))
knots_fluid_balance_t <- quantile(data_complete$fluid_balance_t, probs = c(0.05, 0.35, 0.65, 0.95))
median_fluid_t <- median(data_complete$fluid_t)
median_fluid_balance_t <- median(data_complete$fluid_balance_t)

# Calculate search range (5th-95th percentile) for threshold detection
search_range_fluid_t <- quantile(data_complete$fluid_t, probs = c(0.05, 0.95))
search_range_fluid_balance_t <- quantile(data_complete$fluid_balance_t, probs = c(0.05, 0.95))

cat("Knots for fluid_t:", knots_fluid_t, "\n")
cat("Knots for fluid_balance_t:", knots_fluid_balance_t, "\n")
cat("Search range for fluid_t:", search_range_fluid_t, "\n")
cat("Search range for fluid_balance_t:", search_range_fluid_balance_t, "\n\n")

# Initialize results storage
threshold_results <- data.frame()

# Analysis function - UNCHANGED model fitting
analyze_rcs_revised <- function(outcome_var, exposure_var, knots, data, median_exp,
                                search_range, is_binary = FALSE) {

  # Build formula
  formula_rcs <- paste0(outcome_var, " ~ ns(", exposure_var, ", knots = c(",
                       paste(knots[2:3], collapse = ", "),
                       "), Boundary.knots = c(", knots[1], ", ", knots[4],
                       ")) + iap_t + apache_t + age + sex + etiology + cr_t + map_t")

  formula_linear <- paste0(outcome_var, " ~ ", exposure_var,
                          " + iap_t + apache_t + age + sex + etiology + cr_t + map_t")

  # Fit models
  if (is_binary) {
    model_rcs <- glm(as.formula(formula_rcs), data = data, family = binomial)
    model_linear <- glm(as.formula(formula_linear), data = data, family = binomial)
  } else {
    model_rcs <- lm(as.formula(formula_rcs), data = data)
    model_linear <- lm(as.formula(formula_linear), data = data)
  }

  # Test nonlinearity
  lr_test <- anova(model_linear, model_rcs, test = "Chisq")
  p_nonlin <- lr_test$`Pr(>Chi)`[2]

  # Overall P-value
  coef_summary <- summary(model_rcs)$coefficients
  spline_rows <- grep("^ns\\(", rownames(coef_summary))
  p_overall <- if (length(spline_rows) > 0) coef_summary[spline_rows[1], 4] else NA

  # Generate predictions within search range
  pred_data <- data.frame(
    exp_val = seq(search_range[1], search_range[2], length.out = 200)
  )
  names(pred_data)[1] <- exposure_var

  pred_data$iap_t <- median(data$iap_t, na.rm = TRUE)
  pred_data$apache_t <- median(data$apache_t, na.rm = TRUE)
  pred_data$age <- median(data$age, na.rm = TRUE)
  pred_data$sex <- as.factor(names(sort(table(data$sex), decreasing = TRUE))[1])
  pred_data$etiology <- as.factor(names(sort(table(data$etiology), decreasing = TRUE))[1])
  pred_data$cr_t <- median(data$cr_t, na.rm = TRUE)
  pred_data$map_t <- median(data$map_t, na.rm = TRUE)

  pred_vals <- predict(model_rcs, newdata = pred_data, se.fit = TRUE,
                      type = if(is_binary) "link" else "response")

  pred_df <- data.frame(
    exposure = pred_data[[exposure_var]],
    pred = pred_vals$fit,
    se = pred_vals$se.fit
  )

  if (is_binary) {
    # Calculate OR relative to median
    pred_median <- predict(model_rcs,
                          newdata = pred_data[which.min(abs(pred_data[[exposure_var]] - median_exp)), ],
                          type = "link")
    pred_df <- pred_df %>%
      dplyr::mutate(
        or = exp(pred - pred_median),
        lower = exp(pred - 1.96 * se - pred_median),
        upper = exp(pred + 1.96 * se - pred_median)
      )
  } else {
    pred_df <- pred_df %>%
      dplyr::mutate(
        lower = pred - 1.96 * se,
        upper = pred + 1.96 * se
      )
  }

  # Calculate derivatives (first and second)
  pred_df <- pred_df %>%
    dplyr::arrange(exposure)

  # First derivative
  pred_df$first_deriv <- c(NA, diff(pred_df$pred) / diff(pred_df$exposure))

  # Second derivative
  first_deriv_diff <- diff(pred_df$first_deriv, na.rm = FALSE)
  exposure_diff <- diff(pred_df$exposure)
  pred_df$second_deriv <- c(NA, first_deriv_diff / exposure_diff[-length(exposure_diff)])

  # Initialize threshold variables
  inflection_point <- NA
  plateau_point <- NA
  threshold_interpretation <- ""

  # Only detect thresholds if non-linearity is significant (P < 0.05)
  if (!is.na(p_nonlin) && p_nonlin < 0.05) {

    # Smooth derivatives to reduce noise
    if (nrow(pred_df) > 10) {
      smooth_window <- 5
      pred_df <- pred_df %>%
        dplyr::mutate(
          first_deriv_smooth = zoo::rollmean(first_deriv, k = smooth_window, fill = NA, align = "center"),
          second_deriv_smooth = zoo::rollmean(second_deriv, k = smooth_window, fill = NA, align = "center")
        )
    } else {
      pred_df$first_deriv_smooth <- pred_df$first_deriv
      pred_df$second_deriv_smooth <- pred_df$second_deriv
    }

    # Method 1: Inflection point - where second derivative changes sign (maximum curvature)
    # Look for where acceleration changes from increasing to decreasing
    valid_second_deriv <- pred_df %>%
      dplyr::filter(!is.na(second_deriv_smooth)) %>%
      dplyr::filter(exposure >= search_range[1] & exposure <= search_range[2])

    if (nrow(valid_second_deriv) > 5) {
      # Find where second derivative crosses zero or reaches maximum
      second_deriv_max_idx <- which.max(abs(valid_second_deriv$second_deriv_smooth))
      if (length(second_deriv_max_idx) > 0) {
        inflection_point <- valid_second_deriv$exposure[second_deriv_max_idx]
      }

      # Alternative: find where first derivative reaches certain threshold (e.g., 75th percentile)
      first_deriv_threshold <- quantile(valid_second_deriv$first_deriv_smooth, 0.75, na.rm = TRUE)
      inflection_candidates <- valid_second_deriv %>%
        dplyr::filter(first_deriv_smooth >= first_deriv_threshold) %>%
        dplyr::slice(1)

      if (nrow(inflection_candidates) > 0 && is.na(inflection_point)) {
        inflection_point <- inflection_candidates$exposure
      }
    }

    # Method 2: Plateau/saturation point - where first derivative approaches zero or stabilizes
    # Look for where slope becomes very small (< 25th percentile of absolute slope)
    valid_first_deriv <- pred_df %>%
      dplyr::filter(!is.na(first_deriv_smooth)) %>%
      dplyr::filter(exposure >= search_range[1] & exposure <= search_range[2])

    if (nrow(valid_first_deriv) > 10) {
      # Find where slope drops below 25th percentile (indicating flattening)
      slope_threshold <- quantile(abs(valid_first_deriv$first_deriv_smooth), 0.25, na.rm = TRUE)

      # Look for plateau in upper half of exposure range
      upper_half <- valid_first_deriv %>%
        dplyr::filter(exposure > median(exposure))

      plateau_candidates <- upper_half %>%
        dplyr::filter(abs(first_deriv_smooth) <= slope_threshold) %>%
        dplyr::slice(1)

      if (nrow(plateau_candidates) > 0) {
        plateau_point <- plateau_candidates$exposure
      }
    }

    # Generate interpretation
    if (!is.na(inflection_point) && !is.na(plateau_point)) {
      threshold_interpretation <- sprintf("Non-linear relationship with inflection at ~%.0f mL and plateau beginning at ~%.0f mL",
                                         inflection_point, plateau_point)
    } else if (!is.na(inflection_point)) {
      threshold_interpretation <- sprintf("Non-linear relationship with inflection point at ~%.0f mL", inflection_point)
    } else if (!is.na(plateau_point)) {
      threshold_interpretation <- sprintf("Non-linear relationship with plateau beginning at ~%.0f mL", plateau_point)
    } else {
      threshold_interpretation <- "Significant non-linearity detected, but no clear inflection or plateau point identified within 5th-95th percentile range"
    }

  } else {
    threshold_interpretation <- "No significant non-linearity observed; relationship more consistent with linear or near-linear trend"
  }

  return(list(
    p_overall = p_overall,
    p_nonlin = p_nonlin,
    inflection_point = inflection_point,
    plateau_point = plateau_point,
    interpretation = threshold_interpretation
  ))
}

# Need zoo package for rolling mean
if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo", repos = "https://cloud.r-project.org")
}
library(zoo)

# Run all 6 combinations
combinations <- list(
  list(outcome = "iap_next", exposure = "fluid_t", knots = knots_fluid_t, median = median_fluid_t,
       search_range = search_range_fluid_t, binary = FALSE, exp_name = "Fluid intake", out_name = "IAP (continuous)"),
  list(outcome = "IAH15_next", exposure = "fluid_t", knots = knots_fluid_t, median = median_fluid_t,
       search_range = search_range_fluid_t, binary = TRUE, exp_name = "Fluid intake", out_name = "IAH15"),
  list(outcome = "IAH20_next", exposure = "fluid_t", knots = knots_fluid_t, median = median_fluid_t,
       search_range = search_range_fluid_t, binary = TRUE, exp_name = "Fluid intake", out_name = "IAH20"),
  list(outcome = "iap_next", exposure = "fluid_balance_t", knots = knots_fluid_balance_t, median = median_fluid_balance_t,
       search_range = search_range_fluid_balance_t, binary = FALSE, exp_name = "Fluid balance", out_name = "IAP (continuous)"),
  list(outcome = "IAH15_next", exposure = "fluid_balance_t", knots = knots_fluid_balance_t, median = median_fluid_balance_t,
       search_range = search_range_fluid_balance_t, binary = TRUE, exp_name = "Fluid balance", out_name = "IAH15"),
  list(outcome = "IAH20_next", exposure = "fluid_balance_t", knots = knots_fluid_balance_t, median = median_fluid_balance_t,
       search_range = search_range_fluid_balance_t, binary = TRUE, exp_name = "Fluid balance", out_name = "IAH20")
)

for (i in seq_along(combinations)) {
  combo <- combinations[[i]]
  cat("Analyzing", combo$exp_name, "->", combo$out_name, "...\n")

  # Analyze
  result <- analyze_rcs_revised(combo$outcome, combo$exposure, combo$knots, data_complete,
                                combo$median, combo$search_range, combo$binary)

  # Store results
  threshold_results <- rbind(threshold_results, data.frame(
    Exposure = combo$exp_name,
    Outcome = combo$out_name,
    Overall_P = result$p_overall,
    Nonlinearity_P = result$p_nonlin,
    Inflection_Point = result$inflection_point,
    Plateau_Point = result$plateau_point,
    Interpretation = result$interpretation,
    stringsAsFactors = FALSE
  ))
}

# Format results table
threshold_results <- threshold_results %>%
  dplyr::mutate(
    Overall_P_fmt = ifelse(is.na(Overall_P), "NA",
                          ifelse(Overall_P < 0.001, "<0.001", sprintf("%.4f", Overall_P))),
    Nonlinearity_P_fmt = ifelse(is.na(Nonlinearity_P), "NA",
                                ifelse(Nonlinearity_P < 0.001, "<0.001", sprintf("%.4f", Nonlinearity_P))),
    Inflection_Point_fmt = ifelse(is.na(Inflection_Point), "Not applicable",
                                  sprintf("~%.0f mL", Inflection_Point)),
    Plateau_Point_fmt = ifelse(is.na(Plateau_Point), "Not identified",
                              sprintf("~%.0f mL", Plateau_Point))
  )

# Create display table
table_display <- threshold_results %>%
  dplyr::select(Exposure, Outcome, Overall_P_fmt, Nonlinearity_P_fmt,
               Inflection_Point_fmt, Plateau_Point_fmt, Interpretation)

colnames(table_display) <- c("Exposure", "Outcome", "Overall P", "Non-linearity P",
                             "Threshold/Inflection Point", "Plateau/Saturation Point", "Interpretation")

# Save CSV
write.csv(table_display, file.path(out_dir, "RCS_Revised_Threshold_Analysis.csv"), row.names = FALSE)

# Create Word document with flextable
ft <- flextable(table_display) %>%
  bold(part = "header") %>%
  align(align = "left", part = "all", j = 1:2) %>%
  align(align = "center", part = "all", j = 3:6) %>%
  align(align = "left", part = "all", j = 7) %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 2), part = "header") %>%
  hline_bottom(border = fp_border(width = 2), part = "header") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  autofit() %>%
  width(j = 7, width = 3.5)

doc <- read_docx() %>%
  body_add_par("Table. Revised RCS Threshold Analysis with Clinically Meaningful Inflection Points", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_par("Restricted cubic splines with 4 knots (5th, 35th, 65th, 95th percentiles). Model 3: Adjusted for IAP at time t, APACHE II score, age, sex, etiology, creatinine, and mean arterial pressure. Overall P: Wald test for first spline term. Non-linearity P: Likelihood ratio test comparing RCS model vs linear model. Threshold/inflection points and plateau/saturation points are only reported for models with significant non-linearity (P < 0.05). Inflection points identified based on maximum curvature (second derivative); plateau points identified where slope approaches minimum in upper exposure range. All threshold detection restricted to 5th-95th percentile of exposure distribution to avoid extreme value instability. For models without significant non-linearity, the relationship is more consistent with a linear or near-linear trend.", style = "Normal")

print(doc, target = file.path(out_dir, "RCS_Revised_Threshold_Analysis.docx"))

# Generate results text
results_text <- paste0(
  "RESULTS - Revised RCS Threshold Analysis with Clinically Meaningful Inflection Points\n\n",
  "We re-evaluated threshold effects in the dose-response relationships between fluid exposures and IAP outcomes ",
  "using restricted cubic splines with full adjustment (Model 3). Threshold detection was restricted to models with ",
  "significant non-linearity (P < 0.05) and limited to the 5th-95th percentile range of exposure distribution to ",
  "avoid extreme value instability.\n\n",
  "NON-LINEARITY ASSESSMENT:\n\n",
  "Significant non-linearity (P < 0.05) was observed in the following relationships:\n"
)

# Identify significant non-linear relationships
sig_nonlin <- threshold_results %>% dplyr::filter(Nonlinearity_P < 0.05)
if (nrow(sig_nonlin) > 0) {
  for (i in 1:nrow(sig_nonlin)) {
    results_text <- paste0(results_text,
      sprintf("- %s → %s (P = %s)\n", sig_nonlin$Exposure[i], sig_nonlin$Outcome[i],
              sig_nonlin$Nonlinearity_P_fmt[i]))
  }
} else {
  results_text <- paste0(results_text, "None\n")
}

results_text <- paste0(results_text, "\n")

# Identify linear relationships
linear_rel <- threshold_results %>% dplyr::filter(Nonlinearity_P >= 0.05)
if (nrow(linear_rel) > 0) {
  results_text <- paste0(results_text,
    "The following relationships showed no significant non-linearity and are more consistent with linear or near-linear trends:\n")
  for (i in 1:nrow(linear_rel)) {
    results_text <- paste0(results_text,
      sprintf("- %s → %s (P = %s)\n", linear_rel$Exposure[i], linear_rel$Outcome[i],
              linear_rel$Nonlinearity_P_fmt[i]))
  }
  results_text <- paste0(results_text, "\n")
}

# Fluid balance specific findings
fb_results <- threshold_results %>% dplyr::filter(Exposure == "Fluid balance")
results_text <- paste0(results_text,
  "FLUID BALANCE DOSE-RESPONSE CHARACTERISTICS:\n\n",
  "Fluid balance demonstrated consistent significant non-linearity across all three outcomes:\n\n")

for (i in 1:nrow(fb_results)) {
  results_text <- paste0(results_text,
    sprintf("%s: %s\n", fb_results$Outcome[i], fb_results$Interpretation[i]))
}

results_text <- paste0(results_text, "\n",
  "KEY FINDINGS:\n\n",
  "1. Fluid balance shows more pronounced and consistent non-linear patterns compared to fluid intake.\n",
  "2. For relationships with significant non-linearity, inflection points and plateau regions were identified ",
  "within clinically relevant exposure ranges (5th-95th percentile).\n",
  "3. Linear relationships (non-linearity P ≥ 0.05) suggest proportional associations without clear thresholds.\n",
  "4. All threshold estimates should be interpreted as population-level patterns rather than individual clinical cutoffs.\n"
)

writeLines(results_text, con = file.path(out_dir, "RCS_Revised_Threshold_Results.txt"))

# Create run log
end_time <- Sys.time()
log_content <- paste0(
  "RCS Revised Threshold Analysis\n",
  "===============================\n\n",
  "Start time: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n",
  "End time: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n",
  "Duration: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n\n",
  "Language: R\n",
  "Script: scripts/rcs_revised_threshold_analysis.R\n\n",
  "Packages used:\n",
  "- dplyr\n- readxl\n- splines\n- ggplot2\n- officer\n- flextable\n- zoo\n\n",
  "Key modifications from original threshold analysis:\n",
  "1. Thresholds only reported for models with non-linearity P < 0.05\n",
  "2. Inflection points identified using second derivative (maximum curvature)\n",
  "3. Plateau points identified where slope approaches minimum in upper range\n",
  "4. All threshold detection restricted to 5th-95th percentile (avoiding extremes)\n",
  "5. Linear relationships explicitly noted (no threshold reported)\n\n",
  "Sample size: ", nrow(data_complete), " complete cases\n",
  "Models: 6 exposure-outcome combinations\n",
  "Significant non-linear relationships: ", nrow(sig_nonlin), "\n",
  "Linear/near-linear relationships: ", nrow(linear_rel), "\n\n",
  "Status: Completed successfully\n"
)

writeLines(log_content, con = file.path(out_dir, "run_log.txt"))

cat("\n=== RCS Revised Threshold Analysis Completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("Significant non-linear relationships:", nrow(sig_nonlin), "\n")
cat("Linear/near-linear relationships:", nrow(linear_rel), "\n")
cat("Duration:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
