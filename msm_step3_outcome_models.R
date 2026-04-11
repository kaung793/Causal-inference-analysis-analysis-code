# MSM Step 3: Fit Weighted Outcome Models (Primary Analysis Only)
# Date: 2026-03-30
# Task: Estimate causal effect of high fluid balance on IAH15 and IAP (continuous)

library(dplyr)

# Set working directory
setwd(project_root)

# Read baseline IPTW weights
weights_data <- read.csv("data/IPTW_baseline_model_weights.csv")

# Read longitudinal data to get iap_next
library(readxl)
longitudinal_data <- read_excel("data/longitudinal_data_final.xlsx")

# Merge iap_next into weights_data
weights_data <- weights_data %>%
  left_join(
    longitudinal_data %>% select(subject_id, day, iap_next),
    by = c("subject_id", "day")
  )

cat("\n=== MSM STEP 3: WEIGHTED OUTCOME MODELS ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Sample size:", nrow(weights_data), "person-time observations\n")
cat("Exposure: high_fluid_balance (1 = ≥1010 mL, 0 = <1010 mL)\n")
cat("Outcomes: IAH15_next (binary), iap_next (continuous)\n\n")

# Check data structure
cat("Exposure distribution:\n")
print(table(weights_data$high_fluid_balance))
cat("\n")

cat("Outcome distribution (IAH15_next):\n")
print(table(weights_data$IAH15_next, useNA = "ifany"))
cat("\n")

# ============================================================================
# PART 1: BINARY OUTCOME (IAH15_next)
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep="")
cat("PART 1: BINARY OUTCOME (IAH15_next)\n")
cat(rep("=", 70), "\n\n", sep="")

# --- Primary Analysis: Original Weights ---
cat("PRIMARY ANALYSIS: Baseline Model with Original Weights\n")
cat(rep("-", 70), "\n", sep="")

model_primary_binary <- glm(
  IAH15_next ~ high_fluid_balance + factor(day),
  data = weights_data,
  family = binomial(link = "logit"),
  weights = sw
)

summary_primary_binary <- summary(model_primary_binary)

# Extract coefficients for high_fluid_balance
coef_primary <- summary_primary_binary$coefficients["high_fluid_balance", ]
beta_primary <- coef_primary["Estimate"]
se_primary <- coef_primary["Std. Error"]
z_primary <- coef_primary["z value"]
p_primary <- coef_primary["Pr(>|z|)"]

# Calculate OR and 95% CI
or_primary <- exp(beta_primary)
ci_lower_primary <- exp(beta_primary - 1.96 * se_primary)
ci_upper_primary <- exp(beta_primary + 1.96 * se_primary)

cat("\nModel: IAH15_next ~ high_fluid_balance + factor(day)\n")
cat("Weights: Stabilized IPTW (original)\n")
cat("Family: Binomial (logit link)\n\n")

cat("Results for high_fluid_balance:\n")
cat(sprintf("  Beta coefficient: %.4f\n", beta_primary))
cat(sprintf("  Standard error: %.4f\n", se_primary))
cat(sprintf("  Z value: %.4f\n", z_primary))
cat(sprintf("  P value: %.4f\n", p_primary))
cat(sprintf("  Odds Ratio: %.3f\n", or_primary))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", ci_lower_primary, ci_upper_primary))
cat("\n")

# --- Sensitivity Analysis: Truncated Weights ---
cat("\nSENSITIVITY ANALYSIS: Baseline Model with Truncated Weights\n")
cat(rep("-", 70), "\n", sep="")

model_sens_binary <- glm(
  IAH15_next ~ high_fluid_balance + factor(day),
  data = weights_data,
  family = binomial(link = "logit"),
  weights = sw_truncated
)

summary_sens_binary <- summary(model_sens_binary)

# Extract coefficients for high_fluid_balance
coef_sens <- summary_sens_binary$coefficients["high_fluid_balance", ]
beta_sens <- coef_sens["Estimate"]
se_sens <- coef_sens["Std. Error"]
z_sens <- coef_sens["z value"]
p_sens <- coef_sens["Pr(>|z|)"]

# Calculate OR and 95% CI
or_sens <- exp(beta_sens)
ci_lower_sens <- exp(beta_sens - 1.96 * se_sens)
ci_upper_sens <- exp(beta_sens + 1.96 * se_sens)

cat("\nModel: IAH15_next ~ high_fluid_balance + factor(day)\n")
cat("Weights: Stabilized IPTW (truncated at 1st-99th percentile)\n")
cat("Family: Binomial (logit link)\n\n")

cat("Results for high_fluid_balance:\n")
cat(sprintf("  Beta coefficient: %.4f\n", beta_sens))
cat(sprintf("  Standard error: %.4f\n", se_sens))
cat(sprintf("  Z value: %.4f\n", z_sens))
cat(sprintf("  P value: %.4f\n", p_sens))
cat(sprintf("  Odds Ratio: %.3f\n", or_sens))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", ci_lower_sens, ci_upper_sens))
cat("\n")

# ============================================================================
# PART 2: CONTINUOUS OUTCOME (iap_next)
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep="")
cat("PART 2: CONTINUOUS OUTCOME (iap_next)\n")
cat(rep("=", 70), "\n\n", sep="")

# --- Primary Analysis: Original Weights ---
cat("PRIMARY ANALYSIS: Baseline Model with Original Weights\n")
cat(rep("-", 70), "\n", sep="")

model_primary_cont <- lm(
  iap_next ~ high_fluid_balance + factor(day),
  data = weights_data,
  weights = sw
)

summary_primary_cont <- summary(model_primary_cont)

# Extract coefficients for high_fluid_balance
coef_primary_cont <- summary_primary_cont$coefficients["high_fluid_balance", ]
beta_primary_cont <- coef_primary_cont["Estimate"]
se_primary_cont <- coef_primary_cont["Std. Error"]
t_primary_cont <- coef_primary_cont["t value"]
p_primary_cont <- coef_primary_cont["Pr(>|t|)"]

# Calculate 95% CI
ci_lower_primary_cont <- beta_primary_cont - 1.96 * se_primary_cont
ci_upper_primary_cont <- beta_primary_cont + 1.96 * se_primary_cont

cat("\nModel: iap_next ~ high_fluid_balance + factor(day)\n")
cat("Weights: Stabilized IPTW (original)\n")
cat("Family: Gaussian (identity link)\n\n")

cat("Results for high_fluid_balance:\n")
cat(sprintf("  Beta coefficient: %.4f mmHg\n", beta_primary_cont))
cat(sprintf("  Standard error: %.4f\n", se_primary_cont))
cat(sprintf("  t value: %.4f\n", t_primary_cont))
cat(sprintf("  P value: %.4f\n", p_primary_cont))
cat(sprintf("  Mean difference: %.3f mmHg\n", beta_primary_cont))
cat(sprintf("  95%% CI: [%.3f, %.3f] mmHg\n", ci_lower_primary_cont, ci_upper_primary_cont))
cat("\n")

# --- Sensitivity Analysis: Truncated Weights ---
cat("\nSENSITIVITY ANALYSIS: Baseline Model with Truncated Weights\n")
cat(rep("-", 70), "\n", sep="")

model_sens_cont <- lm(
  iap_next ~ high_fluid_balance + factor(day),
  data = weights_data,
  weights = sw_truncated
)

summary_sens_cont <- summary(model_sens_cont)

# Extract coefficients for high_fluid_balance
coef_sens_cont <- summary_sens_cont$coefficients["high_fluid_balance", ]
beta_sens_cont <- coef_sens_cont["Estimate"]
se_sens_cont <- coef_sens_cont["Std. Error"]
t_sens_cont <- coef_sens_cont["t value"]
p_sens_cont <- coef_sens_cont["Pr(>|t|)"]

# Calculate 95% CI
ci_lower_sens_cont <- beta_sens_cont - 1.96 * se_sens_cont
ci_upper_sens_cont <- beta_sens_cont + 1.96 * se_sens_cont

cat("\nModel: iap_next ~ high_fluid_balance + factor(day)\n")
cat("Weights: Stabilized IPTW (truncated at 1st-99th percentile)\n")
cat("Family: Gaussian (identity link)\n\n")

cat("Results for high_fluid_balance:\n")
cat(sprintf("  Beta coefficient: %.4f mmHg\n", beta_sens_cont))
cat(sprintf("  Standard error: %.4f\n", se_sens_cont))
cat(sprintf("  t value: %.4f\n", t_sens_cont))
cat(sprintf("  P value: %.4f\n", p_sens_cont))
cat(sprintf("  Mean difference: %.3f mmHg\n", beta_sens_cont))
cat(sprintf("  95%% CI: [%.3f, %.3f] mmHg\n", ci_lower_sens_cont, ci_upper_sens_cont))
cat("\n")

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep="")
cat("SUMMARY: MSM CAUSAL EFFECT ESTIMATES\n")
cat(rep("=", 70), "\n\n", sep="")

cat("Exposure: High fluid balance (≥1010 mL) vs Non-high (<1010 mL)\n")
cat("Adjustment: Time trend (factor(day))\n")
cat("Confounding control: Stabilized IPTW\n\n")

cat("BINARY OUTCOME (IAH15_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-30s %10s %20s %10s\n", "Analysis", "OR", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-30s %10.3f %20s %10.4f\n",
            "Primary (original weights)",
            or_primary,
            sprintf("[%.3f, %.3f]", ci_lower_primary, ci_upper_primary),
            p_primary))
cat(sprintf("%-30s %10.3f %20s %10.4f\n",
            "Sensitivity (truncated)",
            or_sens,
            sprintf("[%.3f, %.3f]", ci_lower_sens, ci_upper_sens),
            p_sens))
cat("\n")

cat("CONTINUOUS OUTCOME (iap_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-30s %10s %20s %10s\n", "Analysis", "Beta (mmHg)", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-30s %10.3f %20s %10.4f\n",
            "Primary (original weights)",
            beta_primary_cont,
            sprintf("[%.3f, %.3f]", ci_lower_primary_cont, ci_upper_primary_cont),
            p_primary_cont))
cat(sprintf("%-30s %10.3f %20s %10.4f\n",
            "Sensitivity (truncated)",
            beta_sens_cont,
            sprintf("[%.3f, %.3f]", ci_lower_sens_cont, ci_upper_sens_cont),
            p_sens_cont))
cat("\n")

cat(rep("=", 70), "\n", sep="")
cat("INTERPRETATION:\n")
cat(rep("=", 70), "\n\n", sep="")

if (p_primary < 0.05) {
  cat("Binary outcome (IAH15):\n")
  cat(sprintf("  High fluid balance is associated with %.1f%% ", (or_primary - 1) * 100))
  if (or_primary > 1) {
    cat("INCREASED odds of IAH15\n")
  } else {
    cat("DECREASED odds of IAH15\n")
  }
  cat(sprintf("  (OR = %.3f, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              or_primary, ci_lower_primary, ci_upper_primary, p_primary))
} else {
  cat("Binary outcome (IAH15):\n")
  cat("  No statistically significant causal effect detected\n")
  cat(sprintf("  (OR = %.3f, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              or_primary, ci_lower_primary, ci_upper_primary, p_primary))
}

if (p_primary_cont < 0.05) {
  cat("Continuous outcome (IAP):\n")
  cat(sprintf("  High fluid balance causes a mean %.3f mmHg ", beta_primary_cont))
  if (beta_primary_cont > 0) {
    cat("INCREASE in IAP\n")
  } else {
    cat("DECREASE in IAP\n")
  }
  cat(sprintf("  (Beta = %.3f, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              beta_primary_cont, ci_lower_primary_cont, ci_upper_primary_cont, p_primary_cont))
} else {
  cat("Continuous outcome (IAP):\n")
  cat("  No statistically significant causal effect detected\n")
  cat(sprintf("  (Beta = %.3f, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              beta_primary_cont, ci_lower_primary_cont, ci_upper_primary_cont, p_primary_cont))
}

cat("Sensitivity analysis (truncated weights) shows ", sep="")
if (abs(or_sens - or_primary) / or_primary < 0.1 &&
    abs(beta_sens_cont - beta_primary_cont) / abs(beta_primary_cont) < 0.1) {
  cat("CONSISTENT results,\n")
  cat("indicating robustness to extreme weights.\n")
} else {
  cat("DIFFERENT results,\n")
  cat("suggesting sensitivity to extreme weights.\n")
}

cat("\n=== MSM STEP 3 COMPLETED ===\n")
