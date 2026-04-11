# AIPW Secondary Analysis (Extended Model with cr_t and map_t)
# Date: 2026-03-30
# Task: Double-robust estimation with additional confounders

library(dplyr)

# Set working directory
setwd(project_root)

# Read longitudinal data
library(readxl)
longitudinal_data <- read_excel("data/longitudinal_data_final.xlsx")

# Prepare data with extended covariates
data_aipw_ext <- longitudinal_data %>%
  filter(day %in% c(1, 2, 3, 5)) %>%
  mutate(
    high_fluid_balance = ifelse(fluid_balance_t >= 1010, 1, 0),
    IAH15_next = ifelse(iap_next >= 15, 1, 0)
  ) %>%
  select(subject_id, day, high_fluid_balance, fluid_balance_t,
         age, sex, etiology, iap_t, apache_t, cr_t, map_t,
         IAH15_next, iap_next) %>%
  filter(complete.cases(.))

cat("\n=== AIPW SECONDARY ANALYSIS (EXTENDED MODEL) ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Sample size:", nrow(data_aipw_ext), "person-time observations\n")
cat("Note: Reduced sample due to additional covariates (cr_t, map_t)\n")
cat("Exposure: high_fluid_balance (1 = ≥1010 mL, 0 = <1010 mL)\n")
cat("Outcomes: IAH15_next (binary), iap_next (continuous)\n\n")

# ============================================================================
# STEP 1: FIT TREATMENT MODEL (Propensity Score) - EXTENDED
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 1: TREATMENT MODEL (Propensity Score) - EXTENDED\n")
cat(rep("=", 70), "\n\n", sep="")

ps_model_ext <- glm(
  high_fluid_balance ~ age + sex + etiology + iap_t + apache_t + cr_t + map_t,
  data = data_aipw_ext,
  family = binomial(link = "logit")
)

# Get propensity scores
data_aipw_ext$ps <- predict(ps_model_ext, type = "response")

cat("Treatment model fitted successfully\n")
cat("Covariates: age, sex, etiology, iap_t, apache_t, cr_t, map_t\n")
cat("Propensity score range: [", round(min(data_aipw_ext$ps), 4), ", ",
    round(max(data_aipw_ext$ps), 4), "]\n", sep="")
cat("Mean PS in treated (high FB): ", round(mean(data_aipw_ext$ps[data_aipw_ext$high_fluid_balance == 1]), 4), "\n", sep="")
cat("Mean PS in control (non-high FB): ", round(mean(data_aipw_ext$ps[data_aipw_ext$high_fluid_balance == 0]), 4), "\n\n", sep="")

# ============================================================================
# STEP 2: FIT OUTCOME MODELS - EXTENDED
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 2: OUTCOME MODELS - EXTENDED\n")
cat(rep("=", 70), "\n\n", sep="")

# Binary outcome: IAH15_next
cat("Binary outcome model (IAH15_next):\n")
cat(rep("-", 70), "\n", sep="")

outcome_model_binary_ext <- glm(
  IAH15_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t + cr_t + map_t,
  data = data_aipw_ext,
  family = binomial(link = "logit")
)

cat("Model: logistic regression\n")
cat("Covariates: high_fluid_balance, age, sex, etiology, iap_t, apache_t, cr_t, map_t\n\n")

# Continuous outcome: iap_next
cat("Continuous outcome model (iap_next):\n")
cat(rep("-", 70), "\n", sep="")

outcome_model_cont_ext <- lm(
  iap_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t + cr_t + map_t,
  data = data_aipw_ext
)

cat("Model: linear regression\n")
cat("Covariates: high_fluid_balance, age, sex, etiology, iap_t, apache_t, cr_t, map_t\n\n")

# ============================================================================
# STEP 3: COMPUTE AIPW ESTIMATES - EXTENDED
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 3: AIPW DOUBLE-ROBUST ESTIMATION - EXTENDED\n")
cat(rep("=", 70), "\n\n", sep="")

# --- Binary Outcome: IAH15_next ---
cat("BINARY OUTCOME: IAH15_next\n")
cat(rep("-", 70), "\n", sep="")

# Predict potential outcomes under both treatment conditions
data_aipw_ext_treated <- data_aipw_ext
data_aipw_ext_treated$high_fluid_balance <- 1

data_aipw_ext_control <- data_aipw_ext
data_aipw_ext_control$high_fluid_balance <- 0

# Get predicted probabilities
mu1_binary_ext <- predict(outcome_model_binary_ext, newdata = data_aipw_ext_treated, type = "response")
mu0_binary_ext <- predict(outcome_model_binary_ext, newdata = data_aipw_ext_control, type = "response")

# AIPW estimator for binary outcome
A_ext <- data_aipw_ext$high_fluid_balance
Y_binary_ext <- data_aipw_ext$IAH15_next
ps_ext <- data_aipw_ext$ps

# Potential outcome under treatment
po1_binary_ext <- (A_ext / ps_ext) * Y_binary_ext + (1 - A_ext / ps_ext) * mu1_binary_ext

# Potential outcome under control
po0_binary_ext <- ((1 - A_ext) / (1 - ps_ext)) * Y_binary_ext + (1 - (1 - A_ext) / (1 - ps_ext)) * mu0_binary_ext

# AIPW estimates
E_Y1_binary_ext <- mean(po1_binary_ext)
E_Y0_binary_ext <- mean(po0_binary_ext)

# Risk difference
rd_binary_ext <- E_Y1_binary_ext - E_Y0_binary_ext

# Risk ratio
rr_binary_ext <- E_Y1_binary_ext / E_Y0_binary_ext

# Odds ratio
odds1_ext <- E_Y1_binary_ext / (1 - E_Y1_binary_ext)
odds0_ext <- E_Y0_binary_ext / (1 - E_Y0_binary_ext)
or_binary_ext <- odds1_ext / odds0_ext

# Bootstrap for confidence intervals (binary outcome)
set.seed(123)
n_boot <- 1000
boot_rd_ext <- numeric(n_boot)
boot_rr_ext <- numeric(n_boot)
boot_or_ext <- numeric(n_boot)

cat("Computing bootstrap confidence intervals (1000 iterations)...\n")

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(data_aipw_ext), replace = TRUE)
  boot_data <- data_aipw_ext[boot_idx, ]

  # Refit models
  ps_boot <- glm(high_fluid_balance ~ age + sex + etiology + iap_t + apache_t + cr_t + map_t,
                 data = boot_data, family = binomial)$fitted.values

  outcome_boot <- glm(IAH15_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t + cr_t + map_t,
                      data = boot_data, family = binomial)

  boot_data_treated <- boot_data
  boot_data_treated$high_fluid_balance <- 1
  boot_data_control <- boot_data
  boot_data_control$high_fluid_balance <- 0

  mu1_boot <- predict(outcome_boot, newdata = boot_data_treated, type = "response")
  mu0_boot <- predict(outcome_boot, newdata = boot_data_control, type = "response")

  A_boot <- boot_data$high_fluid_balance
  Y_boot <- boot_data$IAH15_next

  po1_boot <- (A_boot / ps_boot) * Y_boot + (1 - A_boot / ps_boot) * mu1_boot
  po0_boot <- ((1 - A_boot) / (1 - ps_boot)) * Y_boot + (1 - (1 - A_boot) / (1 - ps_boot)) * mu0_boot

  E_Y1_boot <- mean(po1_boot)
  E_Y0_boot <- mean(po0_boot)

  boot_rd_ext[i] <- E_Y1_boot - E_Y0_boot
  boot_rr_ext[i] <- E_Y1_boot / E_Y0_boot

  odds1_boot <- E_Y1_boot / (1 - E_Y1_boot)
  odds0_boot <- E_Y0_boot / (1 - E_Y0_boot)
  boot_or_ext[i] <- odds1_boot / odds0_boot
}

# Calculate 95% CI
rd_ci_ext <- quantile(boot_rd_ext, c(0.025, 0.975))
rr_ci_ext <- quantile(boot_rr_ext, c(0.025, 0.975))
or_ci_ext <- quantile(boot_or_ext, c(0.025, 0.975))

# P-value
p_value_binary_ext <- 2 * min(mean(boot_rd_ext >= 0), mean(boot_rd_ext <= 0))

cat("\nAIPW Results (Extended Model):\n")
cat(sprintf("  E[Y(1)] (risk under high FB): %.4f\n", E_Y1_binary_ext))
cat(sprintf("  E[Y(0)] (risk under non-high FB): %.4f\n", E_Y0_binary_ext))
cat(sprintf("\n  Risk Difference: %.4f (95%% CI: [%.4f, %.4f])\n",
            rd_binary_ext, rd_ci_ext[1], rd_ci_ext[2]))
cat(sprintf("  Risk Ratio: %.4f (95%% CI: [%.4f, %.4f])\n",
            rr_binary_ext, rr_ci_ext[1], rr_ci_ext[2]))
cat(sprintf("  Odds Ratio: %.4f (95%% CI: [%.4f, %.4f])\n",
            or_binary_ext, or_ci_ext[1], or_ci_ext[2]))
cat(sprintf("  P-value: %.4f\n\n", p_value_binary_ext))

# --- Continuous Outcome: iap_next ---
cat("CONTINUOUS OUTCOME: iap_next\n")
cat(rep("-", 70), "\n", sep="")

# Get predicted outcomes
mu1_cont_ext <- predict(outcome_model_cont_ext, newdata = data_aipw_ext_treated)
mu0_cont_ext <- predict(outcome_model_cont_ext, newdata = data_aipw_ext_control)

# AIPW estimator for continuous outcome
Y_cont_ext <- data_aipw_ext$iap_next

# Potential outcome under treatment
po1_cont_ext <- (A_ext / ps_ext) * Y_cont_ext + (1 - A_ext / ps_ext) * mu1_cont_ext

# Potential outcome under control
po0_cont_ext <- ((1 - A_ext) / (1 - ps_ext)) * Y_cont_ext + (1 - (1 - A_ext) / (1 - ps_ext)) * mu0_cont_ext

# AIPW estimates
E_Y1_cont_ext <- mean(po1_cont_ext)
E_Y0_cont_ext <- mean(po0_cont_ext)

# Average treatment effect
ate_cont_ext <- E_Y1_cont_ext - E_Y0_cont_ext

# Bootstrap for confidence intervals (continuous outcome)
boot_ate_ext <- numeric(n_boot)

cat("Computing bootstrap confidence intervals (1000 iterations)...\n")

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(data_aipw_ext), replace = TRUE)
  boot_data <- data_aipw_ext[boot_idx, ]

  # Refit models
  ps_boot <- glm(high_fluid_balance ~ age + sex + etiology + iap_t + apache_t + cr_t + map_t,
                 data = boot_data, family = binomial)$fitted.values

  outcome_boot <- lm(iap_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t + cr_t + map_t,
                     data = boot_data)

  boot_data_treated <- boot_data
  boot_data_treated$high_fluid_balance <- 1
  boot_data_control <- boot_data
  boot_data_control$high_fluid_balance <- 0

  mu1_boot <- predict(outcome_boot, newdata = boot_data_treated)
  mu0_boot <- predict(outcome_boot, newdata = boot_data_control)

  A_boot <- boot_data$high_fluid_balance
  Y_boot <- boot_data$iap_next

  po1_boot <- (A_boot / ps_boot) * Y_boot + (1 - A_boot / ps_boot) * mu1_boot
  po0_boot <- ((1 - A_boot) / (1 - ps_boot)) * Y_boot + (1 - (1 - A_boot) / (1 - ps_boot)) * mu0_boot

  boot_ate_ext[i] <- mean(po1_boot) - mean(po0_boot)
}

# Calculate 95% CI
ate_ci_ext <- quantile(boot_ate_ext, c(0.025, 0.975))

# P-value
p_value_cont_ext <- 2 * min(mean(boot_ate_ext >= 0), mean(boot_ate_ext <= 0))

cat("\nAIPW Results (Extended Model):\n")
cat(sprintf("  E[Y(1)] (IAP under high FB): %.4f mmHg\n", E_Y1_cont_ext))
cat(sprintf("  E[Y(0)] (IAP under non-high FB): %.4f mmHg\n", E_Y0_cont_ext))
cat(sprintf("\n  Average Treatment Effect: %.4f mmHg (95%% CI: [%.4f, %.4f])\n",
            ate_cont_ext, ate_ci_ext[1], ate_ci_ext[2]))
cat(sprintf("  P-value: %.4f\n\n", p_value_cont_ext))

# ============================================================================
# SUMMARY
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("SUMMARY: AIPW SECONDARY ANALYSIS (EXTENDED MODEL)\n")
cat(rep("=", 70), "\n\n", sep="")

cat("Exposure: High fluid balance (≥1010 mL) vs Non-high (<1010 mL)\n")
cat("Method: Augmented Inverse Probability Weighting (AIPW)\n")
cat("Treatment model: age, sex, etiology, iap_t, apache_t, cr_t, map_t\n")
cat("Outcome model: high_fluid_balance + age, sex, etiology, iap_t, apache_t, cr_t, map_t\n")
cat("Sample size: N =", nrow(data_aipw_ext), "(reduced due to additional covariates)\n\n")

cat("BINARY OUTCOME (IAH15_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10s %25s %10s\n", "Measure", "Estimate", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10.4f %25s %10.4f\n",
            "Risk Difference",
            rd_binary_ext,
            sprintf("[%.4f, %.4f]", rd_ci_ext[1], rd_ci_ext[2]),
            p_value_binary_ext))
cat(sprintf("%-25s %10.4f %25s %10s\n",
            "Risk Ratio",
            rr_binary_ext,
            sprintf("[%.4f, %.4f]", rr_ci_ext[1], rr_ci_ext[2]),
            "-"))
cat(sprintf("%-25s %10.4f %25s %10s\n",
            "Odds Ratio",
            or_binary_ext,
            sprintf("[%.4f, %.4f]", or_ci_ext[1], or_ci_ext[2]),
            "-"))
cat("\n")

cat("CONTINUOUS OUTCOME (iap_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10s %25s %10s\n", "Measure", "Estimate", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10.4f %25s %10.4f\n",
            "ATE (mmHg)",
            ate_cont_ext,
            sprintf("[%.4f, %.4f]", ate_ci_ext[1], ate_ci_ext[2]),
            p_value_cont_ext))
cat("\n")

cat(rep("=", 70), "\n", sep="")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep="")

if (p_value_binary_ext < 0.05) {
  cat("Binary outcome (IAH15):\n")
  cat(sprintf("  High fluid balance increases the risk of IAH15 by %.1f percentage points\n",
              rd_binary_ext * 100))
  cat(sprintf("  (Risk Difference = %.4f, 95%% CI [%.4f, %.4f], P = %.4f)\n",
              rd_binary_ext, rd_ci_ext[1], rd_ci_ext[2], p_value_binary_ext))
  cat(sprintf("  Odds Ratio = %.3f (95%% CI [%.3f, %.3f])\n\n",
              or_binary_ext, or_ci_ext[1], or_ci_ext[2]))
} else {
  cat("Binary outcome (IAH15):\n")
  cat("  No statistically significant effect detected\n\n")
}

if (p_value_cont_ext < 0.05) {
  cat("Continuous outcome (IAP):\n")
  cat(sprintf("  High fluid balance increases IAP by %.3f mmHg on average\n", ate_cont_ext))
  cat(sprintf("  (ATE = %.3f mmHg, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              ate_cont_ext, ate_ci_ext[1], ate_ci_ext[2], p_value_cont_ext))
} else {
  cat("Continuous outcome (IAP):\n")
  cat("  No statistically significant effect detected\n\n")
}

cat("Note: Extended model includes additional confounders (cr_t, map_t)\n")
cat("for more comprehensive adjustment, at the cost of reduced sample size.\n")

cat("\n=== AIPW SECONDARY ANALYSIS COMPLETED ===\n")
