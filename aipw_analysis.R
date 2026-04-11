# AIPW (Augmented Inverse Probability Weighting) Analysis
# Date: 2026-03-30
# Task: Double-robust estimation for causal effects

library(dplyr)

# Set working directory
setwd(project_root)

# Read baseline IPTW weights data (contains all needed variables)
weights_data <- read.csv("data/IPTW_baseline_model_weights.csv")

# Read longitudinal data to get iap_next
library(readxl)
longitudinal_data <- read_excel("data/longitudinal_data_final.xlsx")

# Merge iap_next
data_aipw <- weights_data %>%
  left_join(
    longitudinal_data %>% select(subject_id, day, iap_next),
    by = c("subject_id", "day")
  )

cat("\n=== AIPW DOUBLE-ROBUST ANALYSIS ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Sample size:", nrow(data_aipw), "person-time observations\n")
cat("Exposure: high_fluid_balance (1 = ≥1010 mL, 0 = <1010 mL)\n")
cat("Outcomes: IAH15_next (binary), iap_next (continuous)\n\n")

# ============================================================================
# STEP 1: FIT TREATMENT MODEL (Propensity Score)
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 1: TREATMENT MODEL (Propensity Score)\n")
cat(rep("=", 70), "\n\n", sep="")

ps_model <- glm(
  high_fluid_balance ~ age + sex + etiology + iap_t + apache_t,
  data = data_aipw,
  family = binomial(link = "logit")
)

# Get propensity scores
data_aipw$ps <- predict(ps_model, type = "response")

cat("Treatment model fitted successfully\n")
cat("Covariates: age, sex, etiology, iap_t, apache_t\n")
cat("Propensity score range: [", round(min(data_aipw$ps), 4), ", ",
    round(max(data_aipw$ps), 4), "]\n", sep="")
cat("Mean PS in treated (high FB): ", round(mean(data_aipw$ps[data_aipw$high_fluid_balance == 1]), 4), "\n", sep="")
cat("Mean PS in control (non-high FB): ", round(mean(data_aipw$ps[data_aipw$high_fluid_balance == 0]), 4), "\n\n", sep="")

# ============================================================================
# STEP 2: FIT OUTCOME MODELS
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 2: OUTCOME MODELS\n")
cat(rep("=", 70), "\n\n", sep="")

# Binary outcome: IAH15_next
cat("Binary outcome model (IAH15_next):\n")
cat(rep("-", 70), "\n", sep="")

outcome_model_binary <- glm(
  IAH15_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t,
  data = data_aipw,
  family = binomial(link = "logit")
)

cat("Model: logistic regression\n")
cat("Covariates: high_fluid_balance, age, sex, etiology, iap_t, apache_t\n\n")

# Continuous outcome: iap_next
cat("Continuous outcome model (iap_next):\n")
cat(rep("-", 70), "\n", sep="")

outcome_model_cont <- lm(
  iap_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t,
  data = data_aipw
)

cat("Model: linear regression\n")
cat("Covariates: high_fluid_balance, age, sex, etiology, iap_t, apache_t\n\n")

# ============================================================================
# STEP 3: COMPUTE AIPW ESTIMATES
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("STEP 3: AIPW DOUBLE-ROBUST ESTIMATION\n")
cat(rep("=", 70), "\n\n", sep="")

# --- Binary Outcome: IAH15_next ---
cat("BINARY OUTCOME: IAH15_next\n")
cat(rep("-", 70), "\n", sep="")

# Predict potential outcomes under both treatment conditions
data_aipw_treated <- data_aipw
data_aipw_treated$high_fluid_balance <- 1

data_aipw_control <- data_aipw
data_aipw_control$high_fluid_balance <- 0

# Get predicted probabilities
mu1_binary <- predict(outcome_model_binary, newdata = data_aipw_treated, type = "response")
mu0_binary <- predict(outcome_model_binary, newdata = data_aipw_control, type = "response")

# AIPW estimator for binary outcome (risk difference)
# E[Y(1)] = mean(A/ps * Y + (1 - A/ps) * mu1)
# E[Y(0)] = mean((1-A)/(1-ps) * Y + (1 - (1-A)/(1-ps)) * mu0)

A <- data_aipw$high_fluid_balance
Y_binary <- data_aipw$IAH15_next
ps <- data_aipw$ps

# Potential outcome under treatment
po1_binary <- (A / ps) * Y_binary + (1 - A / ps) * mu1_binary

# Potential outcome under control
po0_binary <- ((1 - A) / (1 - ps)) * Y_binary + (1 - (1 - A) / (1 - ps)) * mu0_binary

# AIPW estimates
E_Y1_binary <- mean(po1_binary)
E_Y0_binary <- mean(po0_binary)

# Risk difference
rd_binary <- E_Y1_binary - E_Y0_binary

# Risk ratio
rr_binary <- E_Y1_binary / E_Y0_binary

# Odds ratio (from risk estimates)
odds1 <- E_Y1_binary / (1 - E_Y1_binary)
odds0 <- E_Y0_binary / (1 - E_Y0_binary)
or_binary <- odds1 / odds0

# Bootstrap for confidence intervals (binary outcome)
set.seed(123)
n_boot <- 1000
boot_rd <- numeric(n_boot)
boot_rr <- numeric(n_boot)
boot_or <- numeric(n_boot)

cat("Computing bootstrap confidence intervals (1000 iterations)...\n")

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(data_aipw), replace = TRUE)
  boot_data <- data_aipw[boot_idx, ]

  # Refit models
  ps_boot <- glm(high_fluid_balance ~ age + sex + etiology + iap_t + apache_t,
                 data = boot_data, family = binomial)$fitted.values

  outcome_boot <- glm(IAH15_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t,
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

  boot_rd[i] <- E_Y1_boot - E_Y0_boot
  boot_rr[i] <- E_Y1_boot / E_Y0_boot

  odds1_boot <- E_Y1_boot / (1 - E_Y1_boot)
  odds0_boot <- E_Y0_boot / (1 - E_Y0_boot)
  boot_or[i] <- odds1_boot / odds0_boot
}

# Calculate 95% CI
rd_ci <- quantile(boot_rd, c(0.025, 0.975))
rr_ci <- quantile(boot_rr, c(0.025, 0.975))
or_ci <- quantile(boot_or, c(0.025, 0.975))

# P-value (two-sided test for RD)
p_value_binary <- 2 * min(mean(boot_rd >= 0), mean(boot_rd <= 0))

cat("\nAIPW Results:\n")
cat(sprintf("  E[Y(1)] (risk under high FB): %.4f\n", E_Y1_binary))
cat(sprintf("  E[Y(0)] (risk under non-high FB): %.4f\n", E_Y0_binary))
cat(sprintf("\n  Risk Difference: %.4f (95%% CI: [%.4f, %.4f])\n",
            rd_binary, rd_ci[1], rd_ci[2]))
cat(sprintf("  Risk Ratio: %.4f (95%% CI: [%.4f, %.4f])\n",
            rr_binary, rr_ci[1], rr_ci[2]))
cat(sprintf("  Odds Ratio: %.4f (95%% CI: [%.4f, %.4f])\n",
            or_binary, or_ci[1], or_ci[2]))
cat(sprintf("  P-value: %.4f\n\n", p_value_binary))

# --- Continuous Outcome: iap_next ---
cat("CONTINUOUS OUTCOME: iap_next\n")
cat(rep("-", 70), "\n", sep="")

# Get predicted outcomes
mu1_cont <- predict(outcome_model_cont, newdata = data_aipw_treated)
mu0_cont <- predict(outcome_model_cont, newdata = data_aipw_control)

# AIPW estimator for continuous outcome
Y_cont <- data_aipw$iap_next

# Potential outcome under treatment
po1_cont <- (A / ps) * Y_cont + (1 - A / ps) * mu1_cont

# Potential outcome under control
po0_cont <- ((1 - A) / (1 - ps)) * Y_cont + (1 - (1 - A) / (1 - ps)) * mu0_cont

# AIPW estimates
E_Y1_cont <- mean(po1_cont)
E_Y0_cont <- mean(po0_cont)

# Average treatment effect
ate_cont <- E_Y1_cont - E_Y0_cont

# Bootstrap for confidence intervals (continuous outcome)
boot_ate <- numeric(n_boot)

cat("Computing bootstrap confidence intervals (1000 iterations)...\n")

for (i in 1:n_boot) {
  boot_idx <- sample(1:nrow(data_aipw), replace = TRUE)
  boot_data <- data_aipw[boot_idx, ]

  # Refit models
  ps_boot <- glm(high_fluid_balance ~ age + sex + etiology + iap_t + apache_t,
                 data = boot_data, family = binomial)$fitted.values

  outcome_boot <- lm(iap_next ~ high_fluid_balance + age + sex + etiology + iap_t + apache_t,
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

  boot_ate[i] <- mean(po1_boot) - mean(po0_boot)
}

# Calculate 95% CI
ate_ci <- quantile(boot_ate, c(0.025, 0.975))

# P-value (two-sided test)
p_value_cont <- 2 * min(mean(boot_ate >= 0), mean(boot_ate <= 0))

cat("\nAIPW Results:\n")
cat(sprintf("  E[Y(1)] (IAP under high FB): %.4f mmHg\n", E_Y1_cont))
cat(sprintf("  E[Y(0)] (IAP under non-high FB): %.4f mmHg\n", E_Y0_cont))
cat(sprintf("\n  Average Treatment Effect: %.4f mmHg (95%% CI: [%.4f, %.4f])\n",
            ate_cont, ate_ci[1], ate_ci[2]))
cat(sprintf("  P-value: %.4f\n\n", p_value_cont))

# ============================================================================
# SUMMARY
# ============================================================================

cat(rep("=", 70), "\n", sep="")
cat("SUMMARY: AIPW DOUBLE-ROBUST ESTIMATES\n")
cat(rep("=", 70), "\n\n", sep="")

cat("Exposure: High fluid balance (≥1010 mL) vs Non-high (<1010 mL)\n")
cat("Method: Augmented Inverse Probability Weighting (AIPW)\n")
cat("Treatment model: age, sex, etiology, iap_t, apache_t\n")
cat("Outcome model: high_fluid_balance + age, sex, etiology, iap_t, apache_t\n\n")

cat("BINARY OUTCOME (IAH15_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10s %25s %10s\n", "Measure", "Estimate", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10.4f %25s %10.4f\n",
            "Risk Difference",
            rd_binary,
            sprintf("[%.4f, %.4f]", rd_ci[1], rd_ci[2]),
            p_value_binary))
cat(sprintf("%-25s %10.4f %25s %10s\n",
            "Risk Ratio",
            rr_binary,
            sprintf("[%.4f, %.4f]", rr_ci[1], rr_ci[2]),
            "-"))
cat(sprintf("%-25s %10.4f %25s %10s\n",
            "Odds Ratio",
            or_binary,
            sprintf("[%.4f, %.4f]", or_ci[1], or_ci[2]),
            "-"))
cat("\n")

cat("CONTINUOUS OUTCOME (iap_next):\n")
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10s %25s %10s\n", "Measure", "Estimate", "95% CI", "P value"))
cat(rep("-", 70), "\n", sep="")
cat(sprintf("%-25s %10.4f %25s %10.4f\n",
            "ATE (mmHg)",
            ate_cont,
            sprintf("[%.4f, %.4f]", ate_ci[1], ate_ci[2]),
            p_value_cont))
cat("\n")

cat(rep("=", 70), "\n", sep="")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep="")

if (p_value_binary < 0.05) {
  cat("Binary outcome (IAH15):\n")
  cat(sprintf("  High fluid balance increases the risk of IAH15 by %.1f percentage points\n",
              rd_binary * 100))
  cat(sprintf("  (Risk Difference = %.4f, 95%% CI [%.4f, %.4f], P = %.4f)\n",
              rd_binary, rd_ci[1], rd_ci[2], p_value_binary))
  cat(sprintf("  Odds Ratio = %.3f (95%% CI [%.3f, %.3f])\n\n",
              or_binary, or_ci[1], or_ci[2]))
} else {
  cat("Binary outcome (IAH15):\n")
  cat("  No statistically significant effect detected\n\n")
}

if (p_value_cont < 0.05) {
  cat("Continuous outcome (IAP):\n")
  cat(sprintf("  High fluid balance increases IAP by %.3f mmHg on average\n", ate_cont))
  cat(sprintf("  (ATE = %.3f mmHg, 95%% CI [%.3f, %.3f], P = %.4f)\n\n",
              ate_cont, ate_ci[1], ate_ci[2], p_value_cont))
} else {
  cat("Continuous outcome (IAP):\n")
  cat("  No statistically significant effect detected\n\n")
}

cat("Double-robust property: AIPW estimates are consistent if EITHER\n")
cat("the treatment model OR the outcome model is correctly specified.\n")

cat("\n=== AIPW ANALYSIS COMPLETED ===\n")
