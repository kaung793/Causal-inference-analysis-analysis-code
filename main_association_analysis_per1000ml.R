start_time <- Sys.time()
project_root <- "PROJECT_ROOT_PLACEHOLDER"
setwd(project_root)

required_packages <- c("dplyr", "readxl", "lme4", "lmerTest", "officer", "flextable", "ggplot2", "gridExtra")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(dplyr)
library(readxl)
library(lme4)
library(lmerTest)
library(officer)
library(flextable)
library(ggplot2)
library(gridExtra)

ts <- format(start_time, "%Y%m%d_%H%M")
out_dir <- file.path("output", ts)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

warnings_vec <- character()
errors_vec <- character()
convergence_notes <- character()

withCallingHandlers({

  # Read data - use imputed version for consistency
  data <- readxl::read_xlsx("data/longitudinal_data_final_imputed.xlsx")

  # Create binary outcomes
  data <- data %>%
    dplyr::mutate(
      IAH15_next = ifelse(iap_next >= 15, 1, 0),
      IAH20_next = ifelse(iap_next >= 20, 1, 0)
    )

  # Calculate SD for exposures
  sd_fluid_t <- sd(data$fluid_t, na.rm = TRUE)
  sd_fluid_balance_t <- sd(data$fluid_balance_t, na.rm = TRUE)

  # Initialize results storage
  results_list <- list()

  # Function to fit linear mixed model with fallback
  fit_lmm <- function(formula_str, data, model_name) {
    tryCatch({
      model <- lmer(as.formula(formula_str), data = data, REML = FALSE)
      convergence_notes <<- c(convergence_notes,
                              paste0(model_name, ": Mixed model converged"))
      return(list(model = model, type = "mixed", converged = TRUE))
    }, error = function(e) {
      convergence_notes <<- c(convergence_notes,
                              paste0(model_name, ": Mixed model failed, using LM. Reason: ", e$message))
      # Fallback to regular linear model
      formula_lm <- gsub("\\s*\\+\\s*\\(1\\s*\\|\\s*subject_id\\)", "", formula_str)
      model <- lm(as.formula(formula_lm), data = data)
      return(list(model = model, type = "lm", converged = FALSE))
    })
  }

  # Function to fit logistic mixed model with fallback
  fit_glmm <- function(formula_str, data, model_name) {
    tryCatch({
      model <- glmer(as.formula(formula_str), data = data, family = binomial,
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
      convergence_notes <<- c(convergence_notes,
                              paste0(model_name, ": Mixed logistic model converged"))
      return(list(model = model, type = "mixed", converged = TRUE))
    }, error = function(e) {
      convergence_notes <<- c(convergence_notes,
                              paste0(model_name, ": Mixed logistic model failed, using GLM. Reason: ", e$message))
      # Fallback to regular logistic model
      formula_glm <- gsub("\\s*\\+\\s*\\(1\\s*\\|\\s*subject_id\\)", "", formula_str)
      model <- glm(as.formula(formula_glm), data = data, family = binomial)
      return(list(model = model, type = "glm", converged = FALSE))
    })
  }

  # Function to extract results from continuous models
  extract_continuous_results <- function(model_obj, exposure_var, outcome_var, model_label,
                                        sd_exposure, exposure_name) {
    model <- model_obj$model
    model_type <- model_obj$type

    if (model_type == "mixed") {
      coef_summary <- summary(model)$coefficients
    } else {
      coef_summary <- summary(model)$coefficients
    }

    # Extract coefficient for exposure
    beta <- coef_summary[exposure_var, "Estimate"]
    se <- coef_summary[exposure_var, "Std. Error"]
    p_val <- coef_summary[exposure_var, ifelse(model_type == "mixed", "Pr(>|t|)", "Pr(>|t|)")]

    # Calculate CI
    ci_lower <- beta - 1.96 * se
    ci_upper <- beta + 1.96 * se

    # Per SD effect
    beta_per_sd <- beta * sd_exposure
    ci_lower_per_sd <- ci_lower * sd_exposure
    ci_upper_per_sd <- ci_upper * sd_exposure

    # Per 1000 mL effect
    beta_per_1000ml <- beta * 1000
    ci_lower_per_1000ml <- ci_lower * 1000
    ci_upper_per_1000ml <- ci_upper * 1000

    return(data.frame(
      Exposure = exposure_name,
      Outcome = "iap_next",
      Outcome_Label = "IAP (continuous)",
      Model = model_label,
      Effect_per_1000mL = beta_per_1000ml,
      CI1000_Lower = ci_lower_per_1000ml,
      CI1000_Upper = ci_upper_per_1000ml,
      Effect_per_SD = beta_per_sd,
      CISD_Lower = ci_lower_per_sd,
      CISD_Upper = ci_upper_per_sd,
      P_value = p_val,
      Model_type = model_type,
      N = nobs(model),
      stringsAsFactors = FALSE
    ))
  }

  # Function to extract results from binary models
  extract_binary_results <- function(model_obj, exposure_var, outcome_var, model_label,
                                    sd_exposure, exposure_name, outcome_name) {
    model <- model_obj$model
    model_type <- model_obj$type

    if (model_type == "mixed") {
      coef_summary <- summary(model)$coefficients
    } else {
      coef_summary <- summary(model)$coefficients
    }

    # Extract coefficient for exposure
    log_or <- coef_summary[exposure_var, "Estimate"]
    se <- coef_summary[exposure_var, "Std. Error"]
    p_val <- coef_summary[exposure_var, ifelse(model_type == "mixed", "Pr(>|z|)", "Pr(>|z|)")]

    # Per 1000 mL effect
    or_per_1000ml <- exp(log_or * 1000)
    ci_lower_per_1000ml <- exp((log_or - 1.96 * se) * 1000)
    ci_upper_per_1000ml <- exp((log_or + 1.96 * se) * 1000)

    # Per SD effect
    or_per_sd <- exp(log_or * sd_exposure)
    ci_lower_per_sd <- exp((log_or - 1.96 * se) * sd_exposure)
    ci_upper_per_sd <- exp((log_or + 1.96 * se) * sd_exposure)

    return(data.frame(
      Exposure = exposure_name,
      Outcome = outcome_var,
      Outcome_Label = outcome_name,
      Model = model_label,
      Effect_per_1000mL = or_per_1000ml,
      CI1000_Lower = ci_lower_per_1000ml,
      CI1000_Upper = ci_upper_per_1000ml,
      Effect_per_SD = or_per_sd,
      CISD_Lower = ci_lower_per_sd,
      CISD_Upper = ci_upper_per_sd,
      P_value = p_val,
      Model_type = model_type,
      N = nobs(model),
      stringsAsFactors = FALSE
    ))
  }

  # ========== Analysis for fluid_t ==========

  # Continuous outcome: iap_next

  # Model 1: Crude
  cat("Fitting fluid_t -> iap_next, Model 1 (Crude)...\n")
  m1_fluid_iap <- fit_lmm("iap_next ~ fluid_t + (1|subject_id)", data, "fluid_t->iap_next Model1")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m1_fluid_iap, "fluid_t", "iap_next", "Model 1", sd_fluid_t, "Fluid intake"
  )

  # Model 2: Main
  cat("Fitting fluid_t -> iap_next, Model 2 (Main)...\n")
  m2_fluid_iap <- fit_lmm("iap_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                          data, "fluid_t->iap_next Model2")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m2_fluid_iap, "fluid_t", "iap_next", "Model 2", sd_fluid_t, "Fluid intake"
  )

  # Model 3: Extended
  cat("Fitting fluid_t -> iap_next, Model 3 (Extended)...\n")
  m3_fluid_iap <- fit_lmm("iap_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                          data, "fluid_t->iap_next Model3")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m3_fluid_iap, "fluid_t", "iap_next", "Model 3", sd_fluid_t, "Fluid intake"
  )

  # Binary outcome: IAH15_next

  # Model 1: Crude
  cat("Fitting fluid_t -> IAH15_next, Model 1 (Crude)...\n")
  m1_fluid_iah15 <- fit_glmm("IAH15_next ~ fluid_t + (1|subject_id)", data, "fluid_t->IAH15 Model1")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m1_fluid_iah15, "fluid_t", "IAH15_next", "Model 1", sd_fluid_t, "Fluid intake", "IAH15"
  )

  # Model 2: Main
  cat("Fitting fluid_t -> IAH15_next, Model 2 (Main)...\n")
  m2_fluid_iah15 <- fit_glmm("IAH15_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                             data, "fluid_t->IAH15 Model2")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m2_fluid_iah15, "fluid_t", "IAH15_next", "Model 2", sd_fluid_t, "Fluid intake", "IAH15"
  )

  # Model 3: Extended
  cat("Fitting fluid_t -> IAH15_next, Model 3 (Extended)...\n")
  m3_fluid_iah15 <- fit_glmm("IAH15_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                             data, "fluid_t->IAH15 Model3")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m3_fluid_iah15, "fluid_t", "IAH15_next", "Model 3", sd_fluid_t, "Fluid intake", "IAH15"
  )

  # Binary outcome: IAH20_next

  # Model 1: Crude
  cat("Fitting fluid_t -> IAH20_next, Model 1 (Crude)...\n")
  m1_fluid_iah20 <- fit_glmm("IAH20_next ~ fluid_t + (1|subject_id)", data, "fluid_t->IAH20 Model1")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m1_fluid_iah20, "fluid_t", "IAH20_next", "Model 1", sd_fluid_t, "Fluid intake", "IAH20"
  )

  # Model 2: Main
  cat("Fitting fluid_t -> IAH20_next, Model 2 (Main)...\n")
  m2_fluid_iah20 <- fit_glmm("IAH20_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                             data, "fluid_t->IAH20 Model2")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m2_fluid_iah20, "fluid_t", "IAH20_next", "Model 2", sd_fluid_t, "Fluid intake", "IAH20"
  )

  # Model 3: Extended
  cat("Fitting fluid_t -> IAH20_next, Model 3 (Extended)...\n")
  m3_fluid_iah20 <- fit_glmm("IAH20_next ~ fluid_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                             data, "fluid_t->IAH20 Model3")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m3_fluid_iah20, "fluid_t", "IAH20_next", "Model 3", sd_fluid_t, "Fluid intake", "IAH20"
  )

  # ========== Analysis for fluid_balance_t ==========

  # Continuous outcome: iap_next

  # Model 1: Crude
  cat("Fitting fluid_balance_t -> iap_next, Model 1 (Crude)...\n")
  m1_balance_iap <- fit_lmm("iap_next ~ fluid_balance_t + (1|subject_id)", data, "fluid_balance_t->iap_next Model1")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m1_balance_iap, "fluid_balance_t", "iap_next", "Model 1", sd_fluid_balance_t, "Fluid balance"
  )

  # Model 2: Main
  cat("Fitting fluid_balance_t -> iap_next, Model 2 (Main)...\n")
  m2_balance_iap <- fit_lmm("iap_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                            data, "fluid_balance_t->iap_next Model2")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m2_balance_iap, "fluid_balance_t", "iap_next", "Model 2", sd_fluid_balance_t, "Fluid balance"
  )

  # Model 3: Extended
  cat("Fitting fluid_balance_t -> iap_next, Model 3 (Extended)...\n")
  m3_balance_iap <- fit_lmm("iap_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                            data, "fluid_balance_t->iap_next Model3")
  results_list[[length(results_list) + 1]] <- extract_continuous_results(
    m3_balance_iap, "fluid_balance_t", "iap_next", "Model 3", sd_fluid_balance_t, "Fluid balance"
  )

  # Binary outcome: IAH15_next

  # Model 1: Crude
  cat("Fitting fluid_balance_t -> IAH15_next, Model 1 (Crude)...\n")
  m1_balance_iah15 <- fit_glmm("IAH15_next ~ fluid_balance_t + (1|subject_id)", data, "fluid_balance_t->IAH15 Model1")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m1_balance_iah15, "fluid_balance_t", "IAH15_next", "Model 1", sd_fluid_balance_t, "Fluid balance", "IAH15"
  )

  # Model 2: Main
  cat("Fitting fluid_balance_t -> IAH15_next, Model 2 (Main)...\n")
  m2_balance_iah15 <- fit_glmm("IAH15_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                               data, "fluid_balance_t->IAH15 Model2")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m2_balance_iah15, "fluid_balance_t", "IAH15_next", "Model 2", sd_fluid_balance_t, "Fluid balance", "IAH15"
  )

  # Model 3: Extended
  cat("Fitting fluid_balance_t -> IAH15_next, Model 3 (Extended)...\n")
  m3_balance_iah15 <- fit_glmm("IAH15_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                               data, "fluid_balance_t->IAH15 Model3")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m3_balance_iah15, "fluid_balance_t", "IAH15_next", "Model 3", sd_fluid_balance_t, "Fluid balance", "IAH15"
  )

  # Binary outcome: IAH20_next

  # Model 1: Crude
  cat("Fitting fluid_balance_t -> IAH20_next, Model 1 (Crude)...\n")
  m1_balance_iah20 <- fit_glmm("IAH20_next ~ fluid_balance_t + (1|subject_id)", data, "fluid_balance_t->IAH20 Model1")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m1_balance_iah20, "fluid_balance_t", "IAH20_next", "Model 1", sd_fluid_balance_t, "Fluid balance", "IAH20"
  )

  # Model 2: Main
  cat("Fitting fluid_balance_t -> IAH20_next, Model 2 (Main)...\n")
  m2_balance_iah20 <- fit_glmm("IAH20_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + (1|subject_id)",
                               data, "fluid_balance_t->IAH20 Model2")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m2_balance_iah20, "fluid_balance_t", "IAH20_next", "Model 2", sd_fluid_balance_t, "Fluid balance", "IAH20"
  )

  # Model 3: Extended
  cat("Fitting fluid_balance_t -> IAH20_next, Model 3 (Extended)...\n")
  m3_balance_iah20 <- fit_glmm("IAH20_next ~ fluid_balance_t + iap_t + apache_t + age + sex + etiology + cr_t + map_t + (1|subject_id)",
                               data, "fluid_balance_t->IAH20 Model3")
  results_list[[length(results_list) + 1]] <- extract_binary_results(
    m3_balance_iah20, "fluid_balance_t", "IAH20_next", "Model 3", sd_fluid_balance_t, "Fluid balance", "IAH20"
  )

  # Combine all results
  results_df <- do.call(rbind, results_list)

  # Save convergence notes
  writeLines(convergence_notes, con = file.path(out_dir, "model_convergence_notes.txt"))

  # Save detailed results
  write.csv(results_df, file.path(out_dir, "Main_Association_Results_Per1000mL.csv"), row.names = FALSE)

  # ========== Create Forest Plots ==========

  # Prepare data for forest plots (Per SD only)
  forest_data <- results_df %>%
    dplyr::select(Exposure, Outcome_Label, Model, Effect_per_SD, CISD_Lower, CISD_Upper, P_value) %>%
    dplyr::mutate(
      Model_num = as.numeric(gsub("Model ", "", Model)),
      Outcome_order = case_when(
        Outcome_Label == "IAP (continuous)" ~ 1,
        Outcome_Label == "IAH15" ~ 2,
        Outcome_Label == "IAH20" ~ 3
      )
    )

  # Forest plot for Fluid intake
  p1 <- ggplot(forest_data %>% filter(Exposure == "Fluid intake"),
               aes(x = Effect_per_SD, y = interaction(Outcome_Label, Model_num),
                   color = Model)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CISD_Lower, xmax = CISD_Upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Model 1" = "#E74C3C", "Model 2" = "#3498DB", "Model 3" = "#2ECC71")) +
    labs(title = "Fluid Intake → IAP Outcomes (Per SD)",
         x = "Effect Size (β for continuous, OR for binary)",
         y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  # Forest plot for Fluid balance
  p2 <- ggplot(forest_data %>% filter(Exposure == "Fluid balance"),
               aes(x = Effect_per_SD, y = interaction(Outcome_Label, Model_num),
                   color = Model)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CISD_Lower, xmax = CISD_Upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Model 1" = "#E74C3C", "Model 2" = "#3498DB", "Model 3" = "#2ECC71")) +
    labs(title = "Fluid Balance → IAP Outcomes (Per SD)",
         x = "Effect Size (β for continuous, OR for binary)",
         y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  # Save individual plots
  ggsave(file.path(out_dir, "Forest_Plot_Fluid_Intake.png"), p1, width = 10, height = 8, dpi = 300)
  ggsave(file.path(out_dir, "Forest_Plot_Fluid_Balance.png"), p2, width = 10, height = 8, dpi = 300)

  # Combined plot
  p_combined <- grid.arrange(p1, p2, ncol = 1)
  ggsave(file.path(out_dir, "Forest_Plot_Combined.png"), p_combined, width = 10, height = 14, dpi = 300)

  cat("Forest plots generated successfully.\n")

  # ========== Create Formatted Table ==========

  # Create table with per 1000 mL and per SD side by side
  table_data <- results_df %>%
    dplyr::mutate(
      Effect_1000mL_CI = sprintf("%.2f (%.2f, %.2f)", Effect_per_1000mL, CI1000_Lower, CI1000_Upper),
      Effect_SD_CI = sprintf("%.2f (%.2f, %.2f)", Effect_per_SD, CISD_Lower, CISD_Upper),
      P_formatted = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
    ) %>%
    dplyr::select(Exposure, Outcome_Label, Model, Effect_1000mL_CI, Effect_SD_CI, P_formatted, N)

  # Pivot to wide format for easier reading
  table_wide <- table_data %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = c(Effect_1000mL_CI, Effect_SD_CI, P_formatted, N),
      names_sep = "_"
    )

  # Create flextable
  ft <- flextable(table_data) %>%
    set_header_labels(
      Exposure = "Exposure",
      Outcome_Label = "Outcome",
      Model = "Model",
      Effect_1000mL_CI = "Per 1000 mL (95% CI)",
      Effect_SD_CI = "Per SD (95% CI)",
      P_formatted = "P value",
      N = "N"
    ) %>%
    bold(part = "header") %>%
    align(align = "left", part = "all", j = 1:3) %>%
    align(align = "center", part = "all", j = 4:7) %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "header") %>%
    hline_bottom(border = fp_border(width = 2), part = "header") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    autofit()

  # Create Word document
  doc <- read_docx() %>%
    body_add_par("Table 2. Association Between Fluid Exposures and IAP Outcomes (Per 1000 mL and Per SD)",
                 style = "heading 1") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par("Model 1: Unadjusted model with exposure and random intercept for subject. Model 2: Adjusted for IAP at time t, APACHE II score, age, sex, and etiology. Model 3: Model 2 plus creatinine and mean arterial pressure. Effect sizes are β coefficients for continuous outcome (IAP) and odds ratios (OR) for binary outcomes (IAH15, IAH20). Per 1000 mL: effect per 1000 mL increase in exposure. Per SD: effect per 1 standard deviation increase in exposure. Mixed-effects models with random intercept were used when converged; otherwise, standard linear or logistic regression was used.",
                 style = "Normal")

  print(doc, target = file.path(out_dir, "Table2_Main_Association_Per1000mL.docx"))
  write.csv(table_data, file.path(out_dir, "Table2_Main_Association_Per1000mL.csv"), row.names = FALSE)

  cat("Table 2 generated successfully.\n")

  # ========== Write run log ==========
  end_time <- Sys.time()
  run_log <- c(
    paste0("start_time=", format(start_time, "%Y-%m-%d %H:%M:%S")),
    paste0("end_time=", format(end_time, "%Y-%m-%d %H:%M:%S")),
    "language=R",
    "script=scripts/main_association_analysis_per1000ml.R",
    paste0("packages=dplyr ", as.character(packageVersion("dplyr")),
           "; readxl ", as.character(packageVersion("readxl")),
           "; lme4 ", as.character(packageVersion("lme4")),
           "; lmerTest ", as.character(packageVersion("lmerTest")),
           "; officer ", as.character(packageVersion("officer")),
           "; flextable ", as.character(packageVersion("flextable")),
           "; ggplot2 ", as.character(packageVersion("ggplot2")),
           "; gridExtra ", as.character(packageVersion("gridExtra"))),
    paste0("output_dir=", out_dir),
    paste0("warnings=", ifelse(length(warnings_vec) == 0, "none", paste(unique(warnings_vec), collapse = " || "))),
    paste0("errors=", ifelse(length(errors_vec) == 0, "none", paste(unique(errors_vec), collapse = " || ")))
  )

  writeLines(run_log, con = file.path(out_dir, "run_log.txt"))
  cat(paste(run_log, collapse = "\n"), file = "logs/run_log.txt", append = TRUE, sep = "\n\n")

}, warning = function(w) {
  warnings_vec <<- c(warnings_vec, conditionMessage(w))
  invokeRestart("muffleWarning")
})

cat("\n=== Main association analysis completed ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files generated:\n")
cat("  - Main_Association_Results_Per1000mL.csv\n")
cat("  - Table2_Main_Association_Per1000mL.docx\n")
cat("  - Table2_Main_Association_Per1000mL.csv\n")
cat("  - Forest_Plot_Fluid_Intake.png\n")
cat("  - Forest_Plot_Fluid_Balance.png\n")
cat("  - Forest_Plot_Combined.png\n")
cat("  - model_convergence_notes.txt\n")
cat("  - run_log.txt\n")
