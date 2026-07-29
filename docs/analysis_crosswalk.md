# Analysis crosswalk

This crosswalk links each script to its analysis role and corresponding manuscript
or supplement output. Table numbers refer to the current revision package and
should be rechecked if that package is reordered.

| Script | Analysis | Manuscript or supplement output |
|---|---|---|
| `00_validate_inputs.R` | Cohort dimensions, keys, lag sequence, and missingness checks | Reproducibility QA |
| `01_descriptive_missingness.R` | Window and patient summaries | Descriptive and missing-data support |
| `02_primary_mixed_models.R` | Model 1-3 LMM/GLMM associations for both fluid exposures | Main association table inputs |
| `03_external_validation.R` | External-cohort fluid-intake estimates | External-validation supplement inputs |
| `04_rcs_dose_response.R` | Restricted cubic splines and descriptive curve features | Dose-response figure and RCS supplement |
| `05_msm_iptw.R` | Stabilized and truncated IPTW marginal structural models | Weighted-analysis table |
| `06_aipw.R` | Primary and extended augmented-IPW analyses | Weighted-analysis table |
| `07_time_window_heterogeneity.R` | Unequal-window interaction models with full-covariance contrasts | Figure 5 and Table 4 |
| `08_equal_lag_sensitivity.R` | Equal two-day-window complete-case and optional MI analyses | Equal-lag sensitivity supplement |
| `09_mediation_mi_rubin.R` | Patient-level AUC exposure, 20 imputations, per-imputation mediation, and Rubin pooling | Table 5 and mediation figure |
| `10_alternative_iap_outcomes.R` | IAP >=12 mmHg, WSACS ordinal grade, and partial proportional odds | Supplementary Table 14 |
| `11_additional_adjustments.R` | Corrected Model 4, pH, SpO2, GEE, and expanded-comorbidity models | Supplementary Table 15 |
| `12_exclusions_subgroups.R` | PCD/CRRT exclusions, weight normalization, and shock strata/interactions | Supplementary Tables 16-18 |
| `13_standardized_absolute_risks.R` | IPTW-standardized risks, RD, RR, and OR with patient bootstrap | Supplementary absolute-risk table |
| `14_model3_multiple_imputation.R` | Dedicated pre-imputation freeze, patient-wide MICE, and Rubin-pooled Model 3 estimates | Supplementary MI table |
| `15_repeated_threshold.R` | Repeated IAP >=12 and >=15 mmHg outcomes | Repeated-threshold sensitivity supplement |
| `16_figures.R` | Diagnostic/source plots generated from analysis outputs | Figure-development support |

`run_all.R` executes the scripts in this order. Long MI and bootstrap analyses
are controlled by the configuration flags so that an accidental default run
does not launch hours of computation.
