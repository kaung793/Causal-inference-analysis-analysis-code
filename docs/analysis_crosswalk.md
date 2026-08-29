# Analysis crosswalk

This crosswalk links each portable script to its analysis role and aggregate
checkpoint. Supplementary-table numbers refer to the current revision package
and should be rechecked if that package is reordered.

| Script | Analysis | Manuscript or supplement output |
|---|---|---|
| `00_validate_inputs.R` | Cohort dimensions, keys, lag sequence, and missingness checks | Reproducibility QA |
| `01_descriptive_missingness.R` | Window and patient summaries | Table 1 support and missing-data summary |
| `02_primary_mixed_models.R` | Model 1-3 LMM/GLMM associations for both fluid exposures | Main association table and standardized-effect figure inputs |
| `03_external_validation.R` | Limited external validation of fluid-intake associations | External-validation supplement |
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
| `16_figures.R` | Figures generated from reviewed aggregate outputs | Figure source-data workflow |

`run_all.R` executes the scripts in this order. Long MI and bootstrap analyses
are controlled by the configuration flags so that an accidental default run
does not launch hours of computation.

## Second-round manuscript-release entries

| Script | Role | Release output |
|---|---|---|
| `release/01_manuscript_aipw.R` | Reproduces the Table 3 AIPW estimates from the dedicated weight/longitudinal freezes with 1,000 patient-cluster resamples | `results/manuscript_release_v2/aipw.csv` |
| `release/02_manuscript_time_window.R` | Reconstructs the archived Table 4 display intervals and emits full-covariance companion intervals and diagnostics | `results/manuscript_release_v2/time_window.csv` |
| `release/03_check_manuscript_release.R` | Exact AIPW/RCS checks and package-sensitive time-window concordance checks | Generated verification summary |
| `run_manuscript_release.R` | Ordered release runner for the core and release-specific analyses | Run-status log plus generated outputs |

The complete document-to-code mapping and input roles are in
`results/manuscript_release_v2/release_manifest.csv`.
