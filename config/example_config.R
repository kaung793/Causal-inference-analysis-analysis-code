# Copy this file to config/local_config.R and update the local data paths.

config <- list(
  primary_long = "data/longitudinal_primary.csv",
  revision_primary_long = "data/longitudinal_revision.csv",
  model3_mi_long = "data/longitudinal_primary_preimputation.csv",
  external_long = "data/longitudinal_external.csv",
  mediation_patient = "data/mediation_patient.csv",
  equal_lag_long = "data/equal_lag.csv",
  patient_flags = "data/patient_flags.csv",
  revision_covariates = "data/revision_covariates.csv",

  # These may be left NULL when the primary longitudinal file is used.
  rcs_long = NULL,
  time_window_long = NULL,

  # Analysis-ready files used by the AIPW analysis.
  aipw_weights = "data/aipw_weights.csv",
  aipw_long = "data/longitudinal_aipw.csv",

  output_dir = "results/generated",
  figure_dir = "figures/generated",

  expected_primary_patients = 801L,
  expected_primary_windows = 3204L,
  expected_external_patients = 171L,
  expected_external_windows = 684L,

  seed = 20260723L,
  bootstrap_reps = 1000L,
  mediation_sims = 1000L,
  mi_m = 20L,
  mi_maxit = 20L,
  equal_lag_mi_maxit = 100L,
  mi_donors = 5L,
  model3_mi_seed = 20260716L,
  model3_mi_map_reconstruction_day = 3L,
  model3_mi_map_reconstruction_n = 1L,
  aipw_seed = 123L,
  aipw_bootstrap_reps = 1000L,

  run_external_validation = FALSE,
  run_msm = TRUE,
  run_aipw = FALSE,
  run_equal_lag_mi = FALSE,
  run_mediation = FALSE,
  run_model3_mi = FALSE,
  run_revision_sensitivities = FALSE,
  run_figures = TRUE,

  strict_cohort_checks = TRUE
)
