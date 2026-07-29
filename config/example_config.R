# Copy this file to config/local_config.R for local use.
# Paths may be absolute or relative to the repository root.

config <- list(
  primary_long = "data/longitudinal_primary.csv",
  revision_primary_long = "data/longitudinal_revision.csv",
  model3_mi_long = "data/longitudinal_primary_preimputation.csv",
  external_long = "data/longitudinal_external.csv",
  mediation_patient = "data/mediation_patient.csv",
  equal_lag_long = "data/equal_lag.csv",
  patient_flags = "data/patient_flags.csv",
  revision_covariates = "data/revision_covariates.csv",
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
  parallel_workers = 1L,

  run_external_validation = FALSE,
  run_msm = TRUE,
  run_aipw = FALSE,
  run_time_window_mi = FALSE,
  run_equal_lag_mi = FALSE,
  run_mediation = FALSE,
  run_model3_mi = FALSE,
  run_revision_sensitivities = FALSE,
  run_figures = TRUE,

  # The manuscript AIPW intervals used row-level bootstrap resampling. The
  # patient-cluster bootstrap is also emitted as a recommended sensitivity.
  aipw_bootstrap_units = c("row", "patient"),

  # Input data are expected to have already passed source-record adjudication.
  # Strict mode enforces the manuscript cohort dimensions and key invariants.
  strict_cohort_checks = TRUE
)
