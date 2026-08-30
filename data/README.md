# Input data

Participant-level data are not distributed with the code. Place de-identified,
institutionally authorized files under `data/` and set their paths in
`config/local_config.R`. CSV, TSV, XLS/XLSX, and RDS files are supported.

## Primary longitudinal file

The primary file has one row per patient and observation window. The full panel
contains 801 patients and 3,204 windows beginning on Days 1, 2, 3, and 5.

| Variable | Accepted aliases | Description |
|---|---|---|
| `subject_id` | `analysis_id`, `id` | De-identified patient identifier |
| `day` | `start_day` | Start day of the observation window |
| `fluid_intake_ml` | `fluid_t`, `fluid_input_ml` | Recorded 24-hour fluid intake |
| `estimated_fluid_balance_ml` | `fluid_balance_t`, `fluid_balance_ml` | Intake minus urine output |
| `iap_current` | `iap_t`, `iap_start` | Daily maximum IAP at the start of the window |
| `iap_next` | `iap_end` | Daily maximum IAP at the subsequent retained assessment |
| `apache_ii` | `apache_t`, `apache` | APACHE II score |
| `creatinine` | `cr_t`, `cr` | Serum creatinine |
| `map` | `map_t`, `map_repaired` | Mean arterial pressure |
| `age` |  | Age in years |
| `sex` | `sex_code` | Sex category |
| `etiology` | `etiology_code` | Pancreatitis etiology |

Additional analyses may use `urine_output_ml`, `weight_kg`, `heart_rate`,
`albumin`, `ph`, `spo2`, `shock`, `pcd`, and `crrt`. Binary IAP
outcomes are derived from `iap_next`.

The optional `rcs_long` and `time_window_long` files use the same structure.
If their configuration values are `NULL`, the primary longitudinal file is used.

## Other analysis files

- `longitudinal_revision`: longitudinal data used for the additional sensitivity
  analyses; it follows the primary-file structure.
- `revision_covariates`: optional subject-day file containing `subject_id`,
  `day`, and corrected covariate values.
- `longitudinal_external`: external-cohort data with `subject_id`, `day`,
  `fluid_intake_ml`, `iap_current`, `iap_next`, `map`, `age`, `sex`,
  and `etiology`. It contains 171 patients and 684 windows.
- `longitudinal_primary_preimputation`: primary longitudinal data before
  multiple imputation; it also retains systolic and diastolic blood pressure.
- `equal_lag_long`: one row per patient and two-day window, with
  `subject_id`, `window_start`, `window_end`, `fluid_balance_ml`,
  `iap_start`, `iap_end`, `apache_start`, `creatinine_start`,
  `map_start`, `age`, `sex`, and `etiology`.
- `mediation_patient`: one row per patient with `subject_id`, `fluid_auc`,
  `iap_day5`, `iap_day1`, `age`, `sex`, `etiology`, `apache`,
  `creatinine`, `map`, and `death_hospital`.
- `patient_flags`: optional patient-level indicators for subgroup and exclusion
  analyses.

For the mediation analysis, `fluid_auc` is the trapezoidal area under the
Day 1-Day 3 fluid-balance curve:

```text
fluid_auc = (Day1 + Day2) / 2 + (Day2 + Day3) / 2
```

## AIPW files

`aipw_weights` contains one row per eligible patient-window and the variables
`subject_id`, `day`, `high_fluid_balance`, `fluid_balance_t`, `age`,
`sex`, `etiology`, `iap_t`, `apache_t`, and `IAH15_next`.

`aipw_long` contains the corresponding longitudinal records with `subject_id`,
`day`, `fluid_balance_t`, `age`, `sex`, `etiology`, `iap_t`,
`apache_t`, `cr_t`, `map_t`, and `iap_next`.

Local data and configuration files should not be committed to the repository.
