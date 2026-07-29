# Data contract

Participant-level data are not distributed in this repository. To reproduce the
analyses, place locally authorized files under `data/` and copy
`config/example_config.R` to `config/local_config.R`.

Files may be CSV, TSV, XLS/XLSX, or RDS. Column aliases used by legacy
data exports are accepted, but the canonical names below are preferred.

## `longitudinal_primary`

One row per patient and observation window. The validated primary panel contains
801 patients and 3,204 windows before complete-case restrictions.

Required columns:

| Canonical name | Historical aliases | Definition |
|---|---|---|
| `subject_id` | `analysis_id`, `id` | De-identified patient identifier |
| `day` | `start_day` | Start day of the observation window: 1, 2, 3, or 5 |
| `fluid_intake_ml` | `fluid_t`, `fluid_input_ml` | Fluid intake during the exposure window, mL |
| `fluid_balance_ml` | `estimated_fluid_balance_ml`, `fluid_balance_t` | Fluid balance during the exposure window, mL |
| `iap_current` | `iap_t`, `iap_start` | IAP at the start of the window, mmHg |
| `iap_next` | `iap_end` | IAP at the end of the window, mmHg |
| `apache_ii` | `apache_t`, `apache` | APACHE II score at the start of the window |
| `creatinine` | `cr_t`, `cr` | Creatinine at the start of the window |
| `map` | `map_t`, `map_repaired` | Mean arterial pressure at the start of the window |
| `age` | - | Age in years |
| `sex` | `sex_code` | Sex category |
| `etiology` | `etiology_code` | Etiology category |

Optional columns used by revision analyses include `urine_output_ml`,
`weight_kg`, `heart_rate`, `albumin`, `ph`, `spo2`, `shock`, `pcd`, and `crrt`.
Threshold outcomes are derived from `iap_next` and need not be supplied.

`fluid_balance_ml` is the study-specific measure calculated from recorded fluid
intake and urine output. It does not represent complete physiological net fluid
balance because other losses and net ultrafiltration are not fully captured.

## `longitudinal_primary_preimputation`

This is the adjudicated pre-imputation derivation-data freeze used only for the
patient-level wide-format Model 3 MICE sensitivity analysis. It follows the
primary schema and must additionally retain `sbp`/`sbp_t` and `dbp`/`dbp_t`.
The script reconstructs the one pre-identified Day 3 MAP record from observed
SBP and DBP without exposing a subject identifier.

Keep this input distinct from later analysis-ready freezes: the imputation
model uses a nonduplicated Day 1/2/3/5/7 IAP sequence, four etiology categories,
20 imputations, 20 iterations, five predictive-mean-matching donors, and the
analysis-specific seed recorded in `config/example_config.R`.

## `longitudinal_revision` and `revision_covariates`

`longitudinal_revision` is the frozen longitudinal source used by analyses that
were added during peer review. It follows the primary schema. Keeping its path
separate prevents results from an earlier reviewer-analysis freeze from being
silently combined with a later primary-data correction.

`revision_covariates` is an optional subject-day file used for adjudicated
corrections such as albumin, pH, or MAP. It must contain `subject_id` and `day`.
Supported fields include `fluid_balance_ml`, `iap_current`,
`iap_next`, `apache_ii`, `creatinine`, `map`, `age`, `weight_kg`,
`heart_rate`, `albumin`, `ph`, and `spo2`. Supplemental nonmissing values take
precedence, and duplicate subject-day keys cause an error.

## `longitudinal_external`

One row per patient and observation window. The validated external panel contains
171 patients and 684 windows. Required columns are `subject_id`, `day`,
`fluid_intake_ml`, `iap_current`, `iap_next`, `map`, `age`, `sex`, and
`etiology`. The binary IAP outcomes are derived from `iap_next`.

## `mediation_patient`

One row per patient. Required columns are `subject_id`, `fluid_auc`, `iap_day5`,
`iap_day1`, `age`, `sex`, `etiology`, `apache`, `creatinine`, `map`, and
`death_hospital`.

`fluid_auc` is the trapezoidal area under the fluid-balance curve from
Day 1 through Day 3, in mL-day:

```text
fluid_auc = (Day1 + Day2) / 2 + (Day2 + Day3) / 2
```

The exposure, mediator, and outcome are observed rather than imputed. Missing
covariates are multiply imputed inside the mediation script.

## `equal_lag_long`

One row per patient and equal two-day window. This file uses the same canonical
fields as the primary longitudinal data. Required columns are `subject_id`, `day`
(start day: 1, 3, or 5), `fluid_intake_ml`, `fluid_balance_ml`, `iap_current`,
`iap_next`, `apache_ii`, `creatinine`, `map`, `age`, `sex`, and `etiology`.
The corresponding intervals are Day 1 to 3, Day 3 to 5, and Day 5 to 7.

## `patient_flags`

Optional one-row-per-patient file keyed by `subject_id`. It may contain `shock`,
`pcd`, and `crrt` indicators when these variables are not already present in the
primary longitudinal file. Expanded-comorbidity analyses additionally accept
`smoking`, `drinking`, `hypertension`, `diabetes`,
`hyperlipidemia_history`, and `copd`.

## Privacy and validation

Use only de-identified, institutionally authorized data. Do not commit local
data, configuration files, run metadata, fitted models, or generated results.
Run `analysis/00_validate_inputs.R` before any inferential analysis.
