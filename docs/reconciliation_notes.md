# Reconciliation notes

## Purpose

The original root scripts preserve an important historical trace, but they used
machine-specific paths and did not cover analyses added during revision. The
portable pipeline separates two questions:

1. What values were carried forward in the revision documents?
2. What values are produced by a fully specified, diagnostically acceptable
   rerun on a named input version?

Aggregate checkpoints for both questions are stored in `results/expected/`.
Every actual run also records input hashes and software versions. Values from
different input hashes must not be combined.

## Items requiring explicit synchronization

### Primary mixed models

The historical binary models fitted the unscaled adjustment variables and
emitted convergence or identifiability warnings. The portable implementation
uses an algebraically equivalent scaled parameterization and requires a
successful optimizer. Rounded Model 3 conclusions are unchanged, although some
unrounded ORs differ slightly. The diagnostics file, rather than only the
coefficient table, is part of the result.

### External validation

The current corrected 171-patient, 684-window input produces a continuous-IAP
beta of 0.125 and an IAP >=15 mmHg OR of 1.139 per 1,000 mL. An older supplement
checkpoint used 0.112 and 1.079. This is an input-version difference; the
release package must freeze one external file and use its hash consistently.

### Time-window heterogeneity

The day-specific point estimates were generated in the Model 3 complete-case
sample of 771 patients and 2,188 windows. They were not generated in a
20-imputation sample. Historical confidence intervals for later windows added
the main-effect and interaction variances while omitting their covariance.
`07_time_window_heterogeneity.R` uses the full covariance matrix.

The historical binary interaction fit also had convergence warnings. The
portable scaled fit gives ORs of 1.031, 1.155, 1.350, and 1.328 across the four
windows. Figure legends and tables should identify the complete-case sample and
use one internally consistent set of estimates and intervals.

### Model 3 multiple imputation

The validated Model 3 MI sensitivity uses a dedicated adjudicated
pre-imputation freeze rather than the later analysis-ready primary file. Its
predictor matrix contains one Day 1/2/3/5/7 IAP sequence, four etiology
categories, one pre-specified MAP reconstruction, and the recorded seed. With
20 imputations, the portable script reproduces beta 0.317 and ORs 1.236 and
1.266 for estimated fluid balance. All 120 mixed-model fits pass optimizer
checks, 60 boundary fits are retained, and no MICE event is logged.

### Mediation

The historical workflow created 20 completed datasets but estimated mediation
only in the first completed dataset. It therefore did not constitute pooled MI
inference. `09_mediation_mi_rubin.R` estimates ACME, ADE, and total effect in
every completed dataset and pools estimates and covariance with Rubin's rules.

The validated pooled mediated proportions are 14.3% for the continuous
1,000-mL-day exposure and 15.3% for the median-defined binary exposure. The
revision-document checkpoints of 13.5% and 12.1% came from the historical
implementation. The manuscript, abstract, table, figure legend, discussion,
conclusion, and response letter must use the same declared analysis.

### AIPW and RCS

AIPW point estimates are sensitive to the frozen longitudinal input version.
The portable script emits both row-bootstrap intervals, for compatibility, and
patient-cluster bootstrap intervals, as the recommended dependence-aware
sensitivity. RCS curve-derived values are descriptive features, not validated
physiological or treatment thresholds. The three fluid-balance inflection
features near 1,900 mL are reproduced.

## Release checklist

- Freeze each input file and retain its MD5 in the run metadata.
- Run `tests/static_checks.R`.
- Run `analysis/00_validate_inputs.R` before inferential scripts.
- Review every nonempty convergence message and every singular-fit flag.
- Compare generated aggregate outputs with `results/expected/`.
- Synchronize the manuscript, supplement, figures, abstract, and response
  letter before merging or tagging a release.
- Never commit participant-level data, local configuration, fitted models, or
  generated run metadata.
