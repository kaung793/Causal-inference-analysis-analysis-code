# Fluid balance and intra-abdominal pressure: analysis code

This repository contains the R code for the longitudinal, weighted, mediation,
external-validation, and sensitivity analyses in the accompanying severe acute
pancreatitis manuscript. Participant-level data are not included.

## What this revision adds

The original repository contained seven historical scripts with local paths and
did not include many analyses added during peer review. This revision:

- removes machine-specific paths, `setwd()`, package auto-installation, and
  private filenames;
- adds a documented data contract, portable configuration, ordered runner,
  model diagnostics, run metadata, and aggregate checkpoints;
- adds external validation, time-window heterogeneity, equal-lag models,
  alternative IAP outcomes, additional adjustments, exclusions/subgroups,
  standardized absolute risks, multiple imputation, repeated-threshold analyses,
  and figure-source workflows;
- implements mediation in every completed dataset and pools estimates and
  covariance with Rubin's rules;
- calculates window-specific intervals from the full coefficient covariance
  matrix;
- adds a versioned second-round manuscript result contract, analysis-specific
  input hashes, a dedicated patient-cluster AIPW rerun, and an automated
  release verifier; and
- keeps the seven historical root filenames as compatibility entry points.

The compatibility files dispatch to the reviewed scripts under `analysis/`; the
original implementations remain available in Git history.

## Repository layout

```text
analysis/                 Ordered scripts 00-16
R/                        Input, model, pooling, and QA helpers
config/example_config.R   Portable paths and run settings
data/README.md            Required files, schemas, and privacy boundary
docs/                     Crosswalk, methods, environment, reconciliation notes
release/                  Manuscript-specific reruns and release verification
results/expected/         Aggregate non-identifying checkpoints
results/manuscript_release_v2/  Selected second-round result contract
tests/static_checks.R     Dependency-free repository checks
run_all.R                 Ordered runner
run_manuscript_release.R  Second-round manuscript release runner
```

## Quick start

1. Install R 4.5.3 and the package versions listed in
   [`docs/software_environment.md`](docs/software_environment.md).
2. Prepare only institutionally authorized, de-identified data using
   [`data/README.md`](data/README.md).
3. Copy `config/example_config.R` to `config/local_config.R` and update local
   paths. The local file is ignored by Git.
4. Validate the inputs:

```sh
Rscript analysis/00_validate_inputs.R --config=config/local_config.R
```

5. Run all analyses enabled in the configuration:

```sh
Rscript run_all.R --config=config/local_config.R
```

An individual analysis can be run directly, for example:

```sh
Rscript analysis/07_time_window_heterogeneity.R --config=config/local_config.R
Rscript analysis/09_mediation_mi_rubin.R --config=config/local_config.R
```

To verify the result set used in the second-round submission, configure the
analysis-specific private freezes described in `data/README.md`, set
`run_manuscript_release=TRUE`, and run:

```sh
Rscript run_manuscript_release.R --config=config/local_config.R
```

The three release-specific checks can also be run separately:

```sh
Rscript release/01_manuscript_aipw.R --config=config/local_config.R
Rscript release/02_manuscript_time_window.R --config=config/local_config.R
Rscript release/03_check_manuscript_release.R --config=config/local_config.R
```

Long MI and bootstrap analyses are disabled in the example configuration. The
`--force` option runs one selected disabled script after input validation.

## Reproducibility boundary

Every completed analysis records input MD5 hashes, dimensions, seeds, and R
session information. Scripts fail on missing required columns and do not
silently replace a mixed model with ordinary regression. Singular fits and
convergence messages are retained in diagnostic output.

`results/expected/` contains the broader audit history. The exact set selected
for the second-round documents is stated separately in
[`results/manuscript_release_v2/`](results/manuscript_release_v2/), including
the analysis entry point, private-input hash, and verification class for each
component. This makes analysis-specific freezes explicit instead of asking a
reader to infer a release by combining historical labels.

The time-window document retained its archived variance-sum display intervals;
the release script also emits full-covariance intervals and model diagnostics.
The manuscript AIPW estimates are reproduced from their dedicated freeze using
1,000 patient-cluster bootstrap resamples. Review
[`docs/reconciliation_notes.md`](docs/reconciliation_notes.md) for the complete
scientific and computational boundary.

## Privacy

Do not commit participant-level data, local configurations, run metadata,
fitted models, or generated results. The `.gitignore` excludes those locations
by default. Data access is governed by the manuscript's Data Availability
statement and applicable institutional approvals.

## License

No open-source license has been assigned. The code is publicly viewable for
research transparency; reuse permissions remain with the authors.
