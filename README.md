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
  matrix; and
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
results/expected/         Aggregate non-identifying checkpoints
tests/static_checks.R     Dependency-free repository checks
run_all.R                 Ordered runner
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

Long MI and bootstrap analyses are disabled in the example configuration. The
`--force` option runs one selected disabled script after input validation.

## Reproducibility boundary

Every completed analysis records input MD5 hashes, dimensions, seeds, and R
session information. Scripts fail on missing required columns and do not
silently replace a mixed model with ordinary regression. Singular fits and
convergence messages are retained in diagnostic output.

`results/expected/` contains aggregate checkpoints only; it is not a substitute
for source data. Some files deliberately contain both a historical document
checkpoint and a validated recomputation because several differences arose from
input versions, covariance calculations, or the earlier mediation workflow. Do
not mix values across those labels. Review
[`docs/reconciliation_notes.md`](docs/reconciliation_notes.md) before changing
manuscript numbers, figures, or legends.

## Privacy

Do not commit participant-level data, local configurations, run metadata,
fitted models, or generated results. The `.gitignore` excludes those locations
by default. Data access is governed by the manuscript's Data Availability
statement and applicable institutional approvals.

## License

No open-source license has been assigned. The code is publicly viewable for
research transparency; reuse permissions remain with the authors.
