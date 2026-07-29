# Fluid balance and intra-abdominal pressure: analysis code

This repository provides the R analysis workflows supporting the accompanying
severe acute pancreatitis manuscript. Participant-level data are not included.

## Scope

The ordered scripts cover input validation, descriptive summaries, primary
mixed-effects models, external validation, dose-response analyses, weighting,
AIPW, time-window analyses, equal-lag sensitivity analyses, multiple imputation,
mediation, additional outcome definitions, subgroup analyses, and figure-source
workflows.

## Repository layout

```text
analysis/                 Ordered analysis scripts 00-16
R/                        Input, model, pooling, and QA helpers
config/example_config.R   Portable paths and run settings
data/README.md            Required files, schemas, and privacy boundary
docs/                     Analysis crosswalk, methods, and software environment
tests/static_checks.R     Dependency-free repository checks
run_all.R                 Ordered runner
```

Generated results, fitted models, run logs, and figures are written to ignored
local directories and are not part of the public source tree.

## Data terminology

`fluid_balance_ml` is the canonical public variable name. It represents the
study-specific fluid-balance measure defined from recorded fluid intake and
urine output. It should not be interpreted as complete physiological net fluid
balance because other losses and net ultrafiltration are not fully captured.
The input layer continues to accept the historical aliases
`estimated_fluid_balance_ml` and `fluid_balance_t`.

## Quick start

1. Install R 4.5.3 and the package versions listed in
   [`docs/software_environment.md`](docs/software_environment.md).
2. Prepare only institutionally authorized, de-identified data using
   [`data/README.md`](data/README.md).
3. Copy `config/example_config.R` to `config/local_config.R` and update the local
   paths. The local file is ignored by Git.
4. Validate the inputs:

```sh
Rscript analysis/00_validate_inputs.R --config=config/local_config.R
```

5. Run the analyses enabled in the configuration:

```sh
Rscript run_all.R --config=config/local_config.R
```

An individual analysis can also be run directly:

```sh
Rscript analysis/07_time_window_heterogeneity.R --config=config/local_config.R
Rscript analysis/09_mediation_mi_rubin.R --config=config/local_config.R --force
```

Long multiple-imputation and bootstrap analyses are disabled in the example
configuration. The `--force` option runs one selected disabled analysis after
input preparation.

## Reproducibility boundary

Completed analyses record input hashes, dimensions, seeds, and R session
information. Scripts fail on missing required columns and retain mixed-model
convergence and singularity diagnostics. Numerical reproduction additionally
requires the authorized frozen input files described in the data contract.

For a formal release, cite a fixed Git tag rather than the moving `main` branch.
Release archives and their checksums should be attached to the corresponding
GitHub Release rather than committed to the source tree.

## Privacy

Do not commit participant-level data, local configurations, run metadata,
fitted models, generated results, or figures. The `.gitignore` excludes the
standard local locations for these artifacts. Data access remains governed by
the manuscript's Data Availability statement and applicable institutional
approvals.

## License

No open-source license has been assigned. The code is publicly viewable for
research transparency; reuse permissions remain with the authors.
