# Fluid balance and intra-abdominal pressure: analysis code

This repository contains the R code used for the longitudinal analyses of fluid
balance and subsequent intra-abdominal pressure in severe acute pancreatitis.
Participant-level data are not included.

## Files

- `analysis/`: numbered scripts for the main, weighted, mediation, external-cohort,
  and sensitivity analyses
- `R/`: functions shared by the analysis scripts
- `config/example_config.R`: data paths and analysis settings
- `data/README.md`: required input files and variable names
- `run_all.R`: runs the enabled analyses in order

## Software

The analyses were run with R 4.5.3. Required packages are `lme4`, `lmerTest`,
`sandwich`, `mice`, `mediation`, `geepack`, `ordinal`, `zoo`, and
`ggplot2`. The `readxl` package is also needed when Excel input files are used.

## Running the analyses

1. Copy `config/example_config.R` to `config/local_config.R`.
2. Update the paths in the local configuration file.
3. Run the enabled analyses:

```sh
Rscript run_all.R --config=config/local_config.R
```

An individual script can also be run directly, for example:

```sh
Rscript analysis/02_primary_mixed_models.R --config=config/local_config.R
Rscript analysis/07_time_window_heterogeneity.R --config=config/local_config.R
```

Computationally intensive bootstrap and multiple-imputation analyses are disabled
in the example configuration and can be enabled when required.

## Data

Only de-identified, institutionally authorized data should be used. Local data,
configuration files, and generated results are excluded by `.gitignore`. The
external cohort is analyzed separately from the primary cohort and is not pooled
with it for the main estimates.

No open-source license has been assigned; reuse permission remains with the authors.
