# Second-round manuscript result contract

This directory is the non-identifying aggregate result contract for the
second-round revision package frozen in August 2026. It identifies the result
set actually used by the manuscript, tables, figures, supplement, and response
letter. Participant-level data are not included.

The release does **not** assume that one longitudinal file generated every
analysis. Several analyses were developed at different revision stages and had
documented, analysis-specific data freezes. `input_checksums.csv` records the
SHA-256 hash of each authorized private input without disclosing its contents.
`release_manifest.csv` maps each manuscript component to its code entry, input
role, and verification class.
`validation_summary.csv` records the completed local release check; rerunning
the verifier writes the full metric-level comparison to the configured output
directory.

## Verification classes

- **Exact deterministic rerun:** the released script, named input freeze, seed,
  and software environment reproduce the stored aggregate values.
- **Current-runtime exact:** the portable pipeline reproduces the aggregate
  result on the named current freeze.
- **Archived exact document output plus runtime concordance:** the exact
  aggregate file that generated the document is retained, while the current
  `lme4` runtime produces very small binary-model differences. The verifier
  checks a prespecified tolerance and reports the current diagnostics.
- **Manuscript rounded display:** the public file records the exact numbers as
  displayed in the revision package; higher-precision current estimates remain
  in `results/expected/` and generated outputs.

The time-window document used the historical display interval calculation that
added the main-effect and interaction variances without their covariance.
`release/02_manuscript_time_window.R` reproduces that calculation and also
emits the statistically complete covariance-matrix interval as a companion.
This distinction is explicit; the document interval is not relabelled as the
full-covariance interval.

Run `release/01_manuscript_aipw.R`,
`release/02_manuscript_time_window.R`, and then
`release/03_check_manuscript_release.R` with an authorized local configuration.
The root `run_manuscript_release.R` provides the ordered full entry point.
