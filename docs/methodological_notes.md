# Methodological implementation notes

## Longitudinal models

The primary models relate exposure in one scheduled window to IAP at the next
scheduled assessment and include a patient random intercept. Continuous
adjustment variables are standardized only to improve numerical conditioning;
the exposure remains expressed per 1,000 mL. GLMMs are tried with declared
optimizers and do not silently fall back to ordinary logistic regression.

Singular random-intercept fits are retained and reported because a zero
estimated random-intercept variance is a model diagnostic, not an optimizer
failure. The expanded-comorbidity sensitivity retains its documented
identifiability warning and labels that exception explicitly.

## Window contrasts

Window-specific effects are linear combinations of the exposure main effect and
interaction coefficients. Standard errors use the complete model covariance
matrix, including the covariance between those coefficients.

## Multiple imputation

Longitudinal MI is performed in patient-wide form to preserve within-patient
structure. The Model 3 sensitivity uses its adjudicated pre-imputation data
freeze, a nonduplicated Day 1/2/3/5/7 IAP sequence, four etiology categories,
one pre-specified deterministic MAP reconstruction, and an analysis-specific
seed. The data are converted back to long form after imputation. Every model is
fitted in every completed dataset, and estimates and variances are pooled with
Rubin's rules; MICE events and mixed-model diagnostics are retained.

For mediation, early cumulative estimated fluid balance is a trapezoidal AUC
over Days 1-3. ACME, ADE, total effect, and their covariance are estimated
within each completed dataset. The mediated proportion uses the pooled
ACME-total covariance and a logit-scale delta interval to keep bounds between
zero and one.

## Dependence-aware uncertainty

Repeated-window GEE analyses use patient-clustered sandwich variance.
Standardized absolute risks and partial proportional-odds sensitivity analyses
use patient-cluster bootstrap resampling. All nuisance models and weights are
re-estimated inside each applicable bootstrap replicate.
