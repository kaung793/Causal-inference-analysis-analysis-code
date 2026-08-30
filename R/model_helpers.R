scale_covariates <- function(data, vars) {
  d <- data
  for (v in intersect(vars, names(d))) {
    if (!is.numeric(d[[v]])) next
    s <- stats::sd(d[[v]], na.rm = TRUE)
    if (!is.finite(s) || s == 0) stop("Cannot scale constant covariate: ", v)
    d[[paste0(v, "_z")]] <- (d[[v]] - mean(d[[v]], na.rm = TRUE)) / s
  }
  d
}

mixed_convergence <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages
  opt <- fit@optinfo$conv$opt %||% 0L
  hard <- messages
  if (!is.null(hard)) hard <- hard[!grepl("boundary \\(singular\\) fit", hard, ignore.case = TRUE)]
  data.frame(
    optimizer_ok = all(opt == 0L) && (is.null(hard) || !length(hard)),
    optimizer_code = paste(opt, collapse = ";"),
    convergence_message = if (is.null(messages)) "" else paste(messages, collapse = " | "),
    singular = lme4::isSingular(fit, tol = 1e-5),
    n = stats::nobs(fit),
    stringsAsFactors = FALSE
  )
}

fit_lmm_strict <- function(formula, data, reml = FALSE) {
  require_packages(c("lme4", "lmerTest"))
  fit <- lmerTest::lmer(
    formula, data = data, REML = reml,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 300000))
  )
  status <- mixed_convergence(fit)
  if (!status$optimizer_ok) stop("LMM did not pass strict convergence checks: ", status$convergence_message)
  attr(fit, "analysis_diagnostics") <- status
  fit
}

fit_glmm_strict <- function(formula, data) {
  require_packages("lme4")
  attempts <- list()
  for (optimizer in c("bobyqa", "nloptwrap", "Nelder_Mead")) {
    fit <- tryCatch(
      lme4::glmer(
        formula, data = data, family = stats::binomial(), nAGQ = 1,
        control = lme4::glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 300000))
      ),
      error = identity
    )
    if (inherits(fit, "error")) {
      attempts[[optimizer]] <- conditionMessage(fit)
      next
    }
    status <- mixed_convergence(fit)
    status$optimizer <- optimizer
    if (status$optimizer_ok) {
      attr(fit, "analysis_diagnostics") <- status
      return(fit)
    }
    attempts[[optimizer]] <- status$convergence_message
  }
  stop("GLMM did not pass strict convergence checks: ", paste(unlist(attempts), collapse = " || "))
}

fit_glmm_with_warnings <- function(formula, data) {
  require_packages("lme4")
  fit <- lme4::glmer(
    formula, data = data, family = stats::binomial(), nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 300000))
  )
  status <- mixed_convergence(fit); status$optimizer <- "bobyqa"
  attr(fit, "analysis_diagnostics") <- status
  fit
}
fixed_coef <- function(model) {
  if (inherits(model, "merMod")) return(lme4::fixef(model))
  stats::coef(model)
}

extract_effect <- function(model, term, multiplier = 1, exponentiate = FALSE, level = 0.95) {
  b <- fixed_coef(model)
  if (!term %in% names(b)) stop("Term not found in fitted model: ", term)
  V <- as.matrix(stats::vcov(model))
  estimate <- unname(b[[term]]) * multiplier
  se <- sqrt(V[term, term]) * abs(multiplier)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  p <- 2 * stats::pnorm(-abs(estimate / se))
  lower <- estimate - crit * se
  upper <- estimate + crit * se
  if (exponentiate) c(estimate = exp(estimate), ci_lower = exp(lower), ci_upper = exp(upper), p_value = p,
                      linear_estimate = estimate, standard_error = se)
  else c(estimate = estimate, ci_lower = lower, ci_upper = upper, p_value = p,
         linear_estimate = estimate, standard_error = se)
}

contrast_vector <- function(model, terms) {
  b <- fixed_coef(model)
  L <- setNames(numeric(length(b)), names(b))
  unknown <- setdiff(names(terms), names(b))
  if (length(unknown)) stop("Contrast term(s) not found: ", paste(unknown, collapse = ", "))
  L[names(terms)] <- unname(terms)
  L
}

linear_contrast <- function(model, terms, exponentiate = FALSE, level = 0.95) {
  b <- fixed_coef(model)
  V <- as.matrix(stats::vcov(model))
  L <- contrast_vector(model, terms)
  estimate <- as.numeric(crossprod(L, b))
  variance <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(variance) || variance < 0) stop("Invalid contrast variance")
  se <- sqrt(variance)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  p <- 2 * stats::pnorm(-abs(estimate / se))
  lower <- estimate - crit * se
  upper <- estimate + crit * se
  if (exponentiate) c(estimate = exp(estimate), ci_lower = exp(lower), ci_upper = exp(upper),
                      p_value = p, linear_estimate = estimate, standard_error = se, variance = variance)
  else c(estimate = estimate, ci_lower = lower, ci_upper = upper,
         p_value = p, linear_estimate = estimate, standard_error = se, variance = variance)
}

joint_wald_test <- function(model, terms) {
  b <- fixed_coef(model)
  V <- as.matrix(stats::vcov(model))
  missing <- setdiff(terms, names(b))
  if (length(missing)) stop("Joint-test term(s) missing: ", paste(missing, collapse = ", "))
  q <- b[terms]
  U <- V[terms, terms, drop = FALSE]
  stat <- as.numeric(t(q) %*% qr.solve(U, q))
  data.frame(statistic = stat, df = length(terms), p_value = stats::pchisq(stat, length(terms), lower.tail = FALSE))
}

model_diagnostics <- function(model, label) {
  if (inherits(model, "merMod")) {
    d <- attr(model, "analysis_diagnostics") %||% mixed_convergence(model)
    if (!"optimizer" %in% names(d)) d$optimizer <- "bobyqa"
    vc <- as.data.frame(lme4::VarCorr(model))
    ri <- vc$vcov[vc$var1 == "(Intercept)" & is.na(vc$var2)]
    return(data.frame(model = label, d, random_intercept_variance = if (length(ri)) ri[[1L]] else NA_real_))
  }
  data.frame(model = label, optimizer_ok = TRUE, optimizer_code = "not_applicable", optimizer = "not_applicable",
             convergence_message = "", singular = NA, n = stats::nobs(model),
             random_intercept_variance = NA_real_)
}

model3_rhs <- function(exposure, random_intercept = TRUE, scaled = FALSE) {
  covars <- if (scaled) c("iap_current_z", "apache_ii_z", "age_z", "sex", "etiology", "creatinine_z", "map_z")
  else c("iap_current", "apache_ii", "age", "sex", "etiology", "creatinine", "map")
  rhs <- paste(c(exposure, covars), collapse = " + ")
  if (random_intercept) paste(rhs, "+ (1 | subject_id)") else rhs
}
