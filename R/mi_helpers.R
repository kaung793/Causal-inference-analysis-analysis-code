rubin_pool_scalar <- function(estimates, variances, conf.level = 0.95) {
  keep <- is.finite(estimates) & is.finite(variances) & variances >= 0
  q <- estimates[keep]
  u <- variances[keep]
  m <- length(q)
  if (m < 2L) stop("Rubin pooling requires at least two valid estimates")
  qbar <- mean(q)
  ubar <- mean(u)
  between <- stats::var(q)
  total <- ubar + (1 + 1 / m) * between
  se <- sqrt(total)
  r <- if (ubar > 0) (1 + 1 / m) * between / ubar else Inf
  df <- if (is.finite(r) && r > 1e-12) (m - 1) * (1 + 1 / r)^2 else Inf
  alpha <- 1 - conf.level
  crit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df) else stats::qnorm(1 - alpha / 2)
  stat <- qbar / se
  p <- if (is.finite(df)) 2 * stats::pt(-abs(stat), df) else 2 * stats::pnorm(-abs(stat))
  data.frame(
    m = m, estimate = qbar, within_variance = ubar, between_variance = between,
    total_variance = total, standard_error = se, df = df,
    ci_lower = qbar - crit * se, ci_upper = qbar + crit * se,
    p_value = p, fraction_missing_information = if (total > 0) (1 + 1 / m) * between / total else 0
  )
}

rubin_pool_grouped <- function(data, group_vars = c("exposure", "outcome"),
                               estimate = "estimate", variance = "variance") {
  key <- interaction(data[group_vars], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(data, key), function(z) {
    pooled <- rubin_pool_scalar(z[[estimate]], z[[variance]])
    cbind(z[1L, group_vars, drop = FALSE], pooled)
  })
  out <- do.call(rbind, pieces)
  row.names(out) <- NULL
  out
}

pool_vector_rubin <- function(Q, U) {
  Q <- as.matrix(Q)
  if (length(dim(U)) != 3L || dim(U)[3L] != nrow(Q)) stop("U must be p x p x m")
  qbar <- colMeans(Q)
  ubar <- apply(U, c(1, 2), mean)
  between <- stats::cov(Q)
  total <- ubar + (1 + 1 / nrow(Q)) * between
  list(estimate = qbar, within = ubar, between = between, total = total)
}

bounded_proportion_ci <- function(estimate, se, conf.level = 0.95) {
  if (!is.finite(estimate) || !is.finite(se) || estimate <= 0 || estimate >= 1) return(c(NA_real_, NA_real_))
  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  se_logit <- se / (estimate * (1 - estimate))
  stats::plogis(stats::qlogis(estimate) + c(-1, 1) * crit * se_logit)
}

make_patient_wide_imputation <- function(wide, targets, m, maxit, seed, donors = 5L,
                                         id_column = "subject_id") {
  require_packages("mice")
  assert_columns(wide, c(id_column, targets), "wide imputation data")
  ini <- mice::mice(wide, m = 1, maxit = 0, printFlag = FALSE,
                    remove.constant = FALSE, remove.collinear = FALSE)
  method <- ini$method
  method[] <- ""
  method[targets] <- "pmm"
  pred <- ini$predictorMatrix
  pred[,] <- 0
  predictors <- setdiff(names(wide), id_column)
  for (target in targets) pred[target, setdiff(predictors, target)] <- 1
  pred[, id_column] <- 0
  pred[id_column, ] <- 0
  mice::mice(
    wide, m = m, maxit = maxit, method = method, predictorMatrix = pred,
    visitSequence = targets, donors = donors, seed = seed, printFlag = FALSE,
    remove.constant = FALSE, remove.collinear = FALSE
  )
}
