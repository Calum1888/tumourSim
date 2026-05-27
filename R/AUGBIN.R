# =============================================================================
# AugBin (augmented binary) estimator and tests
# =============================================================================
# Implementation of the Lin and Wason (2020) parametric model for PFS:
#
#   S_PFS(t) = S_T(t) * S_D(t),
#
# where S_T(t) is computed from the multivariate normal model on tumour
# log-ratios (no progression up to t = all pairwise differences below the
# log-threshold) and S_D(t) is computed from a longitudinal logistic
# regression for new lesions with cumulative-product survival.
#
# Both bootstrap and delta-method variance estimators are provided.
# =============================================================================


#' Pairwise transform matrix for the tumour-growth survival
#'
#' Builds the linear transform \eqn{A_t} such that \eqn{A_t Y_{1:t}} stacks
#' the \eqn{t} log-ratios versus baseline followed by all
#' \eqn{\binom{t}{2}} pairwise differences \eqn{Y_a - Y_b} (\eqn{1 \le b < a \le t}).
#' The joint "no progression up to visit \eqn{t}" event is then
#' \eqn{\{\text{all rows of } A_t Y_{1:t} < \log c\}}, which is
#' multivariate normal with mean \eqn{A_t \mu_{1:t}} and covariance
#' \eqn{A_t \Sigma_{1:t,1:t} A_t^\top}.
#'
#' @param t Integer. Number of visits up to which to build the transform.
#'
#' @return A numeric matrix with \eqn{t + \binom{t}{2}} rows and \eqn{t}
#'   columns.
#'
build_pairwise_A <- function(t) {
  n_baseline <- t
  pairs <- if (t >= 2) {
    do.call(rbind, lapply(2:t, function(a) cbind(a, 1:(a - 1))))
  } else {
    matrix(integer(0), 0, 2)
  }
  n_pairs <- nrow(pairs)
  A <- matrix(0, nrow = n_baseline + n_pairs, ncol = t)
  for (i in 1:t) A[i, i] <- 1
  if (n_pairs > 0) {
    for (k in seq_len(n_pairs)) {
      a <- pairs[k, 1]
      b <- pairs[k, 2]
      A[n_baseline + k, a] <-  1
      A[n_baseline + k, b] <- -1
    }
  }
  A
}


#' Tumour-side survival from MVN parameters
#'
#' Computes \eqn{S_T(t)}, the probability of no tumour progression by visit
#' \eqn{t}, given the mean vector and covariance matrix of the log-ratio
#' process and a multiplicative progression threshold. Uses the
#' multivariate-normal CDF via \code{\link[mvtnorm]{pmvnorm}}.
#'
#' @param mu Numeric vector of length \code{T} giving the log-ratio mean.
#' @param Sigma Numeric \code{T} x \code{T} covariance matrix.
#' @param threshold Numeric (> 1). Multiplicative progression threshold.
#'
#' @return Numeric vector of length \code{T} giving \eqn{S_T(t)} for
#'   \eqn{t = 1, \ldots, T}.
#'
#' @importFrom mvtnorm pmvnorm
#'
S_tumour_from_params <- function(mu, Sigma, threshold) {
  T_visits <- length(mu)
  log_c <- log(threshold)
  S <- numeric(T_visits)
  for (t in 1:T_visits) {
    A_t     <- build_pairwise_A(t)
    mu_t    <- as.numeric(A_t %*% mu[1:t])
    Sigma_t <- A_t %*% Sigma[1:t, 1:t, drop = FALSE] %*% t(A_t)
    Sigma_t <- (Sigma_t + t(Sigma_t)) / 2
    S[t] <- as.numeric(pmvnorm(lower = rep(-Inf, nrow(A_t)),
                               upper = rep(log_c, nrow(A_t)),
                               mean  = mu_t, sigma = Sigma_t))
  }
  S
}


#' Lesion-side survival from logistic-regression parameters
#'
#' Computes \eqn{S_D(t)}, the probability of no new lesion by visit
#' \eqn{t}, given fitted logistic-regression coefficients \code{beta} and
#' the model design implied by \code{fit}. The per-visit non-event
#' probability is averaged across patients to obtain a population-level
#' survival via cumulative product.
#'
#' @param beta Numeric vector of logistic-regression coefficients.
#' @param fit A fitted \code{\link[stats]{glm}} object (used only to
#'   reconstruct the design matrix via \code{model.matrix} on its
#'   formula).
#' @param lesion_out A list with components \code{events} (lesion event
#'   matrix) and \code{all_tumour_sizes} (per-visit tumour sizes,
#'   including baseline), e.g. from \code{\link{lesion_events}}.
#' @param treatment Numeric vector of length \code{nrow(events)} giving
#'   the patient-level treatment indicator.
#'
#' @return Numeric vector of length \code{T} giving \eqn{S_D(t)} for
#'   \eqn{t = 1, \ldots, T}.
#'
#' @importFrom dplyr arrange group_by mutate ungroup summarise pull "%>%"
#' @importFrom stats model.matrix plogis delete.response terms
#'
S_lesion_from_params <- function(beta, fit, lesion_out, treatment) {
  D          <- lesion_out$events
  tumour_lag <- lesion_out$all_tumour_sizes
  n          <- nrow(D)
  T_visits   <- ncol(D)

  df <- data.frame(
    id         = rep(1:n, each = T_visits),
    time       = rep(1:T_visits, times = n),
    D          = as.vector(t(D)),
    tumour_lag = as.vector(t(tumour_lag[, 1:T_visits])),
    treatment  = rep(treatment, each = T_visits)
  )

  X <- model.matrix(delete.response(terms(fit)), data = df)
  eta <- as.numeric(X %*% beta)
  df$prob_no <- 1 - plogis(eta)

  # Bind variables to silence R CMD check NOTEs
  id <- time <- prob_no <- S_indiv <- NULL

  df %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(S_indiv = cumprod(prob_no)) %>%
    ungroup() %>%
    group_by(time) %>%
    summarise(S = mean(S_indiv), .groups = "drop") %>%
    pull(S)
}


#' Fit the lesion-side longitudinal logistic regression
#'
#' Fits a longitudinal logistic regression of the lesion event indicator
#' on \code{factor(time)}, optionally \code{treatment}, and the lagged
#' tumour size. Patients with \code{NA} events (post-censoring rows) are
#' dropped; treatment is dropped from the right-hand side when only one
#' arm is present.
#'
#' @param lesion_out A list with components \code{events} and
#'   \code{all_tumour_sizes}.
#' @param treatment Numeric vector of length \code{nrow(events)} giving
#'   the patient-level treatment indicator.
#'
#' @return A fitted \code{\link[stats]{glm}} object.
#'
#' @importFrom stats glm binomial as.formula
#'
fit_lesion <- function(lesion_out, treatment) {
  D          <- lesion_out$events
  tumour_lag <- lesion_out$all_tumour_sizes
  n          <- nrow(D)
  T_visits   <- ncol(D)
  df <- data.frame(
    id         = rep(1:n, each = T_visits),
    time       = rep(1:T_visits, times = n),
    D          = as.vector(t(D)),
    tumour_lag = as.vector(t(tumour_lag[, 1:T_visits])),
    treatment  = rep(treatment, each = T_visits)
  )
  df_obs <- df[!is.na(df$D), ]
  use_trt <- length(unique(df_obs$treatment)) > 1
  rhs <- if (use_trt) "factor(time) + treatment + tumour_lag"
  else         "factor(time) + tumour_lag"
  glm(as.formula(paste("D ~", rhs)),
      family = binomial(), data = df_obs)
}


#' MVN-based estimator of tumour-side survival
#'
#' Convenience wrapper that fits the MVN mean and covariance to the
#' observed log-ratio matrix and then evaluates \code{\link{S_tumour_from_params}}.
#'
#' @param Y Numeric \code{n} x \code{T} matrix of log-ratios versus
#'   baseline.
#' @param threshold Numeric (> 1). Multiplicative progression threshold.
#'
#' @return Numeric vector of length \code{T} giving the estimated
#'   tumour-side survival on the visit grid.
#'
#' @importFrom stats cov
#'
estimate_S_tumour_mvn <- function(Y, threshold) {
  mu_hat    <- colMeans(Y)
  Sigma_hat <- cov(Y)
  S_tumour_from_params(mu_hat, Sigma_hat, threshold)
}


#' Logistic-regression estimator of lesion-side survival
#'
#' Convenience wrapper that fits the lesion-side longitudinal logistic
#' regression with \code{\link{fit_lesion}} and evaluates
#' \code{\link{S_lesion_from_params}} at the fitted coefficients.
#'
#' @param lesion_out A list with components \code{events} and
#'   \code{all_tumour_sizes}.
#' @param treatment Numeric vector of length \code{nrow(events)}.
#'
#' @return Numeric vector of length \code{T} giving the estimated
#'   lesion-side survival on the visit grid.
#'
#' @importFrom stats coef
#'
estimate_S_lesion <- function(lesion_out, treatment) {
  fit <- fit_lesion(lesion_out, treatment)
  S_lesion_from_params(coef(fit), fit, lesion_out, treatment)
}


#' Single-dataset AugBin PFS point estimate
#'
#' Computes the augmented-binary progression-free survival curve on a
#' single simulated dataset by combining the tumour-side multivariate-normal
#' survival with the lesion-side logistic-regression survival:
#' \deqn{S_{\mathrm{PFS}}(t) = S_T(t) \, S_D(t).}
#' The product form assumes independence between tumour and lesion event
#' processes; relaxing that assumption is the motivation for the
#' copula-based estimator in \code{\link{pfs_copula}}. Implements the
#' method of Lin and Wason (2020).
#'
#' @param tumour A list as returned by \code{\link{tumour_events}}, with at
#'   least element \code{log_ratios_vs_baseline}.
#' @param lesion A list as returned by \code{\link{lesion_events}}, with
#'   elements \code{events} and \code{all_tumour_sizes}.
#' @param treatment Numeric vector of length
#'   \code{nrow(tumour$log_ratios_vs_baseline)} giving the per-patient
#'   treatment indicator (typically 0 for control, 1 for treatment).
#' @param threshold Numeric (> 1). Multiplicative tumour-progression
#'   threshold. Defaults to 1.2 (RECIST v1.1).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{S_pfs}}{Numeric vector of length \code{T} giving the AugBin
#'     PFS estimate.}
#'   \item{\code{S_tumour}}{Tumour-side survival \eqn{S_T(t)}.}
#'   \item{\code{S_lesion}}{Lesion-side survival \eqn{S_D(t)}.}
#' }
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(50, mu = mu, covariance = Sigma, threshold = 1.2)
#' les   <- lesion_events(50, tum$log_ratios_vs_baseline,
#'                        alpha = -2.5, beta = 0, gamma = 0.2,
#'                        treatment_arm = 0)
#' fit <- augbin_estimate(tum, les, treatment = rep(0, 50))
#' round(fit$S_pfs, 3)
#'
#' @seealso \code{\link{augbin_bootstrap}} for confidence intervals,
#'   \code{\link{augbin_test_pooled}} for two-arm testing,
#'   \code{\link{pfs_copula}} for the copula alternative.
#'
#' @references
#' Lin, C.-J. and Wason, J. (2020). Efficient analysis of time-to-event
#' endpoints when the event involves a continuous variable crossing a
#' threshold. \emph{Journal of Statistical Planning and Inference}, 208,
#' 119--129.
#'
#' @export
augbin_estimate <- function(tumour, lesion, treatment, threshold = 1.2) {
  stopifnot(
    is.list(tumour), "log_ratios_vs_baseline" %in% names(tumour),
    is.list(lesion), all(c("events", "all_tumour_sizes") %in% names(lesion)),
    length(treatment) == nrow(tumour$log_ratios_vs_baseline),
    nrow(tumour$log_ratios_vs_baseline) == nrow(lesion$events),
    threshold > 1
  )

  S_t <- estimate_S_tumour_mvn(tumour$log_ratios_vs_baseline, threshold)
  S_d <- estimate_S_lesion(lesion, treatment)

  list(S_pfs = S_t * S_d, S_tumour = S_t, S_lesion = S_d)
}


#' Patient-level bootstrap confidence intervals for AugBin PFS
#'
#' Computes pointwise patient-level percentile bootstrap confidence
#' intervals for the AugBin PFS estimator. Each replicate resamples
#' patients with replacement, refits the tumour-side MVN and lesion-side
#' logistic models, and combines them via \eqn{S_T(t) S_D(t)}. Replicates
#' that produce degenerate fits (e.g. when the resampled cohort lacks
#' variability in lesion outcomes) are recorded as \code{NA} rows and
#' excluded from the quantile computation.
#'
#' @param tumour A list as returned by \code{\link{tumour_events}}.
#' @param lesion A list as returned by \code{\link{lesion_events}}.
#' @param treatment Numeric vector of length
#'   \code{nrow(tumour$log_ratios_vs_baseline)} giving per-patient
#'   treatment indicators.
#' @param B Integer (>= 1). Number of bootstrap replicates.
#' @param threshold Numeric (> 1). Tumour progression threshold. Defaults
#'   to 1.2.
#' @param conf_level Numeric in (0, 1). Confidence level for the
#'   percentile intervals. Defaults to 0.95.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{S_boot}}{A \code{B} x \code{T} matrix where row \code{b}
#'     holds the AugBin PFS curve from the \code{b}-th bootstrap replicate.
#'     Rows are \code{NA} for replicates that failed.}
#'   \item{\code{lower, upper}}{Numeric vectors of length \code{T} giving
#'     the \code{conf_level} percentile confidence interval at each visit.}
#'   \item{\code{n_failed}}{Integer count of bootstrap replicates that
#'     produced degenerate fits.}
#' }
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(30, mu = mu, covariance = Sigma, threshold = 1.2)
#' les   <- lesion_events(30, tum$log_ratios_vs_baseline,
#'                        alpha = -2.5, beta = 0, gamma = 0.2,
#'                        treatment_arm = 0)
#'
#' fit <- augbin_bootstrap(tum, les, treatment = rep(0, 30), B = 20)
#' round(rbind(fit$lower, fit$upper), 3)
#'
#' \donttest{
#' # Closer to the simulation-study defaults (B = 200)
#' fit_full <- augbin_bootstrap(tum, les, treatment = rep(0, 30), B = 200)
#' }
#'
#' @seealso \code{\link{augbin_estimate}}, \code{\link{pfs_copula_boot}}
#'
#' @importFrom stats quantile
#' @export
augbin_bootstrap <- function(tumour, lesion, treatment, B,
                             threshold = 1.2, conf_level = 0.95) {
  stopifnot(
    is.list(tumour), "log_ratios_vs_baseline" %in% names(tumour),
    is.list(lesion), all(c("events", "all_tumour_sizes") %in% names(lesion)),
    length(treatment) == nrow(tumour$log_ratios_vs_baseline),
    B >= 1,
    conf_level > 0, conf_level < 1,
    threshold > 1
  )

  n         <- nrow(tumour$log_ratios_vs_baseline)
  T_visits  <- ncol(tumour$log_ratios_vs_baseline)
  alpha_err <- 1 - conf_level
  S_boot    <- matrix(NA_real_, B, T_visits)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)

    tumour_b <- list(
      log_ratios_vs_baseline = tumour$log_ratios_vs_baseline[idx, , drop = FALSE]
    )
    lesion_b <- list(
      all_tumour_sizes = lesion$all_tumour_sizes[idx, , drop = FALSE],
      events           = lesion$events[idx, , drop = FALSE]
    )

    S_boot[b, ] <- tryCatch(
      augbin_estimate(tumour_b, lesion_b, treatment[idx], threshold)$S_pfs,
      error = function(e) rep(NA_real_, T_visits)
    )
  }

  lower <- apply(S_boot, 2, stats::quantile, probs = alpha_err / 2,
                 na.rm = TRUE)
  upper <- apply(S_boot, 2, stats::quantile, probs = 1 - alpha_err / 2,
                 na.rm = TRUE)

  n_failed <- sum(is.na(S_boot[, 1]))
  if (n_failed > 0.1 * B) {
    warning(sprintf(
      "%d of %d bootstrap replicates failed (%.1f%%); CIs may be unreliable.",
      n_failed, B, 100 * n_failed / B
    ))
  }

  list(
    S_boot   = S_boot,
    lower    = lower,
    upper    = upper,
    n_failed = n_failed
  )
}

# -----------------------------------------------------------------------------
# vech helpers and asymptotic covariance of vech(Sigma_hat)
# -----------------------------------------------------------------------------

#' Vectorise the lower triangle of a symmetric matrix
#'
#' Returns the lower-triangular elements (including the diagonal) of a
#' symmetric matrix as a vector, in column-major order. Used by the
#' delta-method variance machinery for the AugBin tumour-side fit.
#'
#' @param M A symmetric numeric matrix.
#'
#' @return A numeric vector of length \eqn{T(T+1)/2}.
#'
vech <- function(M) M[lower.tri(M, diag = TRUE)]


#' Index pairs for vech ordering
#'
#' Returns the (row, column) index pairs corresponding to the entries of
#' \code{\link{vech}} for a \code{T_visits} x \code{T_visits} symmetric
#' matrix.
#'
#' @param T_visits Integer.
#'
#' @return An integer matrix with two columns (\code{row}, \code{col}) and
#'   \eqn{T(T+1)/2} rows.
#'
#' @export
vech_indices <- function(T_visits) {
  which(lower.tri(matrix(0, T_visits, T_visits), diag = TRUE), arr.ind = TRUE)
}


#' Asymptotic covariance of vech(Sigma_hat) under normality
#'
#' Computes the asymptotic covariance matrix of \code{vech(Sigma_hat)} for
#' a sample of size \code{n} drawn from a multivariate normal with
#' covariance \code{Sigma}, using
#' \eqn{\mathrm{Cov}(s_{ij}, s_{kl}) = (\sigma_{ik}\sigma_{jl} + \sigma_{il}\sigma_{jk})/(n-1)}.
#'
#' @param Sigma A \code{T_visits} x \code{T_visits} covariance matrix.
#' @param n Integer sample size.
#'
#' @return A numeric matrix of dimension \eqn{T(T+1)/2}.
#'
asymp_cov_vech_Sigma <- function(Sigma, n) {
  T_visits <- nrow(Sigma)
  idx <- vech_indices(T_visits)
  K <- nrow(idx)
  V <- matrix(0, K, K)
  for (a in 1:K) {
    i <- idx[a, 1]; j <- idx[a, 2]
    for (b in 1:K) {
      k <- idx[b, 1]; l <- idx[b, 2]
      V[a, b] <- (Sigma[i, k] * Sigma[j, l] + Sigma[i, l] * Sigma[j, k]) / (n - 1)
    }
  }
  V
}


#' Per-arm AugBin delta-method evaluator (point estimate and gradients)
#'
#' Fits the AugBin model on one arm and returns the PFS curve together
#' with finite-difference gradients of each \eqn{S_{\mathrm{PFS}}(t)}
#' with respect to the tumour-side mean, the vech of the tumour-side
#' covariance, and the lesion-side regression coefficients, alongside the
#' parameter-block covariances. Used by the two-arm pooled-hazard test in
#' \code{\link{augbin_test_pooled}}.
#'
#' @param tumour A list with \code{log_ratios_vs_baseline}.
#' @param lesion A list with \code{events} and \code{all_tumour_sizes}.
#' @param treatment Numeric vector of length \code{n}.
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#' @param eps Numeric. Finite-difference step. Defaults to 1e-3.
#'
#' @return A list with elements \code{S_pfs}, \code{grad_mu},
#'   \code{grad_vechSig}, \code{grad_beta}, \code{V_mu}, \code{V_vechSig},
#'   \code{V_beta}.
#'
#' @importFrom stats coef vcov cov
augbin_fit_with_grads <- function(tumour, lesion, treatment, threshold = 1.2,
                                  eps = 1e-3) {
  Y        <- tumour$log_ratios_vs_baseline
  n        <- nrow(Y)
  T_visits <- ncol(Y)

  mu_hat    <- colMeans(Y)
  Sigma_hat <- cov(Y)
  V_mu      <- Sigma_hat / n
  V_vechSig <- asymp_cov_vech_Sigma(Sigma_hat, n)

  fit       <- fit_lesion(lesion, treatment)
  beta_hat  <- coef(fit)
  V_beta    <- vcov(fit)
  p_beta    <- length(beta_hat)

  S_T   <- S_tumour_from_params(mu_hat, Sigma_hat, threshold)
  S_D   <- S_lesion_from_params(beta_hat, fit, lesion, treatment)
  S_pfs <- S_T * S_D

  # d S_pfs(t) / d mu_i
  grad_mu <- matrix(0, T_visits, T_visits)
  for (i in seq_along(mu_hat)) {
    mu_pert <- mu_hat
    mu_pert[i] <- mu_pert[i] + eps
    S_T_pert <- S_tumour_from_params(mu_pert, Sigma_hat, threshold)
    grad_mu[i, ] <- (S_T_pert * S_D - S_pfs) / eps
  }

  # d S_pfs(t) / d vech(Sigma)_k
  vech_Sigma <- vech(Sigma_hat)
  K   <- length(vech_Sigma)
  idx <- vech_indices(T_visits)
  grad_vechSig <- matrix(0, K, T_visits)
  for (k in 1:K) {
    i <- idx[k, 1]; j <- idx[k, 2]
    Sigma_pert <- Sigma_hat
    if (i == j) {
      Sigma_pert[i, i] <- Sigma_pert[i, i] + eps
    } else {
      Sigma_pert[i, j] <- Sigma_pert[i, j] + eps
      Sigma_pert[j, i] <- Sigma_pert[j, i] + eps
    }
    S_T_pert <- tryCatch(S_tumour_from_params(mu_hat, Sigma_pert, threshold),
                         error = function(e) S_T)
    grad_vechSig[k, ] <- (S_T_pert * S_D - S_pfs) / eps
  }

  # d S_pfs(t) / d beta_i
  grad_beta <- matrix(0, p_beta, T_visits)
  for (i in seq_along(beta_hat)) {
    beta_pert <- beta_hat
    beta_pert[i] <- beta_pert[i] + eps
    S_D_pert <- S_lesion_from_params(beta_pert, fit, lesion, treatment)
    grad_beta[i, ] <- (S_T * S_D_pert - S_pfs) / eps
  }

  list(S_pfs        = S_pfs,
       grad_mu      = grad_mu,
       grad_vechSig = grad_vechSig,
       grad_beta    = grad_beta,
       V_mu         = V_mu,
       V_vechSig    = V_vechSig,
       V_beta       = V_beta)
}


#' Per-arm AugBin delta-method standard errors for S_PFS(t)
#'
#' Internal helper. Given a fitted AugBin model and its gradient blocks
#' from \code{\link{augbin_fit_with_grads}}, returns the per-visit PFS
#' estimates together with delta-method standard errors. Used by
#' \code{\link{augbin_test_delta}} and \code{\link{augbin_test_pooled}}.
#'
#' @param tumour,lesion,treatment,threshold As in \code{\link{augbin_estimate}}.
#' @param fd_step Numeric finite-difference step for the gradient blocks.
#'   Defaults to \code{1e-3}.
#'
#' @return A list with \code{S_pfs} (numeric vector of length \code{T}) and
#'   \code{se} (numeric vector of length \code{T} of delta-method standard
#'   errors).
#'
#' @keywords internal
augbin_delta_full <- function(tumour, lesion, treatment, threshold = 1.2,
                              fd_step = 1e-3) {
  fit      <- augbin_fit_with_grads(tumour, lesion, treatment, threshold, fd_step)
  T_visits <- length(fit$S_pfs)
  var_S    <- numeric(T_visits)

  qf <- function(g, V) as.numeric(crossprod(g, V %*% g))

  for (t in seq_len(T_visits)) {
    var_S[t] <- qf(fit$grad_mu[, t],       fit$V_mu) +
      qf(fit$grad_vechSig[, t],  fit$V_vechSig) +
      qf(fit$grad_beta[, t],     fit$V_beta)
  }

  if (any(var_S < -1e-8)) {
    warning("AugBin delta variance had materially negative entries (min = ",
            signif(min(var_S), 3), "); clipped to 0.")
  }

  # Clip at 0 to guard against tiny negative variances from numerical error
  # in finite-difference gradients
  se_S <- sqrt(pmax(var_S, 0))

  list(S_pfs = fit$S_pfs, se = se_S)
}

#' AugBin delta-method Wald test at the final visit
#'
#' Two-arm Wald test on the log-ratio of treatment-to-control PFS at the
#' final visit, with arm-specific variances obtained from
#' \code{\link{augbin_delta_full}}. The test statistic is
#' \deqn{Z = \log S_0(T) - \log S_1(T),
#'   \quad
#'   \mathrm{Var}(Z) = \frac{\mathrm{Var}(S_0(T))}{S_0(T)^2} +
#'                     \frac{\mathrm{Var}(S_1(T))}{S_1(T)^2},}
#' and under the null hypothesis of equal final-visit survival,
#' \eqn{Z^2 / \mathrm{Var}(Z) \sim \chi^2_1}.
#'
#' This test compares survival only at the final visit. For a test that
#' aggregates information across all visits via a pooled hazard contrast,
#' see \code{\link{augbin_test_pooled}}.
#'
#' @param trial A two-arm trial list as returned by
#'   \code{\link{simulate_trial}}.
#' @param threshold Numeric (> 1). Tumour progression threshold. Defaults
#'   to 1.2.
#' @param eps Numeric in (0, 0.5). Boundary tolerance: if either arm's
#'   final-visit survival estimate falls within \code{eps} of 0 or 1, the
#'   test returns \code{NA_real_}. Defaults to \code{1e-6}.
#'
#' @return A scalar two-sided p-value, or \code{NA_real_} on numerical
#'   failure (degenerate arm fit, boundary survival estimate, or
#'   non-positive variance).
#'
#' @examples
#' set.seed(1)
#' mu_ctrl <- c(-0.10, -0.30, -0.46, -0.50, -0.55)
#' mu_trt  <- c(-0.20, -0.40, -0.56, -0.60, -0.65)
#' Sigma <- matrix(0.05, 5, 5); diag(Sigma) <- c(0.05, 0.10, 0.14, 0.16, 0.18)
#'
#' trial <- simulate_trial(80, mu_ctrl, mu_trt, Sigma,
#'                         beta_lesion = -0.5, gamma_lesion = 0.2)
#' augbin_test_delta(trial)
#'
#' @seealso \code{\link{augbin_test_pooled}}, \code{\link{augbin_estimate}},
#'   \code{\link{km_logrank_pvalue}}
#'
#' @importFrom stats pchisq
#' @export
augbin_test_delta <- function(trial, threshold = 1.2, eps = 1e-6) {
  stopifnot(
    is.list(trial),
    all(c("ctrl", "trt") %in% names(trial)),
    threshold > 1,
    eps > 0, eps < 0.5
  )

  treat_c <- rep(0L, nrow(trial$ctrl$tumour$log_ratios_vs_baseline))
  treat_t <- rep(1L, nrow(trial$trt$tumour$log_ratios_vs_baseline))

  ctrl_res <- tryCatch(
    augbin_delta_full(trial$ctrl$tumour, trial$ctrl$lesion, treat_c, threshold),
    error = function(e) {
      warning("AugBin delta fit failed on control arm: ", conditionMessage(e))
      NULL
    }
  )
  trt_res <- tryCatch(
    augbin_delta_full(trial$trt$tumour, trial$trt$lesion, treat_t, threshold),
    error = function(e) {
      warning("AugBin delta fit failed on treatment arm: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ctrl_res) || is.null(trt_res)) return(NA_real_)

  final  <- length(ctrl_res$S_pfs)
  S0     <- ctrl_res$S_pfs[final]
  S1     <- trt_res$S_pfs[final]

  if (S0 < eps || S1 < eps || S0 > 1 - eps || S1 > 1 - eps) {
    return(NA_real_)
  }

  var_S0 <- ctrl_res$se[final]^2
  var_S1 <- trt_res$se[final]^2

  Z     <- log(S0) - log(S1)
  var_Z <- var_S0 / S0^2 + var_S1 / S1^2

  if (!is.finite(var_Z) || var_Z <= 0) return(NA_real_)

  stats::pchisq(Z^2 / var_Z, df = 1, lower.tail = FALSE)
}


#' AugBin pooled-hazard test with delta-method variance
#'
#' Two-arm test on the pooled discrete-hazard contrast
#' \deqn{Z(\tau) = \tfrac{1}{2} \sum_{t=1}^{\tau} \bigl(h_1(t) - h_0(t)\bigr),}
#' with variance propagated by the delta method through
#' \eqn{\theta \to S \to h \to Z}. Under the null hypothesis of equal
#' arm-specific PFS curves, \eqn{Z^2 / \mathrm{Var}(Z) \sim \chi^2_1}.
#'
#' The variance combines three parameter blocks per arm (the tumour-side
#' mean, the vech of the tumour-side covariance, and the lesion-side
#' regression coefficients), each weighted by the chain-rule gradient
#' \eqn{dZ/d\theta = (dS/d\theta)(dZ/dS)} with \eqn{dZ/dS} obtained
#' analytically. Arms are assumed asymptotically independent, which holds
#' under random assignment but not under matched-pair designs.
#'
#' @param trial A two-arm trial list as returned by
#'   \code{\link{simulate_trial}}.
#' @param threshold Numeric (> 1). Tumour progression threshold. Defaults
#'   to 1.2.
#'
#' @return A scalar two-sided p-value, or \code{NA_real_} on numerical
#'   failure (degenerate arm fit, non-positive variance, or no overlapping
#'   at-risk visits between arms).
#'
#' @examples
#' set.seed(1)
#' mu_ctrl <- c(-0.10, -0.30, -0.46, -0.50, -0.55)
#' mu_trt  <- c(-0.20, -0.40, -0.56, -0.60, -0.65)
#' Sigma <- matrix(0.05, 5, 5); diag(Sigma) <- c(0.05, 0.10, 0.14, 0.16, 0.18)
#'
#' trial <- simulate_trial(80, mu_ctrl, mu_trt, Sigma,
#'                         beta_lesion = -0.5, gamma_lesion = 0.2)
#' augbin_test_pooled(trial)
#'
#' @seealso \code{\link{augbin_test_delta}}, \code{\link{km_logrank_pvalue}},
#'   \code{\link{copula_L1_permutation_test}}
#'
#' @references
#' Lin, C.-J. and Wason, J. (2020). Efficient analysis of time-to-event
#' endpoints when the event involves a continuous variable crossing a
#' threshold. \emph{Journal of Statistical Planning and Inference}, 208,
#' 119--129.
#'
#' @importFrom stats pchisq
#' @export
augbin_test_pooled <- function(trial, threshold = 1.2) {
  stopifnot(
    is.list(trial),
    all(c("ctrl", "trt") %in% names(trial)),
    threshold > 1
  )

  treat_c <- rep(0L, nrow(trial$ctrl$tumour$log_ratios_vs_baseline))
  treat_t <- rep(1L, nrow(trial$trt$tumour$log_ratios_vs_baseline))

  ctrl <- tryCatch(
    augbin_fit_with_grads(trial$ctrl$tumour, trial$ctrl$lesion,
                          treat_c, threshold),
    error = function(e) {
      warning("AugBin pooled fit failed on control arm: ", conditionMessage(e))
      NULL
    }
  )
  trt <- tryCatch(
    augbin_fit_with_grads(trial$trt$tumour, trial$trt$lesion,
                          treat_t, threshold),
    error = function(e) {
      warning("AugBin pooled fit failed on treatment arm: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ctrl) || is.null(trt)) return(NA_real_)

  S0 <- ctrl$S_pfs
  S1 <- trt$S_pfs
  pz  <- pooled_hazard_Z(S0, S1)
  Z   <- pz$Z
  tau <- pz$tau

  if (!is.finite(Z) || tau < 1) {
    warning("AugBin pooled test: no overlapping at-risk visits; returning NA.")
    return(NA_real_)
  }

  # Per-arm contribution to Var(Z) under the delta method
  arm_var <- function(fit, g_S) {
    dZ_dmu      <- as.numeric(fit$grad_mu      %*% g_S)
    dZ_dvechSig <- as.numeric(fit$grad_vechSig %*% g_S)
    dZ_dbeta    <- as.numeric(fit$grad_beta    %*% g_S)
    drop(crossprod(dZ_dmu,      fit$V_mu      %*% dZ_dmu)) +
      drop(crossprod(dZ_dvechSig, fit$V_vechSig %*% dZ_dvechSig)) +
      drop(crossprod(dZ_dbeta,    fit$V_beta    %*% dZ_dbeta))
  }

  g_S0 <- gradZ_wrt_S_arm(S0, tau, arm_sign = -1)
  g_S1 <- gradZ_wrt_S_arm(S1, tau, arm_sign = +1)

  varZ <- arm_var(ctrl, g_S0) + arm_var(trt, g_S1)

  if (!is.finite(varZ) || varZ <= 0) return(NA_real_)

  stats::pchisq(Z^2 / varZ, df = 1, lower.tail = FALSE)
}
