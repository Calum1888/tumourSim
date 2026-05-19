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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
estimate_S_lesion <- function(lesion_out, treatment) {
  fit <- fit_lesion(lesion_out, treatment)
  S_lesion_from_params(coef(fit), fit, lesion_out, treatment)
}


#' Single-dataset AugBin PFS point estimate
#'
#' Computes the AugBin PFS curve on one simulated dataset by combining the
#' tumour-side MVN survival with the lesion-side logistic-regression
#' survival: \eqn{S_{\mathrm{PFS}}(t) = S_T(t) S_D(t)}.
#'
#' @param tumour A list as returned by \code{\link{tumour_events}}.
#' @param lesion A list as returned by \code{\link{lesion_events}}.
#' @param treatment Numeric vector of length \code{nrow(events)} giving
#'   per-patient treatment indicator.
#' @param threshold Numeric (> 1). Tumour progression threshold. Defaults
#'   to 1.2.
#'
#' @return A list with elements \code{S_pfs}, \code{S_tumour}, and
#'   \code{S_lesion}, each a numeric vector of length \code{T}.
#'
#' @export
augbin_estimate <- function(tumour, lesion, treatment, threshold = 1.2) {
  S_t <- estimate_S_tumour_mvn(tumour$log_ratios_vs_baseline, threshold)
  S_d <- estimate_S_lesion(lesion, treatment)
  list(S_pfs = S_t * S_d, S_tumour = S_t, S_lesion = S_d)
}


#' Patient-level bootstrap confidence intervals for AugBin PFS
#'
#' Computes pointwise patient-level bootstrap confidence intervals for the
#' AugBin PFS estimator using the percentile method.
#'
#' @param tumour A list with at least \code{log_ratios_vs_baseline}.
#' @param lesion A list with \code{all_tumour_sizes} and \code{events}.
#' @param treatment Numeric vector of length \code{n}.
#' @param B Integer. Number of bootstrap replicates.
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#' @param alpha Numeric in \eqn{(0, 1)}. Nominal two-sided error rate.
#'   Defaults to 0.05.
#'
#' @return A list with elements \code{S_boot} (the \code{B x T} matrix of
#'   bootstrap PFS curves), \code{lower}, and \code{upper} (numeric
#'   vectors of length \code{T} giving the percentile CIs).
#'
#' @details Bootstrap replicates that throw errors (e.g. degenerate fits)
#'   contribute \code{NA_real_} rows and are excluded from the quantile
#'   computation via \code{na.rm = TRUE}.
#'
#' @importFrom stats quantile
#' @export
augbin_bootstrap <- function(tumour, lesion, treatment, B,
                             threshold = 1.2, alpha = 0.05) {
  n        <- nrow(tumour$log_ratios_vs_baseline)
  T_visits <- ncol(tumour$log_ratios_vs_baseline)
  S_boot   <- matrix(NA_real_, B, T_visits)

  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    tumour_b <- list(
      log_ratios_vs_baseline = tumour$log_ratios_vs_baseline[idx, , drop = FALSE]
    )
    lesion_b <- list(
      all_tumour_sizes = lesion$all_tumour_sizes[idx, , drop = FALSE],
      events           = lesion$events[idx, , drop = FALSE]
    )
    trt_b <- treatment[idx]

    S_boot[b, ] <- tryCatch(
      augbin_estimate(tumour_b, lesion_b, trt_b, threshold)$S_pfs,
      error = function(e) rep(NA_real_, T_visits)
    )
  }

  lower <- apply(S_boot, 2, quantile, probs = alpha / 2,     na.rm = TRUE)
  upper <- apply(S_boot, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  list(S_boot = S_boot, lower = lower, upper = upper)
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
#' @export
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
#' @export
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
#' @export
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
#' Computes the AugBin PFS curve and its pointwise delta-method standard
#' errors on a single arm by combining the parameter-block gradients and
#' covariances returned by \code{\link{augbin_fit_with_grads}}.
#'
#' @inheritParams augbin_fit_with_grads
#'
#' @return A list with elements \code{S_pfs} and \code{se}, each of
#'   length \code{T}.
#'
#' @export
augbin_delta_full <- function(tumour, lesion, treatment, threshold = 1.2,
                              eps = 1e-3) {
  fit <- augbin_fit_with_grads(tumour, lesion, treatment, threshold, eps)
  T_visits <- length(fit$S_pfs)

  var_S <- numeric(T_visits)
  for (t in 1:T_visits) {
    g_mu    <- fit$grad_mu[, t]
    g_vechS <- fit$grad_vechSig[, t]
    g_beta  <- fit$grad_beta[, t]
    var_S[t] <- as.numeric(t(g_mu)    %*% fit$V_mu      %*% g_mu)    +
      as.numeric(t(g_vechS) %*% fit$V_vechSig %*% g_vechS) +
      as.numeric(t(g_beta)  %*% fit$V_beta    %*% g_beta)
  }
  se_S <- sqrt(pmax(var_S, 0))
  list(S_pfs = fit$S_pfs, se = se_S)
}


#' AugBin delta-method Wald test at the final visit
#'
#' Two-arm Wald test on the log-ratio of treatment-to-control PFS at the
#' final visit:
#' \eqn{Z = 0.5 (\log S_0(T) - \log S_1(T))},
#' \eqn{\mathrm{Var}(Z) = 0.25 (\mathrm{Var}(S_0)/S_0^2 + \mathrm{Var}(S_1)/S_1^2)},
#' with arm-specific variances from \code{\link{augbin_delta_full}}.
#' Under H_0, \eqn{Z^2 / \mathrm{Var}(Z) \sim \chi^2_1}.
#'
#' @param trial A two-arm trial list (e.g. from \code{\link{simulate_trial}}).
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#'
#' @return A scalar p-value, or \code{NA_real_} on numerical failure.
#'
#' @importFrom stats pchisq
#' @export
augbin_test_delta <- function(trial, threshold = 1.2) {

  treat_c <- rep(0, nrow(trial$ctrl$tumour$log_ratios_vs_baseline))
  treat_t <- rep(0, nrow(trial$trt$tumour$log_ratios_vs_baseline))

  ctrl_res <- tryCatch(
    augbin_delta_full(trial$ctrl$tumour, trial$ctrl$lesion, treat_c, threshold),
    error = function(e) NULL
  )
  trt_res <- tryCatch(
    augbin_delta_full(trial$trt$tumour, trial$trt$lesion, treat_t, threshold),
    error = function(e) NULL
  )
  if (is.null(ctrl_res) || is.null(trt_res)) return(NA_real_)

  T_idx  <- length(ctrl_res$S_pfs)
  S0     <- ctrl_res$S_pfs[T_idx]
  S1     <- trt_res$S_pfs[T_idx]
  var_S0 <- ctrl_res$se[T_idx]^2
  var_S1 <- trt_res$se[T_idx]^2

  eps <- 1e-6
  S0c <- pmin(pmax(S0, eps), 1 - eps)
  S1c <- pmin(pmax(S1, eps), 1 - eps)

  Z     <- 0.5 * (log(S0c) - log(S1c))
  var_Z <- 0.25 * (var_S0 / S0c^2 + var_S1 / S1c^2)

  if (!is.finite(var_Z) || var_Z <= 0) return(NA_real_)
  chisq <- Z^2 / var_Z
  pchisq(chisq, df = 1, lower.tail = FALSE)
}


#' AugBin pooled-hazard test with delta-method variance
#'
#' Two-arm test on the pooled discrete-hazard contrast
#' \eqn{Z(\tau) = 0.5 \sum_{t=1}^{\tau} (h_1(t) - h_0(t))} with variance
#' propagated by the delta method through
#' \eqn{\theta \to S \to h \to Z}. Under H_0,
#' \eqn{Z^2 / \mathrm{Var}(Z) \sim \chi^2_1}.
#'
#' @param trial A two-arm trial list (e.g. from \code{\link{simulate_trial}}).
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#'
#' @return A scalar p-value, or \code{NA_real_} on numerical failure.
#'
#' @details The variance combines the three parameter blocks per arm
#'   (mean, vech of covariance, lesion coefficients), each weighted by
#'   the gradient \eqn{dZ/d\theta = (dS/d\theta) \cdot (dZ/dS)} where
#'   \eqn{dZ/dS} is given analytically by \code{\link{gradZ_wrt_S_arm}}.
#'   Arms are assumed asymptotically independent.
#'
#' @importFrom stats pchisq
#' @export
augbin_test_pooled <- function(trial, threshold = 1.2) {

  treat_c <- rep(0, nrow(trial$ctrl$tumour$log_ratios_vs_baseline))
  treat_t <- rep(1, nrow(trial$trt$tumour$log_ratios_vs_baseline))

  ctrl <- tryCatch(
    augbin_fit_with_grads(trial$ctrl$tumour, trial$ctrl$lesion,
                          treat_c, threshold),
    error = function(e) NULL
  )
  trt  <- tryCatch(
    augbin_fit_with_grads(trial$trt$tumour, trial$trt$lesion,
                          treat_t, threshold),
    error = function(e) NULL
  )
  if (is.null(ctrl) || is.null(trt)) return(NA_real_)

  S0 <- ctrl$S_pfs
  S1 <- trt$S_pfs
  pz <- pooled_hazard_Z(S0, S1)
  Z  <- pz$Z
  tau <- pz$tau
  if (!is.finite(Z) || tau < 1) return(NA_real_)

  g_S0 <- gradZ_wrt_S_arm(S0, tau, arm_sign = -1)
  g_S1 <- gradZ_wrt_S_arm(S1, tau, arm_sign = +1)

  dZ_dmu_0      <- as.numeric(ctrl$grad_mu      %*% g_S0)
  dZ_dvechSig_0 <- as.numeric(ctrl$grad_vechSig %*% g_S0)
  dZ_dbeta_0    <- as.numeric(ctrl$grad_beta    %*% g_S0)

  dZ_dmu_1      <- as.numeric(trt$grad_mu       %*% g_S1)
  dZ_dvechSig_1 <- as.numeric(trt$grad_vechSig  %*% g_S1)
  dZ_dbeta_1    <- as.numeric(trt$grad_beta     %*% g_S1)

  varZ <-
    drop(t(dZ_dmu_0)      %*% ctrl$V_mu      %*% dZ_dmu_0)      +
    drop(t(dZ_dvechSig_0) %*% ctrl$V_vechSig %*% dZ_dvechSig_0) +
    drop(t(dZ_dbeta_0)    %*% ctrl$V_beta    %*% dZ_dbeta_0)    +
    drop(t(dZ_dmu_1)      %*% trt$V_mu       %*% dZ_dmu_1)      +
    drop(t(dZ_dvechSig_1) %*% trt$V_vechSig  %*% dZ_dvechSig_1) +
    drop(t(dZ_dbeta_1)    %*% trt$V_beta     %*% dZ_dbeta_1)

  if (!is.finite(varZ) || varZ <= 0) return(NA_real_)
  chisq <- Z^2 / varZ
  pchisq(chisq, df = 1, lower.tail = FALSE)
}
