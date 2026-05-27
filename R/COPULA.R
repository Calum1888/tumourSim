# =============================================================================
# Copula-based PFS estimators and tests
# =============================================================================
# The copula approach reconstructs PFS as
#   S_PFS(t) = S_D(t) + S_Y(t) - 1 + C(1 - S_D(t), 1 - S_Y(t); theta),
# where S_D and S_Y are the marginal Kaplan-Meier survival functions for
# the lesion and tumour endpoints, and C is a fitted bivariate Archimedean
# copula linking the two failure times.
# =============================================================================


#' Fit a one-parameter Archimedean copula by maximum likelihood
#'
#' Fits a Clayton, Frank, or Gumbel copula to a pair of pseudo-observations
#' on \eqn{(0, 1)^2} by maximum likelihood. Pseudo-observations are typically
#' obtained from \code{\link{km_pseudo}} applied to right-censored marginal
#' time-to-event data.
#'
#' @param u,v Numeric vectors of pseudo-observations in \eqn{(0, 1)}, of
#'   equal length (at least 2).
#' @param family Character. One of \code{"gumbel"} (default; identified in
#'   Regan (2026) as the most appropriate copula for RECIST PFS endpoints),
#'   \code{"clayton"}, or \code{"frank"}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{theta}}{Numeric scalar. The estimated dependence
#'     parameter.}
#'   \item{\code{copula}}{The fitted S4 \code{copula} object, suitable for
#'     passing to \code{\link[copula]{pCopula}}.}
#'   \item{\code{family}}{Character. The fitted family (echoes the input).}
#'   \item{\code{aic}}{Numeric. Akaike information criterion.}
#'   \item{\code{loglik}}{Numeric. Maximised log-likelihood.}
#' }
#'
#' @examples
#' set.seed(1)
#' # Simulate dependent pseudo-observations from a Gumbel copula
#' true_cop <- copula::gumbelCopula(param = 1.5, dim = 2)
#' uv <- copula::rCopula(200, true_cop)
#'
#' fit <- fit_copula(uv[, 1], uv[, 2], family = "gumbel")
#' c(theta = fit$theta, aic = fit$aic)
#'
#' @seealso \code{\link{km_pseudo}}, \code{\link{pfs_copula}}
#'
#' @references
#' Hofert, M., Kojadinovic, I., Maechler, M., and Yan, J. (2018).
#' \emph{Elements of Copula Modeling with R}. Springer.
#'
#' @importFrom copula claytonCopula frankCopula gumbelCopula fitCopula
#' @importFrom stats AIC coef
#' @export
fit_copula <- function(u, v, family = c("gumbel", "clayton", "frank")) {
  family <- match.arg(family)
  stopifnot(
    is.numeric(u), is.numeric(v),
    length(u) == length(v),
    length(u) >= 2L,
    all(u > 0 & u < 1),
    all(v > 0 & v < 1)
  )

  cop <- switch(family,
                clayton = copula::claytonCopula(dim = 2),
                frank   = copula::frankCopula(dim = 2),
                gumbel  = copula::gumbelCopula(dim = 2))

  fit <- copula::fitCopula(cop, cbind(u, v), method = "ml")

  list(
    theta  = as.numeric(stats::coef(fit)),
    copula = fit@copula,
    family = family,
    aic    = stats::AIC(fit),
    loglik = fit@loglik
  )
}


#' Copula-based estimator of progression-free survival
#'
#' Estimates PFS on a user-supplied time grid by combining Kaplan-Meier
#' marginal survival functions for the two component endpoints (typically
#' new lesion and tumour progression) with a fitted bivariate Archimedean
#' copula, using the survival-copula identity
#' \deqn{S_{\mathrm{PFS}}(t) = S_D(t) + S_Y(t) - 1 +
#'       C\bigl(1 - S_D(t),\, 1 - S_Y(t);\,\theta\bigr)}
#' (Regan, 2026, eq. 4.5). Marginal survivals are estimated by Kaplan-Meier
#' and treated as right-continuous step functions. The copula is fit by
#' two-step maximum likelihood (Shih and Louis, 1995) on pseudo-observations
#' from \code{\link{km_pseudo}}.
#'
#' @param d_time,d_status Numeric and 0/1 vectors giving the observed times
#'   and event indicators for the lesion (D) endpoint.
#' @param y_time,y_status Numeric and 0/1 vectors giving the observed times
#'   and event indicators for the tumour-progression (Y) endpoint. Paired
#'   with \code{d_time} and \code{d_status} by patient (one row per
#'   patient).
#' @param grid Numeric vector of non-negative time points at which to
#'   evaluate the PFS curve.
#' @param family Character. Archimedean copula family, passed to
#'   \code{\link{fit_copula}}. Defaults to \code{"gumbel"}.
#' @param pseudo_obs Character. \code{"all"} (default) fits the copula on
#'   all patients via Kaplan-Meier pseudo-observations, matching the
#'   two-step framework of Shih and Louis (1995). \code{"doubly_observed"}
#'   restricts the fit to patients with both endpoint events observed and
#'   is included for backwards compatibility.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{time}}{The input \code{grid}.}
#'   \item{\code{pfs}}{Estimated PFS at each grid point.}
#'   \item{\code{theta}}{Fitted copula dependence parameter.}
#'   \item{\code{aic}}{AIC of the fitted copula.}
#'   \item{\code{S_D, S_Y}}{Marginal Kaplan-Meier survivals evaluated on
#'     \code{grid}.}
#' }
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(80, mu = mu, covariance = Sigma, threshold = 1.2)
#' les   <- lesion_events(80, tum$log_ratios_vs_baseline,
#'                        alpha = -2.5, beta = 0, gamma = 0.2,
#'                        treatment_arm = 0)
#' tum_te <- combined_event(tum$events, matrix(0, nrow = 80, ncol = 5))
#' les_te <- combined_event(matrix(0, nrow = 80, ncol = 5), les$events)
#'
#' fit <- pfs_copula(
#'   d_time = les_te$time, d_status = les_te$status,
#'   y_time = tum_te$time, y_status = tum_te$status,
#'   grid = 1:5
#' )
#' round(fit$pfs, 3)
#'
#' @seealso \code{\link{fit_copula}}, \code{\link{km_pseudo}},
#'   \code{\link{pfs_copula_boot}}, \code{\link{augbin_estimate}}
#'
#' @references
#' Shih, J. H. and Louis, T. A. (1995). Inferences on the association
#' parameter in copula models for bivariate survival data. \emph{Biometrics},
#' 51(4), 1384--1399.
#'
#' @importFrom survival survfit Surv
#' @importFrom copula pCopula
#' @importFrom stats stepfun
#' @export
pfs_copula <- function(d_time, d_status, y_time, y_status,
                       grid, family = "gumbel",
                       pseudo_obs = c("all", "doubly_observed")) {
  pseudo_obs <- match.arg(pseudo_obs)
  stopifnot(
    length(d_time) == length(d_status),
    length(y_time) == length(y_status),
    length(d_time) == length(y_time),
    all(d_status %in% c(0, 1)),
    all(y_status %in% c(0, 1)),
    all(d_time >= 0), all(y_time >= 0),
    is.numeric(grid), all(grid >= 0)
  )

  surv_D <- survival::Surv(d_time, d_status)
  fit_D  <- survival::survfit(surv_D ~ 1)
  S_D    <- stats::stepfun(fit_D$time, c(1, fit_D$surv))

  surv_Y <- survival::Surv(y_time, y_status)
  fit_Y  <- survival::survfit(surv_Y ~ 1)
  S_Y    <- stats::stepfun(fit_Y$time, c(1, fit_Y$surv))

  u <- km_pseudo(d_time, d_status)
  v <- km_pseudo(y_time, y_status)

  idx <- if (pseudo_obs == "doubly_observed") {
    which(d_status == 1 & y_status == 1)
  } else {
    seq_along(u)
  }
  copula_fit <- fit_copula(u[idx], v[idx], family = family)

  s_d <- S_D(grid)
  s_y <- S_Y(grid)

  # Sklar identity for the survival copula
  cop_val <- copula::pCopula(cbind(1 - s_d, 1 - s_y), copula_fit$copula)
  pfs     <- s_d + s_y - 1 + cop_val

  if (any(pfs < 0 | pfs > 1, na.rm = TRUE)) {
    warning("Estimated PFS outside [0, 1]; copula may be misspecified.")
  }

  list(
    time  = grid,
    pfs   = pfs,
    theta = copula_fit$theta,
    aic   = copula_fit$aic,
    S_D   = s_d,
    S_Y   = s_y
  )
}


#' Bootstrap confidence intervals for the copula-based PFS estimator
#'
#' Computes pointwise bootstrap confidence intervals for the copula-based
#' PFS estimator from \code{\link{pfs_copula}}, with a choice of two
#' interval methods.
#'
#' Bootstrap resamples patients with replacement and refits both the
#' marginal Kaplan-Meier estimators and the copula on each resample.
#' Replicates that fail (e.g. degenerate copula fits) are dropped via
#' \code{na.rm = TRUE}.
#'
#' Two methods for the confidence interval are provided:
#' \describe{
#'   \item{\code{"percentile"}}{(default) Pointwise intervals are formed
#'     by taking the \eqn{\alpha/2} and \eqn{1 - \alpha/2} quantiles of
#'     the bootstrap distribution at each grid point. This matches the
#'     textbook percentile bootstrap and is invariant to monotone
#'     transformations of the survival function.}
#'   \item{\code{"cll"}}{Wald-type interval on the complementary log-log
#'     scale: the bootstrap standard error of \eqn{\log(-\log S(t))} is
#'     used to form a symmetric interval on that scale, then
#'     back-transformed. This produces narrower intervals than the
#'     percentile method when the sampling distribution is approximately
#'     normal on the cloglog scale, and is the method used in the
#'     simulation results of Regan (2026).}
#' }
#'
#' @inheritParams pfs_copula
#' @param B Integer (>= 1). Number of bootstrap replicates. Defaults to 200.
#' @param conf_level Numeric in (0, 1). Confidence level. Defaults to 0.95.
#' @param method Character. Either \code{"percentile"} (default) or
#'   \code{"cll"}; see Details.
#'
#' @return A data frame with one row per grid point and columns
#'   \code{time}, \code{pfs} (the point estimate), \code{lower}, and
#'   \code{upper}.
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(60, mu = mu, covariance = Sigma, threshold = 1.2)
#' les   <- lesion_events(60, tum$log_ratios_vs_baseline,
#'                        alpha = -2.5, beta = 0, gamma = 0.2,
#'                        treatment_arm = 0)
#' tum_te <- combined_event(tum$events, matrix(0, nrow = 60, ncol = 5))
#' les_te <- combined_event(matrix(0, nrow = 60, ncol = 5), les$events)
#'
#' ci <- pfs_copula_boot(
#'   d_time = les_te$time, d_status = les_te$status,
#'   y_time = tum_te$time, y_status = tum_te$status,
#'   grid = 1:5, B = 20
#' )
#' ci
#'
#' \donttest{
#' # Closer to the thesis defaults (B = 200, cll method)
#' ci_thesis <- pfs_copula_boot(
#'   d_time = les_te$time, d_status = les_te$status,
#'   y_time = tum_te$time, y_status = tum_te$status,
#'   grid = 1:5, B = 200, method = "cll"
#' )
#' }
#'
#' @seealso \code{\link{pfs_copula}}, \code{\link{augbin_bootstrap}}
#'
#' @references
#' Regan, C. (2026). \emph{Copula Modelling for Mixed Continuous and Binary
#' Variables in Survival Analysis: Applications to Colorectal Cancer}.
#' MMath dissertation.
#'
#' @importFrom stats quantile sd qnorm
#' @export
pfs_copula_boot <- function(d_time, d_status, y_time, y_status,
                            grid, family = "gumbel", B = 200,
                            conf_level = 0.95,
                            method = c("percentile", "cll")) {
  method <- match.arg(method)
  stopifnot(
    B >= 1,
    conf_level > 0, conf_level < 1
  )

  # Compute point estimate first; fail fast if the original data is bad
  point <- pfs_copula(d_time, d_status, y_time, y_status, grid, family)

  n         <- length(d_time)
  alpha_err <- 1 - conf_level
  boot_mat  <- matrix(NA_real_, nrow = B, ncol = length(grid))

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    res <- tryCatch(
      pfs_copula(d_time[idx], d_status[idx],
                 y_time[idx], y_status[idx],
                 grid = grid, family = family),
      error = function(e) NULL
    )
    if (!is.null(res)) boot_mat[b, ] <- res$pfs
  }

  n_failed <- sum(is.na(boot_mat[, 1]))
  if (n_failed > 0.1 * B) {
    warning(sprintf(
      "%d of %d bootstrap replicates failed (%.1f%%); CIs may be unreliable.",
      n_failed, B, 100 * n_failed / B
    ))
  }

  if (method == "percentile") {
    lower <- apply(boot_mat, 2, stats::quantile,
                   probs = alpha_err / 2,     na.rm = TRUE)
    upper <- apply(boot_mat, 2, stats::quantile,
                   probs = 1 - alpha_err / 2, na.rm = TRUE)
  } else {
    # Complementary log-log Wald
    eps      <- 1e-8
    S_hat    <- pmin(pmax(point$pfs, eps), 1 - eps)
    cll_hat  <- log(-log(S_hat))
    cll_boot <- log(-log(pmin(pmax(boot_mat, eps), 1 - eps)))
    se_cll   <- apply(cll_boot, 2, stats::sd, na.rm = TRUE)
    z        <- stats::qnorm(1 - alpha_err / 2)
    # cloglog is decreasing in S, so the upper-cloglog bound gives the lower S bound
    lower <- exp(-exp(cll_hat + z * se_cll))
    upper <- exp(-exp(cll_hat - z * se_cll))
  }

  data.frame(
    time  = grid,
    pfs   = point$pfs,
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}


#' Per-arm copula PFS curve
#'
#' Computes a single-arm copula-based PFS curve on the integer visit grid
#' \code{1, ..., T_visits} given that arm's tumour and lesion outputs. If
#' the copula fit fails or fewer than five patients have both events
#' observed, falls back to the independence copula
#' \eqn{C(u, v) = u v}.
#'
#' @param arm_data A list with components \code{tumour$events} and
#'   \code{lesion$events}, e.g. one arm of \code{\link{simulate_trial}}.
#' @param T_visits Integer. Number of visits; the grid is \code{1:T_visits}.
#' @param family Character copula family. Defaults to \code{"gumbel"}.
#'
#' @return A numeric vector of length \code{T_visits} giving PFS at each
#'   visit.
#'
#' @importFrom survival survfit Surv
#' @importFrom copula claytonCopula frankCopula gumbelCopula fitCopula pCopula
#' @importFrom stats stepfun
#'
pfs_copula_arm_S <- function(arm_data, T_visits, family = "gumbel") {
  ll <- separate_event_times(arm_data$lesion$events)
  tt <- separate_event_times(arm_data$tumour$events)

  fit_D <- survfit(Surv(time, status) ~ 1, data = ll)
  fit_Y <- survfit(Surv(time, status) ~ 1, data = tt)
  S_D_fun <- stepfun(fit_D$time, c(1, fit_D$surv))
  S_Y_fun <- stepfun(fit_Y$time, c(1, fit_Y$surv))

  u <- km_pseudo(ll$time, ll$status)
  v <- km_pseudo(tt$time, tt$status)
  both_obs <- ll$status == 1 & tt$status == 1

  cop <- switch(family,
                clayton = claytonCopula(dim = 2),
                frank   = frankCopula(dim = 2),
                gumbel  = gumbelCopula(dim = 2))

  grid <- 1:T_visits
  s_d <- S_D_fun(grid); s_y <- S_Y_fun(grid)

  if (sum(both_obs) < 5) {
    return(s_d + s_y - 1 + (1 - s_d) * (1 - s_y))
  }

  fit_c <- tryCatch(suppressWarnings(
    fitCopula(cop, cbind(u[both_obs], v[both_obs]), method = "ml")
  ), error = function(e) NULL)

  if (is.null(fit_c)) {
    return(s_d + s_y - 1 + (1 - s_d) * (1 - s_y))
  }

  cop_val <- pCopula(cbind(1 - s_d, 1 - s_y), fit_c@copula)
  s_d + s_y - 1 + cop_val
}


#' Per-arm copula fit with bootstrap standard error
#'
#' Fits a copula-based PFS curve for a single arm and computes pointwise
#' bootstrap standard errors. Used by \code{\link{run_copula_power}} for
#' two-arm Wald tests on the survival difference at a chosen visit.
#'
#' @param t_events,l_events Numeric matrices of tumour and lesion event
#'   indicators (0/1) for one arm, of identical dimensions.
#' @param grid Numeric vector of non-negative evaluation times.
#' @param family Character copula family passed to \code{\link{fit_copula}}.
#'   Defaults to \code{"gumbel"}.
#' @param B Integer (>= 1). Number of bootstrap replicates. Defaults to 200.
#' @param cens_at Numeric. Administrative censoring cutoff; defaults to
#'   \code{Inf} (no censoring).
#'
#' @return A list with elements \code{pfs} (point estimate) and \code{se}
#'   (bootstrap standard error), each of length \code{length(grid)};
#'   \code{NULL} if the fit could not be obtained (fewer than two events
#'   on either margin, or the initial copula fit fails).
#'
#' @seealso \code{\link{pfs_copula}}, \code{\link{run_copula_power}}
#'
#' @importFrom stats sd
#' @export
copula_arm <- function(t_events, l_events, grid,
                       family = "gumbel", B = 200, cens_at = Inf) {
  stopifnot(
    identical(dim(t_events), dim(l_events)),
    is.numeric(grid), all(grid >= 0),
    B >= 1
  )

  t_time <- find_event_time(t_events)
  l_time <- find_event_time(l_events)
  t_time <- apply_admin_censoring(t_time, cens_at)
  l_time <- apply_admin_censoring(l_time, cens_at)

  if (sum(t_time$status) < 2 || sum(l_time$status) < 2) return(NULL)

  pt <- tryCatch(
    pfs_copula(l_time$time, l_time$status,
               t_time$time, t_time$status,
               grid = grid, family = family),
    error = function(e) NULL
  )
  if (is.null(pt)) return(NULL)

  n        <- length(t_time$time)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(grid))

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    res <- tryCatch(
      pfs_copula(l_time$time[idx], l_time$status[idx],
                 t_time$time[idx], t_time$status[idx],
                 grid = grid, family = family),
      error = function(e) NULL
    )
    if (!is.null(res)) boot_mat[b, ] <- res$pfs
  }

  se <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)

  list(
    pfs      = pt$pfs,
    se       = se,
    n_failed = sum(is.na(boot_mat[, 1])),
    n_boot   = B
  )
}

#' L1 mean-absolute-difference statistic between two survival curves
#'
#' Computes the L1 statistic
#' \eqn{L_1 = (1/b) \sum_{t=a}^{b} |S_1(t) - S_0(t)|}, where \eqn{a} is the
#' smallest event/censoring time in the pooled composite endpoint and
#' \eqn{b} is the smaller of the two arms' maximum follow-up times.
#'
#' @param S0,S1 Numeric vectors of length \code{T_visits} giving the
#'   estimated survival on the integer grid \code{1:T_visits} for control
#'   and treatment arms.
#' @param trial A nested trial list (e.g. from \code{\link{simulate_trial}})
#'   used to derive the pooled event-time range \eqn{(a, b)}.
#'
#' @return A scalar L1 statistic; \code{NA_real_} if \eqn{b < a}.
#'
L1_stat_from_S <- function(S0, S1, trial) {
  T_visits <- length(S0)
  ce_c <- first_event_times(trial$ctrl$tumour$events, trial$ctrl$lesion$events)
  ce_t <- first_event_times(trial$trt$tumour$events,  trial$trt$lesion$events)
  a <- max(1L, min(min(ce_c$time), min(ce_t$time)))
  b <- min(max(ce_c$time), max(ce_t$time))
  if (b < a) return(NA_real_)
  idx <- a:b
  sum(abs(S1[idx] - S0[idx])) / b
}


#' Observed copula-based L1 statistic for a two-arm trial
#'
#' Convenience wrapper that fits per-arm copula PFS curves via
#' \code{\link{pfs_copula_arm_S}} and computes the L1 mean-absolute-
#' difference statistic of \code{\link{L1_stat_from_S}}.
#'
#' @param trial A two-arm trial list (e.g. from \code{\link{simulate_trial}}).
#' @param T_visits Integer. Number of visits.
#' @param family Character copula family. Defaults to \code{"gumbel"}.
#'
#' @return A scalar L1 statistic, or \code{NA_real_} if a per-arm fit
#'   contained \code{NA}s.
#'
#' @export
copula_L1_stat <- function(trial, T_visits, family = "gumbel") {
  S0 <- pfs_copula_arm_S(trial$ctrl, T_visits, family)
  S1 <- pfs_copula_arm_S(trial$trt,  T_visits, family)
  if (any(is.na(S0)) || any(is.na(S1))) return(NA_real_)
  L1_stat_from_S(S0, S1, trial)
}


#' Permutation test for the copula L1 statistic
#'
#' One-sided permutation p-value for the L1 mean-absolute-difference
#' statistic between treatment and control survival curves. Arm labels
#' are shuffled without replacement, the copula is refit on each
#' permutation, and the p-value is the proportion of permuted L1*
#' statistics at least as large as the observed value.
#'
#' @param trial A two-arm trial list.
#' @param T_visits Integer. Number of visits.
#' @param B Integer. Number of permutations. Defaults to 200.
#' @param family Character copula family.
#' @param min_success_frac Numeric in \eqn{(0, 1]}. If the fraction of
#'   non-failing permutations is below this threshold, return
#'   \code{NA_real_}. Defaults to 0.2.
#'
#' @return A scalar permutation p-value, or \code{NA_real_} if the
#'   observed statistic was undefined or too few permutations succeeded.
#'
#' @details Under the null of identical joint distributions across arms,
#'   arm labels are exchangeable, so resampling without replacement
#'   generates valid null replicates. Because L1 is unsigned, the test
#'   has non-trivial power against crossings of the survival curves
#'   where signed tests cancel.
#'
#' @export
copula_L1_permutation_test <- function(trial, T_visits, B = 200,
                                       family = "gumbel",
                                       min_success_frac = 0.2) {
  L1_obs <- tryCatch(copula_L1_stat(trial, T_visits, family),
                     error = function(e) NA_real_)
  if (!is.finite(L1_obs)) return(NA_real_)

  n_c <- nrow(trial$ctrl$tumour$log_ratios_vs_baseline)
  n_t <- nrow(trial$trt$tumour$log_ratios_vs_baseline)
  n_total <- n_c + n_t

  pool <- list(
    Y_b   = rbind(trial$ctrl$tumour$log_ratios_vs_baseline,
                  trial$trt$tumour$log_ratios_vs_baseline),
    Y_n   = rbind(trial$ctrl$tumour$log_ratios_vs_nadir,
                  trial$trt$tumour$log_ratios_vs_nadir),
    tum_e = rbind(trial$ctrl$tumour$events,
                  trial$trt$tumour$events),
    les_e = rbind(trial$ctrl$lesion$events,
                  trial$trt$lesion$events),
    les_z = rbind(trial$ctrl$lesion$all_tumour_sizes,
                  trial$trt$lesion$all_tumour_sizes)
  )

  pick <- function(idx) list(
    tumour = list(log_ratios_vs_baseline = pool$Y_b[idx, , drop = FALSE],
                  log_ratios_vs_nadir    = pool$Y_n[idx, , drop = FALSE],
                  events                 = pool$tum_e[idx, , drop = FALSE]),
    lesion = list(events                 = pool$les_e[idx, , drop = FALSE],
                  all_tumour_sizes       = pool$les_z[idx, , drop = FALSE])
  )

  L1_perm <- numeric(B)
  for (b in seq_len(B)) {
    perm <- sample.int(n_total, n_total, replace = FALSE)
    fake_trial <- list(
      ctrl = pick(perm[1:n_c]),
      trt  = pick(perm[(n_c + 1):n_total])
    )
    L1_perm[b] <- tryCatch(copula_L1_stat(fake_trial, T_visits, family),
                           error = function(e) NA_real_)
  }
  n_ok <- sum(!is.na(L1_perm))
  if (n_ok < min_success_frac * B) return(NA_real_)

  mean(L1_perm >= L1_obs, na.rm = TRUE)
}


#' Observed max-difference statistic between two copula-based PFS curves
#'
#' Computes \eqn{\max_t |S_1(t) - S_0(t)|} over the integer visit grid,
#' with per-arm survival curves obtained from \code{\link{pfs_copula_arm_S}}.
#'
#' @param trial A two-arm trial list.
#' @param grid Numeric vector of evaluation visits. Defaults to \code{1:5}.
#' @param family Character copula family. Defaults to \code{"gumbel"}.
#'
#' @return A scalar max-difference statistic.
#'
#' @export
copula_delta_stat <- function(trial, grid = 1:5, family = "gumbel") {
  S_c <- pfs_copula_arm_S(trial$ctrl, length(grid), family)
  S_t <- pfs_copula_arm_S(trial$trt,  length(grid), family)
  max(abs(S_t - S_c))
}


#' Permutation test for the copula max-difference statistic
#'
#' One-sided permutation p-value for the max-difference statistic
#' \code{\link{copula_delta_stat}}. Arm labels are shuffled without
#' replacement and the copula is refit on each permutation.
#'
#' @param trial A two-arm trial list.
#' @param B Integer. Number of permutations. Defaults to 200.
#' @param grid Numeric vector of evaluation visits. Defaults to \code{1:5}.
#' @param family Character copula family.
#'
#' @return A scalar permutation p-value.
#'
#' @export
copula_test <- function(trial, B = 200, grid = 1:5, family = "gumbel") {

  Delta <- copula_delta_stat(trial, grid, family)

  n_c <- nrow(trial$ctrl$tumour$log_ratios_vs_baseline)
  n_t <- nrow(trial$trt$tumour$log_ratios_vs_baseline)
  n_total <- n_c + n_t

  pool <- list(
    Y     = rbind(trial$ctrl$tumour$log_ratios_vs_baseline,
                  trial$trt$tumour$log_ratios_vs_baseline),
    tum_e = rbind(trial$ctrl$tumour$events,
                  trial$trt$tumour$events),
    les_e = rbind(trial$ctrl$lesion$events,
                  trial$trt$lesion$events),
    les_z = rbind(trial$ctrl$lesion$all_tumour_sizes,
                  trial$trt$lesion$all_tumour_sizes)
  )

  Delta_b <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n_total, n_total, replace = FALSE)
    pseudo_ctrl_idx <- idx[1:n_c]
    pseudo_trt_idx  <- idx[(n_c + 1):n_total]
    fake_trial <- list(
      ctrl = list(
        tumour = list(log_ratios_vs_baseline = pool$Y[pseudo_ctrl_idx, , drop = FALSE],
                      events                 = pool$tum_e[pseudo_ctrl_idx, , drop = FALSE]),
        lesion = list(events                 = pool$les_e[pseudo_ctrl_idx, , drop = FALSE],
                      all_tumour_sizes       = pool$les_z[pseudo_ctrl_idx, , drop = FALSE])
      ),
      trt = list(
        tumour = list(log_ratios_vs_baseline = pool$Y[pseudo_trt_idx, , drop = FALSE],
                      events                 = pool$tum_e[pseudo_trt_idx, , drop = FALSE]),
        lesion = list(events                 = pool$les_e[pseudo_trt_idx, , drop = FALSE],
                      all_tumour_sizes       = pool$les_z[pseudo_trt_idx, , drop = FALSE])
      )
    )
    Delta_b[b] <- tryCatch(copula_delta_stat(fake_trial, grid, family),
                           error = function(e) NA_real_)
  }
  mean(abs(Delta_b) >= abs(Delta), na.rm = TRUE)
}


#' Patient-level bootstrap resampling of a single arm
#'
#' Resamples patients (rows) with replacement from one arm of a trial,
#' preserving all per-arm matrices in lock step (tumour log-ratios,
#' nadirs, events, lesion events, lesion tumour sizes). Useful for
#' patient-level bootstrap variance estimation that needs to keep
#' tumour and lesion components aligned within patients.
#'
#' @param arm_data One arm of a trial list, with \code{tumour} and
#'   \code{lesion} components.
#'
#' @return An arm list of the same shape with rows resampled with
#'   replacement.
#'
#' @export
resample_arm <- function(arm_data) {
  n <- nrow(arm_data$tumour$events)
  idx <- sample.int(n, n, replace = TRUE)
  list(
    tumour = list(
      log_ratios_vs_baseline = arm_data$tumour$log_ratios_vs_baseline[idx, , drop = FALSE],
      log_ratios_vs_nadir    = arm_data$tumour$log_ratios_vs_nadir[idx, , drop = FALSE],
      events                 = arm_data$tumour$events[idx, , drop = FALSE]
    ),
    lesion = list(
      all_tumour_sizes = arm_data$lesion$all_tumour_sizes[idx, , drop = FALSE],
      events           = arm_data$lesion$events[idx, , drop = FALSE]
    )
  )
}
