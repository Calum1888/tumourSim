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
#' on \eqn{(0,1)^2} by maximum likelihood.
#'
#' @param u,v Numeric vectors of pseudo-observations in \eqn{(0, 1)}, of
#'   equal length.
#' @param family Character. One of \code{"clayton"}, \code{"frank"}, or
#'   \code{"gumbel"}. Defaults to \code{"clayton"}.
#'
#' @return A list with elements \code{theta} (estimated dependence parameter),
#'   \code{copula} (fitted \code{copula} object), \code{family}, and
#'   \code{aic}.
#'
#' @importFrom copula claytonCopula frankCopula gumbelCopula fitCopula
#' @importFrom stats AIC coef
#' @export
fit_copula <- function(u, v, family = c("clayton", "frank", "gumbel")) {
  family <- match.arg(family)
  cop <- switch(family,
                clayton = claytonCopula(dim = 2),
                frank   = frankCopula(dim = 2),
                gumbel  = gumbelCopula(dim = 2))
  fit <- fitCopula(cop, cbind(u, v), method = "ml")
  list(theta = coef(fit), copula = fit@copula, family = family, aic = AIC(fit))
}


#' Copula-based estimator of progression-free survival
#'
#' Estimates PFS on a user-supplied time grid by combining Kaplan-Meier
#' marginal survival functions for the two component endpoints (typically
#' new lesion and tumour progression) with a fitted bivariate copula.
#'
#' @param d_time,d_status Vectors of times and (0/1) event indicators for
#'   the first endpoint (e.g. new lesion).
#' @param y_time,y_status Vectors of times and (0/1) event indicators for
#'   the second endpoint (e.g. tumour progression).
#' @param grid Numeric vector of time points at which to evaluate PFS.
#' @param family Character copula family passed to \code{\link{fit_copula}}.
#'   Defaults to \code{"gumbel"}.
#'
#' @return A list with elements \code{time} (the input grid), \code{pfs}
#'   (estimated survival at each grid point), \code{theta}, \code{aic}, and
#'   \code{S_D}, \code{S_Y} (marginal survival evaluated on the grid).
#'
#' @details Marginal survivals are estimated by Kaplan-Meier and read off
#'   as right-continuous step functions. The copula is fit by maximum
#'   likelihood on pseudo-observations using only patients with both events
#'   observed, then PFS is reconstructed via the inclusion-exclusion
#'   identity.
#'
#' @importFrom survival survfit Surv
#' @importFrom copula pCopula
#' @importFrom stats stepfun
#' @export
pfs_copula <- function(d_time, d_status, y_time, y_status,
                       grid, family = "gumbel") {

  fit_D <- survfit(Surv(d_time, d_status) ~ 1)
  fit_Y <- survfit(Surv(y_time, y_status) ~ 1)
  S_D   <- stepfun(fit_D$time, c(1, fit_D$surv))
  S_Y   <- stepfun(fit_Y$time, c(1, fit_Y$surv))

  u <- km_pseudo(d_time, d_status)
  v <- km_pseudo(y_time, y_status)
  both_obs <- d_status == 1 & y_status == 1
  cop_fit  <- fit_copula(u[both_obs], v[both_obs], family = family)

  s_d <- S_D(grid)
  s_y <- S_Y(grid)
  cop_val <- pCopula(cbind(1 - s_d, 1 - s_y), cop_fit$copula)
  pfs <- s_d + s_y - 1 + cop_val

  list(time = grid, pfs = pfs, theta = cop_fit$theta,
       aic = cop_fit$aic, S_D = s_d, S_Y = s_y)
}


#' Bootstrap confidence intervals for the copula-based PFS estimator
#'
#' Computes pointwise bootstrap confidence intervals for the copula-based
#' PFS estimator from \code{\link{pfs_copula}} using the complementary
#' log-log transform.
#'
#' @inheritParams pfs_copula
#' @param B Integer. Number of bootstrap replicates. Defaults to 200.
#' @param alpha Numeric in \eqn{(0, 1)}. Nominal two-sided error rate;
#'   resulting CIs have nominal level \eqn{1 - \alpha}. Defaults to 0.05.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{pfs},
#'   \code{lower}, \code{upper}.
#'
#' @details Bootstrap resamples are drawn with replacement at the patient
#'   level. The standard error is computed on the complementary log-log
#'   scale and back-transformed to preserve the (0, 1) range of the
#'   survival function. Failed bootstrap replicates are dropped silently.
#'
#' @importFrom stats qnorm sd
#' @export
pfs_copula_boot <- function(d_time, d_status, y_time, y_status,
                            grid, family = "gumbel", B = 200,
                            alpha = 0.05) {
  n <- length(d_time)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(grid))

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    res <- try(pfs_copula(d_time[idx], d_status[idx],
                          y_time[idx], y_status[idx],
                          grid = grid, family = family),
               silent = TRUE)
    if (!inherits(res, "try-error")) boot_mat[b, ] <- res$pfs
  }

  point <- pfs_copula(d_time, d_status, y_time, y_status, grid, family)

  eps <- 1e-8
  S_hat   <- pmin(pmax(point$pfs, eps), 1 - eps)
  cll_hat <- log(-log(S_hat))

  cll_boot <- log(-log(pmin(pmax(boot_mat, eps), 1 - eps)))
  se_cll   <- apply(cll_boot, 2, sd, na.rm = TRUE)

  z <- qnorm(1 - alpha / 2)
  cll_lo <- cll_hat - z * se_cll
  cll_hi <- cll_hat + z * se_cll

  lower <- exp(-exp(cll_hi))
  upper <- exp(-exp(cll_lo))

  data.frame(time = grid, pfs = point$pfs, lower = lower, upper = upper)
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
#' @export
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
#'   indicators (0/1) for one arm.
#' @param grid Numeric vector of evaluation times.
#' @param family Character copula family.
#' @param B Integer. Number of bootstrap replicates.
#' @param cens_at Numeric. Administrative censoring cutoff; defaults to
#'   \code{Inf} (no censoring).
#'
#' @return A list with elements \code{pfs} (point estimate) and \code{se}
#'   (bootstrap standard error), each of length \code{length(grid)};
#'   \code{NULL} if the fit could not be obtained (e.g. fewer than two
#'   events in either marginal).
#'
#' @export
copula_arm <- function(t_events, l_events, grid, family, B, cens_at = Inf) {

  t_time <- find_event_time(t_events)
  l_time <- find_event_time(l_events)

  t_time <- apply_admin_censoring(t_time, cens_at)
  l_time <- apply_admin_censoring(l_time, cens_at)

  if (sum(t_time$status) < 2 || sum(l_time$status) < 2) return(NULL)

  pt <- try(pfs_copula(l_time$time, l_time$status,
                       t_time$time, t_time$status,
                       grid = grid, family = family),
            silent = TRUE)
  if (inherits(pt, "try-error")) return(NULL)

  n <- length(t_time$time)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(grid))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    res <- try(pfs_copula(l_time$time[idx], l_time$status[idx],
                          t_time$time[idx], t_time$status[idx],
                          grid = grid, family = family),
               silent = TRUE)
    if (!inherits(res, "try-error")) boot_mat[b, ] <- res$pfs
  }

  se <- apply(boot_mat, 2, sd, na.rm = TRUE)
  list(pfs = pt$pfs, se = se)
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
#' @export
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
