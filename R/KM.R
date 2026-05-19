# =============================================================================
# Kaplan-Meier helpers and the log-rank baseline test
# =============================================================================


#' Kaplan-Meier pseudo-observations of the marginal CDF
#'
#' Computes pseudo-observations on the (0, 1) scale by evaluating the
#' Kaplan-Meier estimate of the marginal CDF \eqn{F(t) = 1 - S(t)} at each
#' observed event or censoring time. These pseudo-observations are used as
#' marginal inputs to bivariate copula fits in \code{\link{pfs_copula}}.
#'
#' @param time Numeric vector of observed event or censoring times.
#' @param status Integer vector (0/1) of event indicators of the same length
#'   as \code{time}.
#'
#' @return A numeric vector of pseudo-observations, truncated to the open
#'   interval \eqn{[10^{-6}, 1 - 10^{-6}]} to avoid boundary issues in
#'   downstream copula fitting.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun
#' @export
km_pseudo <- function(time, status) {
  fit <- survfit(Surv(time, status) ~ 1)
  S_hat <- stepfun(fit$time, c(1, fit$surv))
  pmin(pmax(1 - S_hat(time), 1e-6), 1 - 1e-6)
}


#' Two-arm log-rank test on the composite PFS endpoint
#'
#' Computes the log-rank p-value for a treatment-vs-control comparison on
#' the composite first-event time built from tumour and lesion events.
#' Used as the Kaplan-Meier baseline in power studies.
#'
#' @param trial A nested list as returned by \code{\link{simulate_trial}},
#'   with \code{ctrl} and \code{trt} components each containing
#'   \code{tumour$events} and \code{lesion$events}.
#'
#' @return A scalar p-value from \code{\link[survival]{survdiff}} on the
#'   log-rank chi-squared statistic with 1 degree of freedom.
#'
#' @importFrom survival survdiff Surv
#' @importFrom stats pchisq
#' @export
km_logrank_pvalue <- function(trial) {
  ce_c <- first_event_times(trial$ctrl$tumour$events, trial$ctrl$lesion$events)
  ce_t <- first_event_times(trial$trt$tumour$events,  trial$trt$lesion$events)
  df <- rbind(cbind(ce_c, arm = 0L), cbind(ce_t, arm = 1L))
  sd <- survdiff(Surv(time, status) ~ arm, data = df)
  pchisq(sd$chisq, df = 1, lower.tail = FALSE)
}


#' Approximate the true PFS curve by large-sample Monte Carlo
#'
#' Approximates the true progression-free survival curve under the
#' data-generating model by simulating a very large single-arm cohort and
#' applying Kaplan-Meier to the composite endpoint. Used as ground truth
#' for coverage and bias evaluation in
#' \code{\link{simulate_metrics_all}}.
#'
#' @param grid Numeric vector of time points at which to evaluate the true
#'   PFS curve.
#' @param n_big Integer. Size of the simulated cohort. Defaults to
#'   \code{1e5}; the Monte Carlo error in the truth is
#'   \eqn{O(1/\sqrt{n_{\mathrm{big}}})}.
#' @param mean,covariance,threshold Parameters of the tumour model passed
#'   to \code{\link{tumour_events}}.
#' @param alpha,beta,gamma,treatment_arm Parameters of the lesion model
#'   passed to \code{\link{lesion_events}}. Defaults give the
#'   single-arm baseline used in the included scenarios.
#'
#' @return A numeric vector of length \code{length(grid)} giving Monte Carlo
#'   approximations of true PFS at each grid point.
#'
#' @importFrom survival survfit Surv
#' @export
true_pfs <- function(grid, n_big = 1e5,
                     mean, covariance, threshold,
                     alpha = -2.5, beta = 0.0, gamma = 0.2, treatment_arm = 0) {

  t_big <- tumour_events(n_big, mean, covariance, threshold)
  l_big <- lesion_events(n_big, t_big$log_ratios_vs_baseline,
                         alpha = alpha, beta = beta, gamma = gamma,
                         treatment_arm = treatment_arm)
  ce_big <- combined_event(l_big$events, t_big$events)
  km_big <- survfit(Surv(time, status) ~ 1, data = ce_big)
  summary(km_big, times = grid, extend = TRUE)$surv
}
