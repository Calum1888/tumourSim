# =============================================================================
# Kaplan-Meier helpers and the log-rank baseline test
# =============================================================================


#' Probability-integral transform via Kaplan-Meier marginals
#'
#' Computes pseudo-observations on the open unit interval by evaluating the
#' Kaplan-Meier estimate of the marginal cumulative distribution function
#' \eqn{F(t) = 1 - S(t)} at each observed time. These are marginal inputs to
#' bivariate copula fits in the two-step estimation scheme of Shih and Louis
#' (1995), used internally by \code{\link{fit_copula}} and
#' \code{\link{pfs_copula}}.
#'
#' These are pseudo-observations in the copula sense of Genest, Ghoudi, and
#' Rivest (1995) - rank-based uniformised marginals for a parametric copula
#' likelihood - not the jackknife pseudo-observations of Andersen, Klein,
#' and Rosthoj.
#'
#' @param time Numeric vector of observed event or censoring times.
#' @param status Integer vector (0/1) of event indicators of the same length
#'   as \code{time}; 1 indicates an event, 0 censoring.
#' @param eps Numeric in (0, 0.5). Truncation tolerance; returned values are
#'   clipped to \eqn{[\code{eps}, 1 - \code{eps}]} to avoid boundary
#'   problems in copula likelihoods. Defaults to \code{1e-6}.
#'
#' @return A numeric vector of length \code{length(time)} with values in
#'   \eqn{[\code{eps}, 1 - \code{eps}]}.
#'
#' @details When all observations are censored, every returned value
#'   collapses to \code{eps}; the downstream copula fit will be
#'   uninformative.
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(50, mu = mu, covariance = Sigma, threshold = 1.2)
#' les   <- lesion_events(50, tum$log_ratios_vs_baseline,
#'                        alpha = -2.5, beta = 0, gamma = 0.2,
#'                        treatment_arm = 0)
#' pfs <- combined_event(tum$events, les$events)
#' u   <- km_pseudo(pfs$time, pfs$status)
#' summary(u)
#'
#' @references
#' Shih, J. H. and Louis, T. A. (1995). Inferences on the association
#' parameter in copula models for bivariate survival data. \emph{Biometrics},
#' 51(4), 1384--1399.
#'
#' Genest, C., Ghoudi, K., and Rivest, L.-P. (1995). A semiparametric
#' estimation procedure of dependence parameters in multivariate families of
#' distributions. \emph{Biometrika}, 82(3), 543--552.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun
#' @export
km_pseudo <- function(time, status, eps = 1e-6) {
  stopifnot(
    length(time) == length(status),
    all(status %in% c(0, 1)),
    all(time >= 0),
    eps > 0, eps < 0.5
  )

  fit   <- survival::survfit(survival::Surv(time, status) ~ 1)
  S_hat <- stats::stepfun(fit$time, c(1, fit$surv))

  pmin(pmax(1 - S_hat(time), eps), 1 - eps)
}


#' Two-arm log-rank test on the composite PFS endpoint
#'
#' Computes the log-rank p-value for a treatment-vs-control comparison on
#' the composite first-event time built from tumour progression and new
#' lesion events. Provides the Kaplan-Meier baseline against which the
#' copula-based permutation test is benchmarked in
#' \code{\link{run_all_scenarios}}.
#'
#' @param trial A nested list as returned by \code{\link{simulate_trial}},
#'   with \code{ctrl} and \code{trt} components each containing
#'   \code{tumour$events} and \code{lesion$events}.
#'
#' @return A scalar two-sided p-value for the null hypothesis that the
#'   composite PFS curves are equal between arms, computed from the log-rank
#'   chi-squared statistic on 1 degree of freedom via
#'   \code{\link[survival]{survdiff}}.
#'
#' @examples
#' set.seed(1)
#' mu_ctrl <- c(-0.10, -0.30, -0.46, -0.50, -0.55)
#' mu_trt  <- c(-0.20, -0.40, -0.56, -0.60, -0.65)
#' Sigma <- matrix(0.05, 5, 5); diag(Sigma) <- c(0.05, 0.10, 0.14, 0.16, 0.18)
#'
#' trial <- simulate_trial(50, mu_ctrl, mu_trt, Sigma,
#'                         beta_lesion = -0.5, gamma_lesion = 0.2)
#' km_logrank_pvalue(trial)
#'
#' @seealso \code{\link{copula_L1_permutation_test}},
#'   \code{\link{simulate_trial}}
#'
#' @importFrom survival survdiff Surv
#' @importFrom stats pchisq
#' @export
km_logrank_pvalue <- function(trial) {
  ce_c <- combined_event(trial$ctrl$tumour$events, trial$ctrl$lesion$events)
  ce_t <- combined_event(trial$trt$tumour$events,  trial$trt$lesion$events)

  pooled <- rbind(
    cbind(ce_c, arm = 0L),
    cbind(ce_t, arm = 1L)
  )

  surv_obj <- survival::Surv(pooled$time, pooled$status)
  lr <- survival::survdiff(surv_obj ~ pooled$arm)

  stats::pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
}

#' Monte Carlo approximation of the true progression-free survival curve
#'
#' Approximates the true PFS curve under the data-generating model by
#' simulating a large single-arm cohort and applying Kaplan-Meier to the
#' composite tumour-plus-lesion endpoint. The result is used as ground
#' truth for coverage and bias evaluation in
#' \code{\link{simulate_metrics_all}}. Monte Carlo error decays as
#' \eqn{O(1/\sqrt{n_{\mathrm{big}}})}, so for accuracy at the third decimal
#' place use \code{n_big} of order \eqn{10^5} or larger.
#'
#' @param grid Numeric vector of non-negative time points at which to
#'   evaluate the true PFS curve.
#' @param mu Numeric vector of length \code{T} giving the tumour log-ratio
#'   mean, passed to \code{\link{tumour_events}}.
#' @param covariance Numeric \code{T} x \code{T} covariance matrix for the
#'   tumour trajectory, passed to \code{\link{tumour_events}}.
#' @param threshold Numeric (> 1). Multiplicative tumour-progression
#'   threshold, passed to \code{\link{tumour_events}}.
#' @param n_big Integer. Size of the simulated cohort. Defaults to
#'   \code{1e5}.
#' @param alpha,beta,gamma,treatment_arm Lesion-model coefficients passed
#'   to \code{\link{lesion_events}}. Defaults give the single-arm baseline
#'   used in the simulation studies of Regan (2026).
#'
#' @return A numeric vector of length \code{length(grid)} giving the Monte
#'   Carlo approximation of true PFS at each grid point. Grid points beyond
#'   the largest simulated event time take the value of the survival
#'   function at the rightmost observed time (\code{extend = TRUE}).
#'
#' @examples
#' \donttest{
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#'
#' # Small n_big for example speed; use 1e5 in practice
#' set.seed(1)
#' s <- true_pfs(grid = 1:5, mu = mu, covariance = Sigma,
#'               threshold = 1.2, n_big = 2000)
#' round(s, 3)
#' }
#'
#' @seealso \code{\link{tumour_events}}, \code{\link{lesion_events}},
#'   \code{\link{combined_event}}, \code{\link{simulate_metrics_all}}
#'
#' @importFrom survival survfit Surv
#' @export
true_pfs <- function(grid, mu, covariance, threshold,
                     n_big = 1e5,
                     alpha = -2.5, beta = 0.0, gamma = 0.2,
                     treatment_arm = 0) {
  stopifnot(
    is.numeric(grid), all(grid >= 0),
    n_big >= 100
  )

  t_big <- tumour_events(n_big, mu = mu, covariance = covariance,
                         threshold = threshold)
  l_big <- lesion_events(n_big, t_big$log_ratios_vs_baseline,
                         alpha = alpha, beta = beta, gamma = gamma,
                         treatment_arm = treatment_arm)
  ce_big <- combined_event(t_big$events, l_big$events)

  surv_obj <- survival::Surv(ce_big$time, ce_big$status)
  km_big   <- survival::survfit(surv_obj ~ 1)

  summary(km_big, times = grid, extend = TRUE)$surv
}
