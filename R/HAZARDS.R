# =============================================================================
# Discrete-hazard helpers used by the pooled-hazard test
# =============================================================================
# Given a survival vector S(1), ..., S(T) (with S(0) = 1 by convention),
# these functions compute discrete hazards h(t) = 1 - S(t)/S(t-1),
# pooled-hazard contrast Z(tau) = 0.5 * sum_{t=1..tau} (h_1(t) - h_0(t)),
# and the gradient of Z with respect to an arm's survival vector (needed
# for the delta-method variance in augbin_test_pooled).
# =============================================================================


#' Largest visit at which both arms are still at risk
#'
#' Returns the largest visit \eqn{\tau} such that \eqn{S_0(t-1) > 0} and
#' \eqn{S_1(t-1) > 0} for all \eqn{t \le \tau}, i.e. the largest visit at
#' which the discrete hazards \eqn{h_r(t) = 1 - S_r(t)/S_r(t-1)} are
#' well-defined for both arms.
#'
#' @param S0,S1 Numeric vectors of length \code{T_visits} giving the
#'   estimated survival on the integer visit grid for control and
#'   treatment arms.
#'
#' @return Integer \eqn{\tau} in \eqn{\{0, 1, \ldots, T_{\mathrm{visits}}\}}.
#'   Returns 0 if either arm has \eqn{S_r(0) = 1} immediately followed by
#'   a zero survival at visit 1 (degenerate case).
#'
#' @export
choose_tau <- function(S0, S1) {
  T_visits <- length(S0)
  tau <- 0
  S0_prev <- 1
  S1_prev <- 1
  for (t in 1:T_visits) {
    if (S0_prev > 0 && S1_prev > 0) tau <- t else break
    S0_prev <- S0[t]
    S1_prev <- S1[t]
  }
  tau
}


#' Discrete hazards from a survival vector
#'
#' Computes the discrete hazards \eqn{h(t) = 1 - S(t)/S(t-1)} for
#' \eqn{t = 1, \ldots, T}, with \eqn{S(0) = 1} by convention.
#'
#' @param S Numeric vector of length \code{T_visits} giving the estimated
#'   survival on the integer visit grid.
#'
#' @return Numeric vector of discrete hazards, with \code{NA_real_} for
#'   any visit where the prior survival was already zero.
#'
#' @export
discrete_hazards <- function(S) {
  T_visits <- length(S)
  h <- numeric(T_visits)
  S_prev <- 1
  for (t in 1:T_visits) {
    if (S_prev <= 0) h[t] <- NA_real_
    else             h[t] <- 1 - S[t] / S_prev
    S_prev <- S[t]
  }
  h
}


#' Pooled discrete-hazard contrast \eqn{Z(\tau)}
#'
#' Computes the pooled-hazard test statistic
#' \eqn{Z(\tau) = 0.5 \sum_{t=1}^{\tau} (h_1(t) - h_0(t))},
#' where \eqn{h_r(t) = 1 - S_r(t)/S_r(t-1)} is the discrete hazard in arm
#' \eqn{r} and \eqn{\tau} is the largest visit at which both arms are
#' still at risk (\code{\link{choose_tau}}).
#'
#' @param S0,S1 Numeric vectors of length \code{T_visits} giving the
#'   estimated survival on the integer visit grid for control and
#'   treatment arms.
#'
#' @return A list with elements \code{Z} (the statistic, or
#'   \code{NA_real_} if any hazard is non-finite) and \code{tau}.
#'
#' @export
pooled_hazard_Z <- function(S0, S1) {
  tau <- choose_tau(S0, S1)
  if (tau < 1) return(list(Z = NA_real_, tau = tau))
  h0 <- discrete_hazards(S0)[1:tau]
  h1 <- discrete_hazards(S1)[1:tau]
  if (any(!is.finite(h0)) || any(!is.finite(h1))) {
    return(list(Z = NA_real_, tau = tau))
  }
  list(Z = 0.5 * sum(h1 - h0), tau = tau)
}


#' Gradient of \eqn{Z(\tau)} with respect to one arm's survival vector
#'
#' Returns \eqn{dZ/dS} for a single arm, holding the other arm fixed:
#' \eqn{Z = (\mathrm{sign}/2) \sum_{t=1}^{\tau} h_{\mathrm{arm}}(t)},
#' \eqn{h_{\mathrm{arm}}(t) = 1 - S(t)/S(t-1)}. The s-th entry of \eqn{S}
#' (for \eqn{1 \le s \le \tau}) enters in two places: as the numerator
#' of \eqn{h(s)} (contributes \eqn{-1/S(s-1)}) and as the denominator of
#' \eqn{h(s+1)} (contributes \eqn{S(s+1)/S(s)^2}, only when
#' \eqn{s + 1 \le \tau}).
#'
#' @param S Numeric vector of length \code{T_visits}; the arm's survival
#'   on the integer visit grid.
#' @param tau Integer in \eqn{\{1, \ldots, T_{\mathrm{visits}}\}} from
#'   \code{\link{choose_tau}}.
#' @param arm_sign Numeric, \code{+1} for the treatment arm and \code{-1}
#'   for the control arm (matching the sign of that arm's contribution to
#'   \eqn{Z}).
#'
#' @return Numeric vector of length \code{T_visits}; the gradient is zero
#'   in positions beyond \eqn{\tau}.
#'
#' @export
gradZ_wrt_S_arm <- function(S, tau, arm_sign) {
  T_visits <- length(S)
  g <- numeric(T_visits)
  # S_full[k] = S(k-1); S_full[1] = S(0) = 1
  S_full <- c(1, S)
  for (s in 1:tau) {
    term1 <- -1 / S_full[s]
    term2 <- if (s + 1 <= tau) S_full[s + 2] / S_full[s + 1]^2 else 0
    g[s]  <- 0.5 * arm_sign * (term1 + term2)
  }
  g
}
