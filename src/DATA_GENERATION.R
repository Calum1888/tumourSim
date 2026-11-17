# =============================================================================
# Data-generating functions
# =============================================================================
# These functions simulate the longitudinal tumour-size and binary new-lesion
# processes that underlie the composite PFS endpoint used throughout the
# package. They are intentionally simple and self-contained so that the
# estimation machinery in copula.R, augbin.R, and power.R can be tested
# against a known data-generating process.
# =============================================================================


#' Simulate longitudinal tumour progression events
#'
#' Simulates tumour progression events for a cohort of patients under the
#' RECIST framework. Log tumour sizes relative to baseline,
#' \eqn{Y_t = \log(z_t / z_0)} for visits \eqn{t = 1, \ldots, T}, are drawn
#' from a multivariate normal or (for heavy-tailed robustness studies) a
#' shifted multivariate-t distribution. A progression event is declared at
#' the first visit where \eqn{z_t} exceeds the running nadir
#' \eqn{\min(z_0, \ldots, z_{t-1})} by a multiplicative factor greater than
#' \code{threshold}; visits strictly after the first event are censored to
#' \code{NA}. The baseline log-ratio is fixed at zero by construction.
#'
#' @param n_patients Integer. Number of patients to simulate.
#' @param mu Numeric vector of length \code{T} giving the mean of the
#'   multivariate distribution of \eqn{Y_{1:T}}.
#' @param covariance Numeric symmetric positive-definite \code{T} x \code{T}
#'   covariance (or scale) matrix.
#' @param threshold Numeric (> 1). Multiplicative progression threshold on
#'   tumour size relative to the running nadir. The RECIST default
#'   corresponds to \code{threshold = 1.2}.
#' @param dist Character. Either \code{"mvnorm"} (default) for a multivariate
#'   normal log-ratio distribution, or \code{"mvt"} for a shifted
#'   multivariate-t distribution.
#' @param df Integer. Degrees of freedom for the multivariate-t distribution.
#'   Ignored when \code{dist = "mvnorm"}. Defaults to 4.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{log_ratios_vs_baseline}}{An \code{n_patients} x \code{T}
#'     matrix of simulated log-ratios \eqn{\log(z_t / z_0)}.}
#'   \item{\code{log_ratios_vs_nadir}}{An \code{n_patients} x \code{T} matrix
#'     of log-ratios versus the running prior nadir.}
#'   \item{\code{events}}{An \code{n_patients} x \code{T} matrix of event
#'     indicators (0/1), with \code{NA} for all visits strictly after the
#'     first event for each patient.}
#' }
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' sim   <- tumour_events(n_patients = 20, mu = mu, covariance = Sigma,
#'                        threshold = 1.2)
#' head(sim$events)
#'
#' # Heavy-tailed alternative for robustness studies
#' sim_t <- tumour_events(n_patients = 20, mu = mu, covariance = Sigma,
#'                        threshold = 1.2, dist = "mvt", df = 4)
#' colMeans(sim_t$events == 1, na.rm = TRUE)
#'
#' @seealso \code{\link{lesion_events}}, \code{\link{simulate_trial}}
#' @references
#' Eisenhauer, E. A. et al. (2009). New response evaluation criteria in solid
#' tumours: revised RECIST guideline (version 1.1). \emph{European Journal
#' of Cancer}, 45(2), 228--247. \doi{10.1016/j.ejca.2008.10.026}
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @export
tumour_events <- function(n_patients, mu, covariance, threshold,
                          dist = c("mvnorm", "mvt"), df = 4) {
  dist <- match.arg(dist)
  stopifnot(
    length(mu) == nrow(covariance),
    nrow(covariance) == ncol(covariance),
    threshold > 1,
    n_patients >= 1
  )

  Y <- if (dist == "mvt") {
    mvtnorm::rmvt(n = n_patients, sigma = covariance, df = df,
                  delta = mu, type = "shifted")
  } else {
    MASS::mvrnorm(n = n_patients, mu = mu, Sigma = covariance)
  }

  # Augment with baseline log-ratio (= 0) for the running-minimum computation
  log_ratio_with_baseline <- cbind(0, Y)

  # Running minimum of log(z_s / z_0) over s = 0, ..., t-1
  running_min_prior <- t(apply(log_ratio_with_baseline, 1, function(r) {
    head(cummin(r), -1)
  }))

  log_ratio_vs_nadir <- Y - running_min_prior
  event_indicators   <- (log_ratio_vs_nadir > log(threshold)) + 0L

  # Censor visits strictly after the first event
  cum_events   <- t(apply(event_indicators == 1, 1, cumsum))
  prior_events <- cum_events - (event_indicators == 1)
  event_indicators[prior_events > 0] <- NA

  list(
    log_ratios_vs_baseline = Y,
    log_ratios_vs_nadir    = log_ratio_vs_nadir,
    events                 = event_indicators
  )
}

#' Simulate new-lesion events conditional on tumour trajectory
#'
#' Simulates the occurrence of new lesions over time as Bernoulli events
#' whose probability at visit \eqn{t} depends on treatment arm and the lagged
#' (previous-visit) tumour size through a logistic or probit linear predictor
#' \deqn{\eta_{it} = \alpha_t + \beta_t \, \mathrm{arm}_i + \gamma_t \, z_{i,t-1},}
#' with \eqn{p_{it} = g^{-1}(\eta_{it})}. Setting \code{nonlinear = TRUE}
#' replaces \eqn{z_{i,t-1}} with \eqn{z_{i,t-1}^2}, which is useful for
#' studying robustness to model misspecification. Visits strictly after the
#' first lesion event are censored to \code{NA}.
#'
#' @param n_patients Integer. Number of patients to simulate.
#' @param log_tumour Numeric \code{n_patients} x \code{T} matrix of log-ratios
#'   of tumour size versus baseline, typically the
#'   \code{log_ratios_vs_baseline} element of \code{\link{tumour_events}}.
#' @param alpha Numeric scalar or vector of length \code{T}. Intercept of the
#'   linear predictor at each visit. A scalar is recycled to length \code{T};
#'   any other length must equal \code{T}.
#' @param beta Numeric scalar or vector of length \code{T}. Treatment-arm
#'   coefficient at each visit. Recycled as for \code{alpha}.
#' @param gamma Numeric scalar or vector of length \code{T}. Coefficient on
#'   the lagged tumour size at each visit. Recycled as for \code{alpha}.
#' @param treatment_arm Numeric scalar or vector of length \code{n_patients}.
#'   Treatment-arm indicator (typically 0/1). A scalar is recycled to length
#'   \code{n_patients}.
#' @param baseline Numeric vector of length \code{n_patients} giving baseline
#'   tumour sizes. Defaults to \code{runif(n_patients, 0, 1)}, matching the
#'   convention used in the simulation studies of Regan (2026).
#' @param link Character. Either \code{"logit"} (default) or \code{"probit"}.
#' @param nonlinear Logical. If \code{TRUE}, the tumour-size effect enters
#'   the linear predictor as \eqn{z_{i,t-1}^2} instead of \eqn{z_{i,t-1}}.
#'   Defaults to \code{FALSE}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{all_tumour_sizes}}{An \code{n_patients} x \code{(T+1)} matrix
#'     of tumour sizes on the original scale, with baseline in the first
#'     column.}
#'   \item{\code{lesion_probability}}{An \code{n_patients} x \code{T} matrix
#'     of per-visit lesion probabilities.}
#'   \item{\code{events}}{An \code{n_patients} x \code{T} matrix of lesion
#'     event indicators (0/1), with \code{NA} for visits strictly after the
#'     first event.}
#' }
#'
#' @examples
#' set.seed(1)
#' mu    <- c(0.00, 0.04, 0.07, 0.11, 0.14)
#' Sigma <- diag(c(0.25, 0.45, 0.50, 0.75, 1.00))
#' tum   <- tumour_events(20, mu = mu, covariance = Sigma, threshold = 1.2)
#'
#' les <- lesion_events(
#'   n_patients   = 20,
#'   log_tumour   = tum$log_ratios_vs_baseline,
#'   alpha        = -2.5,
#'   beta         = 0,
#'   gamma        = 0.2,
#'   treatment_arm = rep(0, 20)
#' )
#' colMeans(les$events == 1, na.rm = TRUE)
#'
#' @seealso \code{\link{tumour_events}}, \code{\link{simulate_trial}}
#'
#' @importFrom stats plogis pnorm rbinom runif
#' @export
lesion_events <- function(n_patients, log_tumour, alpha, beta, gamma,
                          treatment_arm,
                          baseline = runif(n_patients, 0, 1),
                          link = c("logit", "probit"),
                          nonlinear = FALSE) {
  link     <- match.arg(link)
  T_visits <- ncol(log_tumour)

  stopifnot(
    nrow(log_tumour) == n_patients,
    length(alpha) %in% c(1L, T_visits),
    length(beta)  %in% c(1L, T_visits),
    length(gamma) %in% c(1L, T_visits),
    length(treatment_arm) %in% c(1L, n_patients),
    length(baseline) == n_patients
  )

  alpha         <- rep_len(alpha, T_visits)
  beta          <- rep_len(beta,  T_visits)
  gamma         <- rep_len(gamma, T_visits)
  treatment_arm <- rep_len(treatment_arm, n_patients)

  actual_tumour <- baseline * exp(log_tumour)
  all_tumour    <- cbind(baseline, actual_tumour)
  tumour_lag    <- all_tumour[, seq_len(T_visits), drop = FALSE]
  tumour_term   <- if (nonlinear) tumour_lag^2 else tumour_lag

  linpred <-
    matrix(alpha, n_patients, T_visits, byrow = TRUE) +
    matrix(beta,  n_patients, T_visits, byrow = TRUE) * treatment_arm +
    matrix(gamma, n_patients, T_visits, byrow = TRUE) * tumour_term

  lesion_prob <- switch(link,
                        logit  = plogis(linpred),
                        probit = pnorm(linpred))

  # rbinom unfolds the prob matrix in column-major order; reshape matches
  lesion_indicator <- matrix(
    rbinom(length(lesion_prob), size = 1, prob = lesion_prob),
    nrow = n_patients,
    ncol = T_visits
  )

  # Censor visits strictly after the first event
  cum_events   <- t(apply(lesion_indicator == 1, 1, cumsum))
  prior_events <- cum_events - (lesion_indicator == 1)
  lesion_indicator[prior_events > 0] <- NA

  list(
    all_tumour_sizes   = all_tumour,
    lesion_probability = lesion_prob,
    events             = lesion_indicator
  )
}


#' Combine tumour and lesion events into a composite time-to-event outcome
#'
#' Combines two event-indicator matrices into a single time-to-first-event
#' outcome per patient, defined as the earliest visit at which either a
#' tumour progression or a new lesion event occurs. This is the standard
#' RECIST composite progression-free survival (PFS) endpoint. Patients with
#' no event by the final visit are administratively censored at the last
#' visit.
#'
#' @param lesion Numeric \code{n_patients} x \code{T} matrix of lesion event
#'   indicators (0/1), typically the \code{events} component returned by
#'   \code{\link{lesion_events}}.
#' @param tumour Numeric \code{n_patients} x \code{T} matrix of tumour event
#'   indicators, with the same dimensions as \code{lesion}, typically the
#'   \code{events} component returned by \code{\link{tumour_events}}.
#'
#' @return A \code{data.frame} with \code{n_patients} rows and two columns:
#' \describe{
#'   \item{\code{time}}{Integer visit index of the first event of either
#'     type, or \code{ncol(tumour)} if no event occurs.}
#'   \item{\code{status}}{Integer (0/1) event indicator; 0 if censored.}
#' }
#'
#' @details \code{NA} entries (post-event censoring rows produced by
#' \code{\link{tumour_events}} and \code{\link{lesion_events}}) are skipped
#' when finding the first event in each row.
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
#' head(pfs)
#' table(pfs$status)
#'
#' @seealso \code{\link{tumour_events}}, \code{\link{lesion_events}},
#'   \code{\link{first_event_times}}
#'
#' @export
combined_event <- function(tumour, lesion) {
  stopifnot(identical(dim(tumour), dim(lesion)))
  n_visits <- ncol(tumour)

  first_visit <- function(row) {
    idx <- which(row == 1)
    if (length(idx) == 0L) Inf else idx[1L]
  }

  t_time     <- apply(tumour, 1, first_visit)
  l_time     <- apply(lesion, 1, first_visit)
  first_time <- pmin(t_time, l_time)

  status <- as.integer(is.finite(first_time))
  time   <- as.integer(ifelse(is.finite(first_time), first_time, n_visits))

  data.frame(time = time, status = status)
}

#' Extract time-to-first-event from an event-indicator matrix
#'
#' Computes a per-patient time-to-event outcome from a single matrix of event
#' indicators by visit, returning the visit index of the first event and a
#' status flag indicating whether any event occurred. Patients with no event
#' are administratively censored at the final visit.
#'
#' @param events Numeric \code{n_patients} x \code{T} matrix of event
#'   indicators (0/1), typically with \code{NA} for visits after the first
#'   event.
#'
#' @return A \code{data.frame} with columns \code{time} (visit index of first
#'   event or \code{ncol(events)} if censored) and \code{status} (0/1).
#'
find_event_time <- function(events) {

  n_visits <- ncol(events)

  time <- apply(events, 1, function(row) {
    idx <- which(row == 1)
    if (length(idx) == 0) n_visits else idx[1]
  })

  status <- as.integer(apply(events, 1, function(row) any(row == 1, na.rm = TRUE)))

  data.frame(time = time, status = status)
}


#' Convert paired tumour and lesion event matrices to a composite (time, status) table
#'
#' Internal alias for \code{\link{combined_event}}; retained for backwards
#' compatibility within the package's internals.
#'
#' @param tumour,lesion Event-indicator matrices of identical dimension.
#' @return A data frame as in \code{\link{combined_event}}.
#' @keywords internal
first_event_times <- function(tumour, lesion) {
  combined_event(tumour, lesion)
}


#' Convert a single event-indicator matrix to a (time, status) table
#'
#' Loop-based per-endpoint analogue of \code{\link{find_event_time}}: takes
#' one event-indicator matrix (tumour or lesion) and returns first-event
#' times, with administrative censoring at the final visit for patients
#' with no event.
#'
#' @param events_mat Numeric \code{n_patients} x \code{T} matrix of event
#'   indicators (0/1).
#'
#' @return A \code{data.frame} with columns \code{time} and \code{status}.
#'
#'
separate_event_times <- function(events_mat) {
  T_visits <- ncol(events_mat)
  n        <- nrow(events_mat)
  time   <- integer(n)
  status <- integer(n)
  for (i in seq_len(n)) {
    hit <- which(events_mat[i, ] == 1)
    if (length(hit) == 0) { time[i] <- T_visits; status[i] <- 0 }
    else                  { time[i] <- min(hit); status[i] <- 1 }
  }
  data.frame(time = time, status = status)
}

#' Apply administrative censoring at a cutoff visit
#'
#' Truncates a composite (\code{time}, \code{status}) data frame at a
#' user-specified administrative cutoff: events occurring strictly after
#' \code{cens_at} are marked as censored at \code{cens_at}. Useful for
#' simulating short-follow-up scenarios or sensitivity analyses where the
#' study is treated as if it had ended earlier than the simulated horizon.
#'
#' @param data A \code{data.frame} with columns \code{time} and
#'   \code{status}, e.g. as returned by \code{\link{combined_event}} or
#'   \code{\link{first_event_times}}.
#' @param cens_at Numeric cutoff. Observations with \code{time > cens_at}
#'   have their \code{status} set to 0 and \code{time} replaced by
#'   \code{cens_at}. If \code{NULL} or \code{Inf}, \code{data} is returned
#'   unchanged.
#'
#' @return A data frame of the same shape as \code{data}, with \code{time}
#'   capped at \code{cens_at} and \code{status} set to 0 for any event
#'   after the cutoff.
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
#'
#' # Censor at visit 3
#' pfs_short <- apply_admin_censoring(pfs, cens_at = 3)
#' table(pfs$status, pfs_short$status)
#'
#' @seealso \code{\link{combined_event}}
#'
#' @export
apply_admin_censoring <- function(data, cens_at) {
  if (is.null(cens_at) || is.infinite(cens_at)) return(data)

  stopifnot(
    is.data.frame(data),
    all(c("time", "status") %in% names(data)),
    is.numeric(cens_at),
    length(cens_at) == 1L
  )

  if (cens_at < 1) {
    warning("cens_at = ", cens_at,
            " censors all observations; this is rarely intended.")
  }

  censored <- data$time > cens_at
  data$status[censored] <- 0L
  data$time[censored]   <- cens_at

  data
}

#' Simulate a single two-arm trial
#'
#' Convenience wrapper that simulates control and treatment arms under the
#' same data-generating process and returns the per-arm tumour and lesion
#' simulation outputs in a single nested list. Tumour log-ratios are drawn
#' with arm-specific means \code{mu_ctrl} and \code{mu_trt} but a shared
#' covariance matrix; lesion events are then generated conditional on the
#' tumour trajectory in each arm. Setting \code{gamma_lesion_ctrl} and
#' \code{gamma_lesion_trt} separately allows the tumour-to-lesion coupling
#' to differ between arms, which is the structure used to study power under
#' arm-asymmetric dependence in Regan (2026, Section 5.2.2).
#'
#' @param n_per_arm Integer. Number of patients in each arm.
#' @param mu_ctrl,mu_trt Numeric vectors of length \code{T} giving the
#'   tumour-size log-ratio means for control and treatment arms.
#' @param cov_mat Numeric \code{T} x \code{T} covariance matrix for the
#'   tumour trajectory, shared by both arms. For arm-specific covariances,
#'   call \code{\link{tumour_events}} and \code{\link{lesion_events}}
#'   directly.
#' @param threshold Numeric (> 1). Multiplicative tumour-progression
#'   threshold. Defaults to 1.2 (RECIST v1.1).
#' @param alpha_lesion,beta_lesion,gamma_lesion Lesion-model coefficients
#'   passed to \code{\link{lesion_events}}. \code{beta_lesion} is the
#'   treatment-arm coefficient (control arm coded as 0, treatment arm as 1).
#' @param gamma_lesion_ctrl,gamma_lesion_trt Optional per-arm overrides for
#'   the tumour-to-lesion coefficient. When non-\code{NULL}, these take
#'   precedence over \code{gamma_lesion} for the respective arm.
#' @param tumour_dist,tumour_df Tumour-distribution arguments passed to
#'   \code{\link{tumour_events}}.
#' @param lesion_link,lesion_nonlinear Lesion-model arguments passed to
#'   \code{\link{lesion_events}}.
#'
#' @return A list with elements \code{ctrl} and \code{trt}, each containing
#'   sub-lists \code{tumour} (from \code{\link{tumour_events}}) and
#'   \code{lesion} (from \code{\link{lesion_events}}) for that arm.
#'
#' @examples
#' set.seed(1)
#' mu_ctrl <- c(-0.10, -0.30, -0.46, -0.50, -0.55)
#' mu_trt  <- c(-0.20, -0.40, -0.56, -0.60, -0.65)
#' Sigma <- matrix(0.05, 5, 5)
#' diag(Sigma) <- c(0.05, 0.10, 0.14, 0.16, 0.18)
#'
#' # Scenario 1 from Regan (2026): null treatment effect, no dependence
#' trial <- simulate_trial(
#'   n_per_arm    = 50,
#'   mu_ctrl      = mu_ctrl,
#'   mu_trt       = mu_trt,
#'   cov_mat      = Sigma,
#'   beta_lesion  = 0,
#'   gamma_lesion = 0
#' )
#' str(trial, max.level = 2)
#'
#' \donttest{
#' # Scenario 6: asymmetric dependence, small treatment effect
#' trial2 <- simulate_trial(
#'   n_per_arm         = 200,
#'   mu_ctrl           = mu_ctrl,
#'   mu_trt            = mu_trt,
#'   cov_mat           = Sigma,
#'   beta_lesion       = -0.5,
#'   gamma_lesion_ctrl = 1.0,
#'   gamma_lesion_trt  = 0.2
#' )
#' }
#'
#' @seealso \code{\link{tumour_events}}, \code{\link{lesion_events}},
#'   \code{\link{pfs_copula}}
#' @references
#' Regan, C. (2026). \emph{Copula Modelling for Mixed Continuous and Binary
#' Variables in Survival Analysis: Applications to Colorectal Cancer}.
#' MMath dissertation.
#'
#' @export
simulate_trial <- function(n_per_arm, mu_ctrl, mu_trt, cov_mat,
                           threshold = 1.2,
                           alpha_lesion = -2.5, beta_lesion = 0,
                           gamma_lesion = 0.2,
                           gamma_lesion_ctrl = NULL,
                           gamma_lesion_trt  = NULL,
                           tumour_dist = c("mvnorm", "mvt"),
                           tumour_df   = 4,
                           lesion_link = c("logit", "probit"),
                           lesion_nonlinear = FALSE) {
  tumour_dist <- match.arg(tumour_dist)
  lesion_link <- match.arg(lesion_link)

  stopifnot(
    length(mu_ctrl) == length(mu_trt),
    length(mu_ctrl) == nrow(cov_mat),
    n_per_arm >= 1
  )

  if (is.null(gamma_lesion_ctrl)) gamma_lesion_ctrl <- gamma_lesion
  if (is.null(gamma_lesion_trt))  gamma_lesion_trt  <- gamma_lesion

  t_ctrl <- tumour_events(n_per_arm, mu_ctrl, cov_mat, threshold,
                          dist = tumour_dist, df = tumour_df)
  l_ctrl <- lesion_events(n_per_arm, t_ctrl$log_ratios_vs_baseline,
                          alpha_lesion, beta_lesion, gamma_lesion_ctrl,
                          treatment_arm = 0,
                          link = lesion_link, nonlinear = lesion_nonlinear)

  t_trt  <- tumour_events(n_per_arm, mu_trt, cov_mat, threshold,
                          dist = tumour_dist, df = tumour_df)
  l_trt  <- lesion_events(n_per_arm, t_trt$log_ratios_vs_baseline,
                          alpha_lesion, beta_lesion, gamma_lesion_trt,
                          treatment_arm = 1,
                          link = lesion_link, nonlinear = lesion_nonlinear)

  list(
    ctrl = list(tumour = t_ctrl, lesion = l_ctrl),
    trt  = list(tumour = t_trt,  lesion = l_trt)
  )
}
