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
#' Simulates tumour progression events for a cohort of patients based on
#' multivariate-normal (or multivariate-t) log-ratios of tumour size versus
#' baseline. An event is declared at the first visit where the tumour size
#' exceeds the running nadir (minimum prior tumour size) by a multiplicative
#' factor greater than \code{threshold}. Visits after the first event are
#' censored (set to \code{NA}).
#'
#' @param n_patients Integer. Number of patients to simulate.
#' @param mean Numeric vector of length \code{T} giving the mean of the
#'   multivariate distribution for \eqn{Y_t = \log(z_t / z_0)}, the log-ratio
#'   of tumour size at visit \eqn{t} relative to baseline, for visits
#'   \eqn{t = 1, \ldots, T}.
#' @param covariance Numeric \code{T} x \code{T} covariance (or scale) matrix.
#' @param threshold Numeric (> 1). Multiplicative threshold on tumour size
#'   relative to the running nadir. An event occurs at visit \eqn{t} if
#'   \eqn{z_t / \mathrm{nadir}_{<t} > \mathrm{threshold}}.
#' @param dist Character. Either \code{"mvnorm"} (default) for a multivariate
#'   normal log-ratio distribution, or \code{"mvt"} for a shifted
#'   multivariate-t distribution (useful for studying robustness to
#'   heavy-tailed tumour trajectories).
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
#' @details The baseline log-ratio is fixed at zero by construction
#'   (\eqn{Y_0 = 0}), and the running nadir at visit \eqn{t} is taken over
#'   visits \eqn{0, 1, \ldots, t-1}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @export
tumour_events <- function(n_patients, mean, covariance, threshold,
                          dist = "mvnorm", df = 4) {

  if (dist == "mvt") {
    Y <- mvtnorm::rmvt(n = n_patients, sigma = covariance, df = df,
                       delta = mean, type = "shifted")
  } else {
    Y <- MASS::mvrnorm(n = n_patients, mu = mean, Sigma = covariance)
  }

  # log(z_t / z_0) for t = 0, 1, ..., T. By definition Y_00 = 0.
  log_ratio_vs_baseline <- cbind(0, Y)

  # Running minimum of log(z_s / z_0) over s = 0, ..., t-1
  running_min_prior <- t(apply(log_ratio_vs_baseline, 1, function(r) {
    cummin(r)[-length(r)]
  }))

  log_ratio_vs_nadir <- Y - running_min_prior
  event_indicators   <- ifelse(log_ratio_vs_nadir > log(threshold), 1, 0)

  # Censor after first event
  cum_ones <- t(apply(event_indicators == 1, 1, cumsum))
  mask <- cbind(0, cum_ones[, -ncol(cum_ones), drop = FALSE]) > 0
  event_indicators[mask] <- NA

  list(
    log_ratios_vs_baseline = Y,
    log_ratios_vs_nadir    = log_ratio_vs_nadir,
    events                 = event_indicators
  )
}


#' Simulate new-lesion events conditional on tumour trajectory
#'
#' Simulates the occurrence of new lesions over time as Bernoulli events
#' whose probability depends on treatment arm and the lagged (previous-visit)
#' tumour size through a logistic (or probit) linear predictor. Visits after
#' the first lesion event are censored (set to \code{NA}).
#'
#' @param n_patients Integer. Number of patients to simulate.
#' @param log_tumour Numeric \code{n_patients} x \code{T} matrix of log-ratios
#'   of tumour size versus baseline, typically the
#'   \code{log_ratios_vs_baseline} element of \code{\link{tumour_events}}.
#' @param alpha Numeric scalar or vector of length \code{T}. Intercept of the
#'   linear predictor at each visit. Recycled to length \code{T}.
#' @param beta Numeric scalar or vector of length \code{T}. Treatment-arm
#'   coefficient at each visit. Recycled to length \code{T}.
#' @param gamma Numeric scalar or vector of length \code{T}. Coefficient on
#'   the lagged tumour size at each visit. Recycled to length \code{T}.
#' @param treatment_arm Numeric scalar or vector of length \code{n_patients}.
#'   Treatment-arm indicator (typically 0/1). Recycled to length
#'   \code{n_patients}.
#' @param link Character. Either \code{"logit"} (default) or \code{"probit"}.
#' @param nonlinear Logical. If \code{TRUE}, the tumour-size effect enters
#'   the linear predictor as \code{tumour_lag^2} instead of \code{tumour_lag}
#'   (useful for misspecification studies). Defaults to \code{FALSE}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{all_tumour_sizes}}{An \code{n_patients} x \code{(T+1)} matrix
#'     of absolute tumour sizes, with the (Uniform(0,1)) baseline size in the
#'     first column.}
#'   \item{\code{lesion_probability}}{An \code{n_patients} x \code{T} matrix
#'     of per-visit lesion probabilities from the linear predictor.}
#'   \item{\code{events}}{An \code{n_patients} x \code{T} matrix of lesion
#'     event indicators (0/1), with \code{NA} for visits strictly after the
#'     first event.}
#' }
#'
#' @details The linear predictor at visit \eqn{t} for patient \eqn{i} is
#'   \eqn{\eta_{it} = \alpha_t + \beta_t \cdot \mathrm{arm}_i + \gamma_t \cdot z_{i,t-1}}
#'   (or \eqn{\gamma_t \cdot z_{i,t-1}^2} when \code{nonlinear = TRUE}). The
#'   lesion probability is then \eqn{p_{it} = g^{-1}(\eta_{it})} with
#'   \eqn{g^{-1}} the logistic or normal CDF.
#'
#' @importFrom stats plogis pnorm rbinom runif
#' @export
lesion_events <- function(n_patients, log_tumour, alpha, beta, gamma,
                          treatment_arm, link = "logit", nonlinear = FALSE) {

  T_visits      <- ncol(log_tumour)
  alpha         <- rep_len(alpha, T_visits)
  beta          <- rep_len(beta,  T_visits)
  gamma         <- rep_len(gamma, T_visits)
  treatment_arm <- rep_len(treatment_arm, n_patients)

  baseline      <- runif(n_patients, min = 0, max = 1)
  actual_tumour <- baseline * exp(log_tumour)
  all_tumour    <- cbind(baseline, actual_tumour)
  tumour_lag    <- all_tumour[, 1:T_visits, drop = FALSE]
  tumour_term   <- if (nonlinear) tumour_lag^2 else tumour_lag

  linpred <-
    matrix(alpha, n_patients, T_visits, byrow = TRUE) +
    matrix(beta,  n_patients, T_visits, byrow = TRUE) * treatment_arm +
    matrix(gamma, n_patients, T_visits, byrow = TRUE) * tumour_term

  lesion_probability_it <- if (link == "probit") pnorm(linpred) else plogis(linpred)
  lesion_probability_it <- as.matrix(lesion_probability_it)

  lesion_indicator <- matrix(
    rbinom(length(lesion_probability_it), size = 1, prob = lesion_probability_it),
    nrow = nrow(lesion_probability_it),
    ncol = ncol(lesion_probability_it)
  )

  cum_ones <- t(apply(lesion_indicator == 1, 1, cumsum))
  mask <- cbind(0, cum_ones[, -ncol(cum_ones), drop = FALSE]) > 0
  lesion_indicator[mask] <- NA

  list(all_tumour_sizes   = all_tumour,
       lesion_probability = lesion_probability_it,
       events             = lesion_indicator)
}


#' Combine tumour and lesion events into a composite time-to-event outcome
#'
#' Combines two event-indicator matrices (one for tumour progression and one
#' for new lesions) into a single time-to-first-event outcome per patient,
#' defined as the earliest visit at which either event occurs. Patients with
#' no event by the final visit are administratively censored.
#'
#' @param lesion_events Numeric \code{n_patients} x \code{T} matrix of lesion
#'   event indicators (0/1), e.g. the \code{events} component of
#'   \code{\link{lesion_events}}. \code{NA} values are treated as non-events.
#' @param tumour_events Numeric \code{n_patients} x \code{T} matrix of tumour
#'   event indicators (0/1), with the same dimensions as \code{lesion_events}.
#'
#' @return A \code{data.frame} with \code{n_patients} rows and two columns:
#' \describe{
#'   \item{\code{time}}{Integer visit index of the first event of either
#'     type, or \code{ncol(tumour_events)} if no event occurs.}
#'   \item{\code{status}}{Integer (0/1) event indicator; 0 if censored.}
#' }
#'
#' @export
combined_event <- function(lesion_events, tumour_events) {
  stopifnot(identical(dim(lesion_events), dim(tumour_events)))

  n_visits <- ncol(tumour_events)

  first_visit <- function(row) {
    idx <- which(row == 1)
    if (length(idx) == 0) Inf else idx[1]
  }

  t_time <- apply(tumour_events, 1, first_visit)
  l_time <- apply(lesion_events, 1, first_visit)

  first_time <- pmin(t_time, l_time)
  status     <- as.integer(is.finite(first_time))
  time       <- ifelse(is.finite(first_time), first_time, n_visits)

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
#' @export
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
#' Loop-based implementation of \code{\link{combined_event}} that does not
#' assume both matrices have identical dimensions. Used internally by the
#' two-arm power simulators where tumour and lesion events arrive as
#' separate matrices and the composite endpoint is needed alongside the
#' per-endpoint endpoints.
#'
#' @param tumour_events_mat Numeric \code{n_patients} x \code{T} matrix of
#'   tumour event indicators (0/1).
#' @param lesion_events_mat Numeric \code{n_patients} x \code{T} matrix of
#'   lesion event indicators (0/1).
#'
#' @return A \code{data.frame} with columns \code{time} and \code{status}.
#'
#' @export
first_event_times <- function(tumour_events_mat, lesion_events_mat) {
  T_visits <- ncol(tumour_events_mat)
  n        <- nrow(tumour_events_mat)
  time   <- integer(n)
  status <- integer(n)
  for (i in seq_len(n)) {
    t_hit <- which(tumour_events_mat[i, ] == 1)
    l_hit <- which(lesion_events_mat[i, ] == 1)
    hits <- c(t_hit, l_hit)
    if (length(hits) == 0) { time[i] <- T_visits; status[i] <- 0 }
    else                   { time[i] <- min(hits); status[i] <- 1 }
  }
  data.frame(time = time, status = status)
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
#' @export
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
#' simulating short follow-up scenarios.
#'
#' @param ce A \code{data.frame} with columns \code{time} and \code{status}
#'   (e.g. as returned by \code{\link{combined_event}} or
#'   \code{\link{find_event_time}}).
#' @param cens_at Numeric cutoff visit. If \code{NULL} or \code{Inf}, the
#'   input is returned unchanged.
#'
#' @return The input data frame with \code{time} capped at \code{cens_at}
#'   and \code{status} set to 0 for any event after the cutoff.
#'
#' @export
apply_admin_censoring <- function(ce, cens_at) {
  if (is.null(cens_at) || is.infinite(cens_at)) return(ce)
  censored <- ce$time > cens_at
  ce$status[censored] <- 0L
  ce$time[censored]   <- cens_at
  ce
}


#' Simulate a single two-arm trial
#'
#' Convenience wrapper that simulates control and treatment arms under the
#' same data-generating process (with optional arm-specific tumour means
#' and dependence parameters) and returns the per-arm tumour and lesion
#' simulation outputs in a single nested list.
#'
#' @param n_per_arm Integer. Number of patients in each arm.
#' @param mean_ctrl,mean_trt Numeric vectors of length \code{T} giving the
#'   tumour-size log-ratio means for control and treatment arms.
#' @param cov_mat Numeric \code{T} x \code{T} covariance matrix for the
#'   tumour trajectory, shared by both arms.
#' @param threshold Numeric (> 1). Multiplicative tumour-progression
#'   threshold passed to \code{\link{tumour_events}}. Defaults to 1.2.
#' @param alpha_lesion,beta_lesion,gamma_lesion Lesion-model coefficients
#'   passed to \code{\link{lesion_events}}.
#' @param gamma_lesion_ctrl,gamma_lesion_trt Optional per-arm overrides for
#'   the tumour-to-lesion coefficient. If \code{NULL}, both arms use
#'   \code{gamma_lesion}.
#' @param tumour_dist,tumour_df Tumour-distribution arguments passed to
#'   \code{\link{tumour_events}}.
#' @param lesion_link,lesion_nonlinear Lesion-model arguments passed to
#'   \code{\link{lesion_events}}.
#'
#' @return A list with elements \code{ctrl} and \code{trt}, each containing
#'   sub-lists \code{tumour} (from \code{\link{tumour_events}}) and
#'   \code{lesion} (from \code{\link{lesion_events}}) for that arm.
#'
#' @export
simulate_trial <- function(n_per_arm, mean_ctrl, mean_trt, cov_mat,
                           threshold = 1.2,
                           alpha_lesion = -2.5, beta_lesion = 0,
                           gamma_lesion = 0.2,
                           gamma_lesion_ctrl = NULL,
                           gamma_lesion_trt  = NULL,
                           tumour_dist = "mvnorm", tumour_df = 4,
                           lesion_link = "logit",
                           lesion_nonlinear = FALSE) {
  if (is.null(gamma_lesion_ctrl)) gamma_lesion_ctrl <- gamma_lesion
  if (is.null(gamma_lesion_trt))  gamma_lesion_trt  <- gamma_lesion

  t_ctrl <- tumour_events(n_per_arm, mean_ctrl, cov_mat, threshold,
                          dist = tumour_dist, df = tumour_df)
  l_ctrl <- lesion_events(n_per_arm, t_ctrl$log_ratios_vs_baseline,
                          alpha_lesion, beta_lesion, gamma_lesion_ctrl,
                          treatment_arm = 0,
                          link = lesion_link, nonlinear = lesion_nonlinear)
  t_trt  <- tumour_events(n_per_arm, mean_trt, cov_mat, threshold,
                          dist = tumour_dist, df = tumour_df)
  l_trt  <- lesion_events(n_per_arm, t_trt$log_ratios_vs_baseline,
                          alpha_lesion, beta_lesion, gamma_lesion_trt,
                          treatment_arm = 1,
                          link = lesion_link, nonlinear = lesion_nonlinear)

  list(ctrl = list(tumour = t_ctrl, lesion = l_ctrl),
       trt  = list(tumour = t_trt,  lesion = l_trt))
}
