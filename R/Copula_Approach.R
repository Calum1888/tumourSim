library(survival)

#' Derive lesion-based event times and status indicators
#'
#' Computes the first time at which a binary lesion indicator equals 1 for each
#' subject, returning both the event time and a censoring indicator. If a subject
#' never experiences a lesion event, the time is set to the final observation
#' time and the status is 0 (censored).
#'
#' @param lesion_data A matrix or data frame where each row corresponds to a
#'   subject and each column corresponds to a follow-up time. Entries must be
#'   binary (0/1), with 1 indicating a lesion event at that time.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{time}{The index of the first time point at which a lesion event
#'       occurs, or the final time point if no event occurs.}
#'     \item{status}{Event indicator: 1 if an event occurs, 0 if censored.}
#'   }
#'
#' @examples
#' lesion_data <- matrix(c(0,0,1,
#'                         0,0,0), nrow = 2, byrow = TRUE)
#' lesion_event(lesion_data)
#' @export
lesion_event <- function(lesion_data) {
  n_times <- ncol(lesion_data)

  event_matrix <- lesion_data == 1

  no_event       <- rowSums(event_matrix) == 0
  idx            <- max.col(event_matrix, ties.method = "first")
  time           <- idx
  time[no_event] <- n_times
  status         <- as.integer(!no_event)

  return(data.frame(time = time, status = status))
}

#' Derive tumour progression event times based on relative growth threshold
#'
#' Identifies the earliest time at which tumour size exceeds a specified
#' multiplicative increase over the subject-specific historical minimum. This
#' implements a rolling-minimum progression rule commonly used in tumour
#' progression modelling.
#'
#' @param tumour_size_data A matrix or data frame where each row corresponds to a
#'   subject and each column corresponds to a tumour size measurement at a given
#'   time. The first column is treated as baseline and excluded from event
#'   detection.
#' @param threshold A numeric multiplier defining progression. An event is
#'   recorded when \eqn{tumour\_size > threshold × rolling\_minimum}. Defaults to
#'   1.2.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{time}{The index (relative to follow-up columns) of the first time
#'       point at which the threshold is exceeded, or the final time point if no
#'       event occurs.}
#'     \item{status}{Event indicator: 1 if progression occurs, 0 if censored.}
#'   }
tumour_event <- function(tumour_size_data, threshold = 1.2) {
  n_times <- ncol(tumour_size_data) - 1

  tumour      <- tumour_size_data[, -1, drop = FALSE]
  prior_cols  <- tumour_size_data[, -ncol(tumour_size_data), drop = FALSE]
  rolling_min <- t(apply(prior_cols, 1, cummin))

  event_matrix <- tumour > threshold * rolling_min

  no_event       <- rowSums(event_matrix) == 0
  idx            <- max.col(event_matrix, ties.method = "first")
  time           <- idx
  time[no_event] <- n_times
  status         <- as.integer(!no_event)

  return(data.frame(time = time, status = status))
}

#' Compute pseudo-observations for copula estimation
#'
#' Converts subject-level event times into pseudo-observations using
#' Kaplan–Meier marginal CDFs.
#'
#' @param events A data frame with columns time and status.
#'
#' @return A numeric vector of pseudo-observations in (0,1).
pseudo_obs <- function(events) {

  fit <- survfit(Surv(time, status) ~ 1, data = events)

  # Evaluate survival at each subject's event time
  S_at_time <- summary(fit, times = events$time, extend = TRUE)$surv

  # Convert to CDF
  U <- 1 - S_at_time

  # Avoid boundary issues
  eps <- 1e-6
  U <- pmin(pmax(U, eps), 1 - eps)

  return(U)
}

#' Estimate marginal survival curves for lesion and tumour progression
#'
#' Computes non‑parametric Kaplan–Meier survival estimates for lesion
#' progression and tumour progression at each discrete follow‑up time
#' point. These marginal survival functions are required inputs for the
#' copula‑based reconstruction of the progression‑free survival (PFS)
#' curve.
#'
#' @param lesion_events A data frame with columns \code{time} and
#'   \code{status}, typically produced by \code{lesion_event()}.
#' @param tumour_events A data frame with columns \code{time} and
#'   \code{status}, typically produced by \code{tumour_event()}.
#' @param n_times Integer. Number of follow‑up time points at which to
#'   evaluate the marginal survival curves.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{LESION}{Kaplan–Meier survival probability for lesion
#'       progression at times \code{1:n_times}.}
#'     \item{TUMOUR}{Kaplan–Meier survival probability for tumour
#'       progression at times \code{1:n_times}.}
#'   }
#'
#' @details
#' The function fits two independent Kaplan–Meier curves—one for lesion
#' progression and one for tumour progression—and evaluates each curve at
#' the discrete time points \code{1, ..., n_times}. These marginal
#' survival estimates are used in the survival‑copula identity for PFS.
#'
#' @examples
#' \dontrun{
#' lesion_events <- lesion_event(matrix(c(0,0,1, 0,0,0), 2, byrow=TRUE))
#' tumour_events <- tumour_event(matrix(c(1,1.1,1.3, 1,1,1), 2, byrow=TRUE))
#' copula_margin_estimation(lesion_events, tumour_events, n_times = 3)
#' }
#'
#' @export
copula_margin_estimation <- function(lesion_events, tumour_events, n_times) {

  fit_lesion <- survfit(Surv(time, status) ~ 1, data = lesion_events)
  fit_tumour <- survfit(Surv(time, status) ~ 1, data = tumour_events)

  times <- seq_len(n_times)

  data.frame(
    LESION = summary(fit_lesion, times = times, extend = TRUE)$surv,
    TUMOUR = summary(fit_tumour, times = times, extend = TRUE)$surv
  )
}

#' Estimate progression‑free survival using a copula model
#'
#' Computes a copula‑based estimate of the progression‑free survival (PFS)
#' curve by combining marginal Kaplan–Meier survival functions for lesion
#' progression and tumour progression with a dependence structure
#' estimated from pseudo‑observations.
#'
#' @param lesion_events A data frame with columns \code{time} and
#'   \code{status}, typically produced by \code{lesion_event()}.
#' @param tumour_events A data frame with columns \code{time} and
#'   \code{status}, typically produced by \code{tumour_event()}.
#' @param n_times Integer. Number of follow‑up time points at which to
#'   evaluate the PFS curve.
#' @param copula_family Integer or vector of integers specifying the
#'   copula family (or families) to be considered by
#'   \code{VineCopula::BiCopSelect()}.
#'
#' @return A numeric vector of length \code{n_times} giving the estimated
#'   progression‑free survival probability at each time point.
#'
#' @details
#' The method proceeds in three steps:
#' \enumerate{
#'   \item Pseudo‑observations are computed from the marginal event times
#'         using Kaplan–Meier CDFs.
#'   \item A bivariate copula is fitted to the pseudo‑observations to
#'         estimate the dependence parameter.
#'   \item The PFS curve is reconstructed using the survival‑copula
#'         identity:
#'         \deqn{
#'           S_{\mathrm{PFS}}(t)
#'           = S_D(t) + S_Y(t) - 1 + C(1 - S_D(t), 1 - S_Y(t); \theta),
#'         }
#'         where \eqn{S_D} and \eqn{S_Y} are the marginal survival
#'         functions and \eqn{C} is the fitted copula.
#' }
#'
#' @examples
#' \dontrun{
#' lesion_events <- lesion_event(matrix(c(0,0,1, 0,0,0), 2, byrow=TRUE))
#' tumour_events <- tumour_event(matrix(c(1,1.1,1.3, 1,1,1), 2, byrow=TRUE))
#' copula_pfs(lesion_events, tumour_events, n_times = 3, copula_family = 1)
#' }
#'
#' @export
copula_pfs <- function(lesion_events, tumour_events, n_times, copula_family) {

  # compute pseudo-observations ----
  U_D <- pseudo_obs(lesion_events)
  U_Y <- pseudo_obs(tumour_events)

  # fit copula to subject-level pseudo-observations ----
  copula_fit <- BiCopSelect(U_D, U_Y, familyset = copula_family)

  # estimate marginal survival curves at each time ----
  margins <- copula_margin_estimation(lesion_events, tumour_events, n_times)
  S_D <- margins$LESION
  S_Y <- margins$TUMOUR

  # convert to CDFs for the PFS identity
  U_D_grid <- 1 - S_D
  U_Y_grid <- 1 - S_Y

  # evaluate copula C(u,v) at each time ----
  C_uv <- BiCopCDF(U_D_grid, U_Y_grid, obj = copula_fit)

  # apply survival-copula identity ----
  S_pfs <- S_D + S_Y - 1 + C_uv

  return(S_pfs)
}
