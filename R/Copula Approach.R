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

#' Estimate marginal survival functions for copula fitting
#'
#' Fits non-parametric Kaplan-Meier survival curves to lesion and tumour
#' progression event data and evaluates the survival probability at each
#' follow-up timepoint. The resulting marginal estimates are intended for
#' use as inputs to a copula model.
#'
#' @param lesion_events A data frame with columns \code{time} and \code{status}
#'   as returned by \code{lesion_event()}.
#' @param tumour_events A data frame with columns \code{time} and \code{status}
#'   as returned by \code{tumour_event()}.
#' @param n_times A positive integer giving the number of follow-up timepoints
#'   T at which to evaluate the survival functions.
#'
#' @return A data frame with \code{n_times} rows and two columns:
#'   \describe{
#'     \item{LESION}{Kaplan-Meier survival probability for lesion progression
#'       at each timepoint t = 1, ..., T.}
#'     \item{TUMOUR}{Kaplan-Meier survival probability for tumour progression
#'       at each timepoint t = 1, ..., T.}
#'   }
copula_margin_estimation <- function(lesion_events, tumour_events, n_times) {
  fit_lesion <- survfit(Surv(time, status) ~ 1, data = lesion_events)
  fit_tumour <- survfit(Surv(time, status) ~ 1, data = tumour_events)

  times <- seq_len(n_times)

  data.frame(
    LESION = summary(fit_lesion, times = times, extend = TRUE)$surv,
    TUMOUR = summary(fit_tumour, times = times, extend = TRUE)$surv
  )
}




