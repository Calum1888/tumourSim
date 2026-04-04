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
#' @importFrom survival survfit Surv
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
#' @importFrom VineCopula BiCopSelect BiCopCDF
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

#' Bootstrap progression-free survival estimates from a copula model
#'
#' Draws \code{B} bootstrap resamples of subjects, re-estimates the copula-based
#' PFS curve on each resample, and returns pointwise percentile confidence
#' intervals, coverage against a known true curve, and average CI width.
#'
#' @param lesion_data A matrix where each row is a subject and each column a
#'   follow-up time. Entries are binary (0/1); passed to \code{lesion_event()}.
#' @param tumour_size_data A matrix where each row is a subject and each column
#'   a tumour size measurement. The first column is baseline; passed to
#'   \code{tumour_event()}.
#' @param n_times Integer. Number of follow-up time points at which to evaluate
#'   the PFS curve.
#' @param copula_family Integer or vector of integers specifying the copula
#'   family (or families) passed to \code{VineCopula::BiCopSelect()}.
#'   1 = Gaussian
#'   2 = t-Student
#'   3 = Clayton
#'   4 = Gumbel
#'   5 = Frank
#' @param B Integer. Number of bootstrap resamples. Defaults to 500.
#' @param alpha_level Numeric. Significance level for the percentile CI. Defaults to
#'   0.05, giving a 95\% CI.
#' @param true_pfs A numeric vector of length \code{n_times} containing the
#'   known true PFS probabilities at each time point, used to compute pointwise
#'   coverage.
#' @param threshold Numeric. Multiplicative threshold for tumour progression,
#'   passed to \code{tumour_event()}. Defaults to 1.2.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility. Defaults
#'   to \code{NULL}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{boot_curves}{A \code{B x n_times} matrix of bootstrap PFS
#'       estimates, one row per resample.}
#'     \item{ci_lower}{Numeric vector of length \code{n_times}. Lower bound of
#'       the pointwise percentile CI at each time point.}
#'     \item{ci_upper}{Numeric vector of length \code{n_times}. Upper bound of
#'       the pointwise percentile CI at each time point.}
#'     \item{coverage}{Numeric vector of length \code{n_times}. Indicator (0/1)
#'       of whether the true PFS value at each time point falls within the CI.}
#'     \item{mean_coverage}{Scalar. Average pointwise coverage across all time
#'       points.}
#'     \item{mean_ci_width}{Scalar. Average pointwise CI width across all time
#'       points.}
#'   }
#'
#' @details
#' On each bootstrap iteration, \code{n} subjects are sampled with replacement
#' (where \code{n} is the number of rows in \code{lesion_data}). The resampled
#' lesion and tumour matrices are passed through \code{lesion_event()},
#' \code{tumour_event()}, and \code{copula_pfs()} to produce a bootstrap PFS
#' curve. Any iteration that throws an error (e.g. due to degenerate resamples)
#' is silently skipped. Coverage is computed pointwise: at each time \code{t},
#' the indicator is 1 if \code{true_pfs[t]} lies within the percentile interval
#' [\code{ci_lower[t]}, \code{ci_upper[t]}].
#'
#' @examples
#' \dontrun{
#' lesion_data      <- matrix(rbinom(100 * 10, 1, 0.2), nrow = 100)
#' tumour_size_data <- matrix(abs(rnorm(100 * 11, mean = 1)), nrow = 100)
#' true_pfs         <- exp(-0.1 * 1:10)
#'
#' result <- bootstrap_copula_pfs(
#'   lesion_data      = lesion_data,
#'   tumour_size_data = tumour_size_data,
#'   n_times          = 10,
#'   copula_family    = 1,
#'   B                = 500,
#'   true_pfs         = true_pfs,
#'   seed             = 42
#' )
#'
#' result$mean_coverage
#' result$mean_ci_width
#' }
#'
#' @export
bootstrap_copula_pfs <- function(lesion_data,
                                 tumour_size_data,
                                 n_times,
                                 copula_family,
                                 B         = 500,
                                 alpha_level     = 0.05,
                                 true_pfs,
                                 threshold = 1.2,
                                 seed      = NULL) {

  # input checks ----
  stopifnot(
    is.matrix(lesion_data) || is.data.frame(lesion_data),
    is.matrix(tumour_size_data) || is.data.frame(tumour_size_data),
    nrow(lesion_data) == nrow(tumour_size_data),
    length(true_pfs) == n_times,
    B > 0, alpha_level > 0, alpha_level < 1
  )

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(lesion_data)

  # run bootstrap resamples ----
  boot_curves <- matrix(NA_real_, nrow = B, ncol = n_times)

  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)

    pfs_b <- tryCatch(
      {
        lesion_events_b <- lesion_event(lesion_data[idx, , drop = FALSE])
        tumour_events_b <- tumour_event(
          tumour_size_data[idx, , drop = FALSE],
          threshold = threshold
        )
        copula_pfs(lesion_events_b, tumour_events_b, n_times, copula_family)
      },
      error = function(e) NULL
    )

    if (!is.null(pfs_b)) boot_curves[b, ] <- pfs_b
  }

  # drop any failed iterations before summarising ----
  boot_curves <- boot_curves[complete.cases(boot_curves), , drop = FALSE]

  n_valid <- nrow(boot_curves)
  if (n_valid == 0) stop("All bootstrap iterations failed. Check your inputs.")
  if (n_valid < B)  warning(sprintf("%d of %d bootstrap iterations failed and were dropped.", B - n_valid, B))

  # pointwise percentile CIs ----
  ci_lower <- apply(boot_curves, 2, quantile, probs = alpha_level / 2)
  ci_upper <- apply(boot_curves, 2, quantile, probs = 1 - alpha_level / 2)

  # pointwise coverage against known true curve ----
  coverage <- as.integer(true_pfs >= ci_lower & true_pfs <= ci_upper)

  # summary scalars ----
  mean_coverage  <- mean(coverage)
  mean_ci_width  <- mean(ci_upper - ci_lower)

  list(
    boot_curves   = boot_curves,
    ci_lower      = ci_lower,
    ci_upper      = ci_upper,
    coverage      = coverage,
    mean_coverage = mean_coverage,
    mean_ci_width = mean_ci_width
  )
}

#' Run bootstrap copula PFS estimation for n iterations
#'
#' Repeatedly generates fresh datasets and applies the bootstrap copula method,
#' returning pointwise coverage, mean CI width, and mean survival estimate
#' across all iterations — analogous to run_iterations() in Simulations.R.
#'
#' @param n_times Integer. Number of follow-up time points.
#' @param n_patients Integer. Number of patients per simulated dataset.
#' @param n_iterations Integer. Number of outer simulation iterations.
#' @param mean Numeric vector. Mean vector for the multivariate normal
#'   distribution of log tumour-size ratios.
#' @param covariance Numeric matrix. Covariance matrix for the
#'   multivariate normal distribution.
#' @param alpha Numeric scalar. Intercept for the logistic new-lesion model.
#' @param beta Numeric scalar. Treatment-effect parameter.
#' @param gamma Numeric scalar. Tumour-size effect parameter.
#' @param R Numeric scalar or vector. Treatment indicator(s).
#' @param copula_family Integer. Copula family passed to BiCopSelect().
#' @param B Integer. Number of bootstrap resamples per iteration. Default 500.
#' @param threshold Numeric. Tumour growth threshold for progression. Default 1.2.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Rate}{Mean estimated PFS at each time point across iterations.}
#'     \item{CI_Width}{Mean bootstrap CI width at each time point.}
#'     \item{Coverage}{Proportion of CIs containing the true rate at each time point.}
#'   }
run_copula_iterations <- function(n_times, n_patients, n_iterations, mean, covariance,
                                  alpha, beta, gamma, R, copula_family,
                                  B = 500, threshold = 1.2) {

  # get true PFS from large cohort
  true_rate <- get_true_rates(n_times, n_true_patients = 100000, mean, covariance,
                              alpha, beta, gamma, R, threshold)$surv

  surv_prob_matrix <- matrix(NA, nrow = n_iterations, ncol = n_times)
  conf_low_matrix  <- matrix(NA, nrow = n_iterations, ncol = n_times)
  conf_high_matrix <- matrix(NA, nrow = n_iterations, ncol = n_times)

  for (i in seq_len(n_iterations)) {

    # generate fresh dataset
    coeffs <- generate_coefficients(n_times, n_patients, alpha, beta, gamma, R)
    cont   <- generate_continuous_data(n_times, n_patients, mean, covariance)
    Z      <- recover_tumour_sizes(cont$baseline_tumour_size,
                                   cont$log_tumour_size_ratio)
    p      <- compute_new_lesion_probability(coeffs, Z[, 1:n_times, drop = FALSE])
    L      <- generate_binary_data(p)

    # run bootstrap copula
    res <- tryCatch(
      bootstrap_copula_pfs(
        lesion_data      = L,
        tumour_size_data = Z,
        n_times          = n_times,
        copula_family    = copula_family,
        B                = B,
        alpha            = 0.05,
        true_pfs         = true_rate,
        threshold        = threshold,
        seed             = NULL
      ),
      error = function(e) NULL
    )

    if (!is.null(res)) {
      surv_prob_matrix[i, ] <- colMeans(res$boot_curves)
      conf_low_matrix[i, ]  <- res$ci_lower
      conf_high_matrix[i, ] <- res$ci_upper
    }
  }

  covered  <- sweep(conf_low_matrix,  2, true_rate, "<=") &
    sweep(conf_high_matrix, 2, true_rate, ">=")
  coverage <- colMeans(covered,  na.rm = TRUE)
  ci_width <- colMeans(conf_high_matrix - conf_low_matrix, na.rm = TRUE)
  mean_surv <- colMeans(surv_prob_matrix, na.rm = TRUE)

  return(data.frame(
    Rate     = round(mean_surv, 3),
    CI_Width = round(ci_width,  3),
    Coverage = round(coverage,  3)
  ))
}
