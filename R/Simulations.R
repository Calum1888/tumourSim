#' Run a single tumour‑progression simulation replicate
#'
#' Executes one full iteration of the tumour‑progression simulation
#' pipeline, including coefficient generation, tumour‑size simulation,
#' tumour‑size recovery, new‑lesion probability computation, binary
#' lesion generation, event definition, and Kaplan–Meier estimation.
#'
#' This function represents one Monte Carlo replicate and is intended
#' to be called repeatedly inside a simulation loop.
#'
#' @param n_times Integer. Number of follow‑up time points.
#' @param n_patients Integer. Number of patients to simulate.
#' @param mean Numeric vector. Mean vector for the multivariate normal
#'   distribution of log tumour‑size ratios.
#' @param covariance Numeric matrix. Covariance matrix for the
#'   multivariate normal distribution.
#' @param alpha Numeric scalar. Baseline intercept parameter for the
#'   logistic new‑lesion model.
#' @param beta Numeric scalar. Treatment‑effect parameter for the
#'   logistic new‑lesion model.
#' @param gamma Numeric scalar. Effect of previous tumour size in the
#'   logistic new‑lesion model.
#' @param R Numeric scalar or vector. Treatment indicator(s) for
#'   patients; defaults to 0 (control arm).
#' @param threshold Numeric. Threshold multiplier for defining tumour
#'   growth progression (default 1.2).
#'
#' @return A \code{summary.survfit} object containing survival
#'   probabilities and confidence intervals at times
#'   \code{1:n_times}.
#'
#' @details
#' The simulation proceeds through the following steps:
#' \enumerate{
#'   \item Generate logistic‑model coefficients.
#'   \item Simulate baseline tumour sizes and log tumour‑size ratios.
#'   \item Recover tumour sizes over time using the recursive nadir rule.
#'   \item Compute new‑lesion probabilities at each time point.
#'   \item Generate binary new‑lesion indicators.
#'   \item Define progression events based on lesions or tumour growth.
#'   \item Fit a Kaplan–Meier curve and extract survival estimates.
#' }
#'
#' @seealso \code{generate_coefficients()},
#'   \code{generate_continuous_data()},
#'   \code{recover_tumour_sizes()},
#'   \code{compute_new_lesion_probability()},
#'   \code{generate_binary_data()},
#'   \code{event_definition()},
#'   \code{survfit()}
run_single_simulation <- function(n_times, n_patients, mean, covariance, alpha, beta, gamma, R, threshold = 1.2) {

  # generate coefficients
  coeffs <- generate_coefficients(n_times, n_patients, alpha, beta, gamma, R)
  # simulate tumour‑size data
  cont <- generate_continuous_data(n_patients, mean, covariance)
  # recover tumour sizes
  tumour_sizes <- recover_tumour_sizes(cont$baseline_tumour_size,
                                       cont$log_tumour_size_ratio)
  # compute new‑lesion probabilities
  p <- compute_new_lesion_probability(coeffs,
                                      tumour_sizes[, 1:n_times, drop = FALSE])
  # generate binary lesions
  lesions <- generate_binary_data(p)
  # define events
  events <- event_definition(lesion_data = lesions,
                             tumour_size_data = tumour_sizes,
                             threshold = threshold)
  # Kaplan–Meier fit
  fit <- survfit(Surv(time, status) ~ 1, data = events)
  return(summary(fit, times = seq_len(n_times), extend = TRUE))
}

#' Executes \code{n_iterations}s of the single Kaplan-Meier iteration in the function \code{run_single_simulation}.
#'
#' @param n_times Integer. Number of follow‑up time points.
#' @param n_patients Integer. Number of patients to simulate.
#' @param mean Numeric vector. Mean vector for the multivariate normal
#'   distribution of log tumour‑size ratios.
#' @param covariance Numeric matrix. Covariance matrix for the
#'   multivariate normal distribution.
#' @param alpha Numeric scalar. Baseline intercept parameter for the
#'   logistic new‑lesion model.
#' @param beta Numeric scalar. Treatment‑effect parameter for the
#'   logistic new‑lesion model.
#' @param gamma Numeric scalar. Effect of previous tumour size in the
#'   logistic new‑lesion model.
#' @param R Numeric scalar or vector. Treatment indicator(s) for
#'   patients; defaults to 0 (control arm).
#' @param threshold Numeric. Threshold multiplier for defining tumour
#'   growth progression (default 1.2).
#'
#' @export
run_iterations <- function(n_times, n_patients, n_iterations, mean, covariance, alpha, beta, gamma, R, threshold = 1.2){

  surv_prob_matrix <- matrix(NA, nrow = n_iterations, ncol = n_times)
  conf_low_matrix  <- matrix(NA, nrow = n_iterations, ncol = n_times)
  conf_high_matrix <- matrix(NA, nrow = n_iterations, ncol = n_times)

  for (i in seq_len(n_iterations)) {
    fit_summary <- run_single_simulation(
      n_times     = n_times,
      n_patients  = n_patients,
      mean        = mean,
      covariance  = covariance,
      alpha       = alpha,
      beta        = beta,
      gamma       = gamma,
      R           = R,
      threshold   = threshold
    )
    surv_prob_matrix[i, ] <- fit_summary$surv
    conf_low_matrix[i, ]  <- fit_summary$lower
    conf_high_matrix[i, ] <- fit_summary$upper
  }

  covered <- (sweep(conf_low_matrix, 2, true_rate, "<=") & sweep(conf_high_matrix, 2, true_rate, ">="))
  coverage <- colMeans(covered, na.rm = TRUE)
  ci_width <- colMeans(conf_high_matrix - conf_low_matrix)
  mean_surv <- colMeans(surv_prob_matrix)
  return(data.frame(Rate = round(mean_surv, 3),
                    Coverage = coverage,
                    CI_Width = round(ci_width, 3)
                    ))
}
