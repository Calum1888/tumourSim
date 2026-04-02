library(MASS)
library(docstring)

generate_coefficients <- function(n_times, n_patients,alpha, beta, gamma, R){
  #'Generate coefficient vectors for logistic regression
  #'
  #' Creates repeated vectors of the parameters \code{alpha}, \code{beta},
  #' \code{gamma}, and \code{R} for use in logistic regression models.
  #' Each of \code{alpha}, \code{beta}, and \code{gamma} is repeated
  #' \code{n_times} times, while \code{R} is repeated \code{n_patients} times.
  #'
  #' @param n_times Integer. Number of time points.
  #' @param n_patients Integer. Number of patients.
  #' @param alpha Numeric scalar. Baseline alpha coefficient.
  #' @param beta Numeric scalar. Baseline beta coefficient.
  #' @param gamma Numeric scalar. Baseline gamma coefficient.
  #' @param R Numeric scalar. Indicator or coefficient repeated per patient.
  #'
  #' @return A list containing:
  #' \describe{
  #'   \item{alpha_r}{Numeric vector of length \code{n_times}.}
  #'   \item{beta_r}{Numeric vector of length \code{n_times}.}
  #'   \item{gamma_r}{Numeric vector of length \code{n_times}.}
  #'   \item{R_indicator}{Numeric vector of length \code{n_patients}.}
  #' }
    
  alpha_vec <- rep(alpha, n_times)
  beta_vec <- rep(beta, n_times)
  gamma_vec <- rep(gamma, n_times)
  R_vec <- rep(R, n_patients)
  
  return(list(alpha_r = alpha_vec,
              beta_r = beta_vec,
              gamma_r = gamma_vec,
              R_indicator = R_vec))
}

generate_continuous_data <- function(n_patients, mean, covariance){
  #' Generate continuous tumour‑size data
  #'
  #' Simulates tumour‑size data for a cohort of patients by drawing
  #' multivariate normal log tumour size ratios and uniform baseline
  #' tumour sizes. The log tumour size ratios are sampled from a
  #' multivariate normal distribution with specified mean vector and
  #' covariance matrix, while baseline tumour sizes are sampled from a
  #' uniform distribution on \eqn{[0, 1]}.
  #'
  #' @param n_patients Integer. Number of patients to simulate.
  #' @param mean Numeric vector. Mean vector for the multivariate normal
  #'   distribution of log tumour size ratios.
  #' @param covariance Numeric matrix. Covariance matrix for the
  #'   multivariate normal distribution.
  #'
  #' @return A list containing:
  #' \describe{
  #'   \item{log_tumour_size_ratio}{An \code{n_patients × length(mean)} matrix
  #'     of multivariate normal samples.}
  #'   \item{baseline_tumour_size}{A numeric vector of length
  #'     \code{n_patients} with values drawn from \eqn{U(0,1)}.}
  #' }
  
  ltsr <- mvrnorm(n = n_patients, mu = mean, Sigma = covariance)
  # min =0, max = 1 remain the same no matter other changes in variables
  bts <- runif(n = n_patients, min = 0, max = 1)
  
  return(list(log_tumour_size_ratio = ltsr,
              baseline_tumour_size = bts))
} 

recover_tumour_sizes <- function(baseline, ratios){
  #' Recover tumour sizes from log-ratio trajectories
  #'
  #' Computes tumour sizes over time using the recursive update
  #' \deqn{z_{it} = \min(z_{i0}, \ldots, z_{i(t-1)}) \exp(\nu_{it})}
  #' where \eqn{z_{i0}} is the baseline tumour size and \eqn{\nu_{it}}
  #' is the log tumour size ratio at time \eqn{t}.
  #'
  #' @param baseline Numeric vector of baseline tumour sizes for each patient.
  #' @param ratios Numeric matrix of log tumour size ratios, with
  #'   \code{nrow(ratios)} equal to the number of patients and
  #'   \code{ncol(ratios)} equal to the number of time points.
  #'
  #' @return A numeric matrix of dimension
  #'   \code{n_patients × (n_times + 1)}, containing tumour sizes at
  #'   baseline and all subsequent time points.
  
  n <- length(baseline)
  T <- ncol(ratios)
  # creating n x (T+1) matrix includes baseline size in Z
  Z <- matrix(NA, n, T + 1)
  Z[, 1] <- baseline
  
  for (t in 1:T) {
    nadir <- apply(Z[, 1:t, drop = FALSE], 1, min)
    Z[, t + 1] <- nadir * exp(ratios[, t])
  }
  return(Z)
}

compute_new_lesion_probability <- function(coefficients, tumour_sizes) {
  #' Compute probability of new lesion appearance
  #'
  #' @param coefficients List returned by generate_coefficients()
  #' @param tumour_sizes Matrix of previous tumour sizes (n × T)
  #'
  #' @return Matrix of probabilities p_it
  
  alpha <- coefficients$alpha_r
  beta  <- coefficients$beta_r
  gamma <- coefficients$gamma_r
  R     <- coefficients$R_indicator
  
  n <- length(R)
  T <- length(alpha)
  
  linpred <- matrix(alpha, n, T, byrow = TRUE) +
    matrix(beta,  n, T, byrow = TRUE) * R +
    matrix(gamma, n, T, byrow = TRUE) * tumour_sizes
  
  p <- 1 / (1 + exp(-linpred))
  return(p)
}

generate_binary_data <- function(probabilities) {
  #' Generate binary new-lesion indicators from probabilities
  #'
  #' @param probabilities Matrix (n × T) of probabilities p_it
  #'
  #' @return Matrix (n × T) of 0/1 lesion indicators
  
  n <- nrow(probabilities)
  T <- ncol(probabilities)
  
  # rbinom needs a vector, so flatten then reshape
  draws <- rbinom(n * T, size = 1, prob = as.vector(probabilities))
  
  lesion <- matrix(draws, nrow = n, ncol = T)
  
  return(lesion)
}
event_definition <- function(lesion_data, tumour_size_data, threshold = 1.2) {
  #' Define progression events based on new lesions or tumour growth
  #'
  #' @param lesion_data Matrix (n × T) of 0/1 new-lesion indicators
  #' @param tumour_size_data Matrix (n × (T+1)) of tumour sizes including baseline
  #' @param threshold Numeric. Growth threshold for progression (default 1.2)
  #'
  #' @return Data frame with time and status for each patient
  
  n_patients <- nrow(tumour_size_data)
  n_times    <- ncol(tumour_size_data) - 1
  
  # tumour sizes at times 1..T
  tumour <- tumour_size_data[, -1, drop = FALSE]
  rolling_min <- t(apply(tumour_size_data, 1, cummin))[, -ncol(tumour_size_data), drop = FALSE]
  event_matrix <- (lesion_data == 1) | (tumour > threshold * rolling_min)
  
  # find first TRUE per row using max.col on reversed matrix
  # max.col returns index of first max; reversing makes TRUE appear last
  idx <- max.col(event_matrix, ties.method = "first")
  
  # if no event, row is all FALSE → max.col returns 1, but event_mat[,1] is FALSE
  no_event <- rowSums(event_matrix) == 0
  
  time   <- idx
  time[no_event] <- n_times
  
  status <- as.integer(!no_event)
  
  data.frame(time = time, status = status)
}

  