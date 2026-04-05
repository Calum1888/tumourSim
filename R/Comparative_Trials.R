#'Perform the log rank test for the control and treatment arms
#'
#' @param control_data time and event data for the control arm in the form for KM estimate
#' @param treatment_data time and event data for the treatment arm in the form for KM estimate
#'
#' @return p-value for the log rank test
#'
#' @importFrom survival survdiff Surv
#' @importFrom stats pchisq
#' @export
km_log_rank_test <- function(control_data, treatment_data){

  test_data <- data.frame(
    time      = c(control_data$time,
                  treatment_data$time),
    status    = c(control_data$status,
                  treatment_data$status),
    treatment = c(rep(0, nrow(control_data)),
                  rep(1, nrow(treatment_data))))

  log_rank_test <- survdiff(Surv(time, status)~treatment,data = test_data)

  p_val <- 1 - pchisq(log_rank_test$chisq, df =1)

  return(p_val)

}

#' Simulate one arm's data using the existing framework
#'
#' @inheritParams run_single_simulation
#' @return list with lesion matrix L and tumour size matrix Z
simulate_arm <- function(n_times, n_patients, mean, covariance,
                         alpha, beta, gamma, R, threshold = 1.2) {
  coeffs <- generate_coefficients(n_times, n_patients, alpha, beta, gamma, R)
  cont   <- generate_continuous_data(n_times, n_patients, mean, covariance)
  Z      <- recover_tumour_sizes(cont$baseline_tumour_size,
                                 cont$log_tumour_size_ratio)
  p      <- compute_new_lesion_probability(coeffs, Z[, 1:n_times, drop = FALSE])
  L      <- generate_binary_data(p)
  list(L = L, Z = Z)
}

#' Compute copula PFS curve for one arm
arm_copula_pfs <- function(arm, n_times, copula_family, threshold = 1.2) {
  le <- lesion_event(arm$L)
  te <- tumour_event(arm$Z, threshold = threshold)
  copula_pfs(le, te, n_times, copula_family)
}

#' Power of the Frank copula method for detecting survival differences
#'
#' Simulates control and treatment arms using the existing data-generation
#' framework, estimates PFS via the Frank copula for each arm, and
#' computes power as the proportion of iterations where the permutation
#' p-value (based on AUC difference) is below alpha.
#'
#' @param n_times Integer. Number of follow-up time points.
#' @param n_patients Integer. Number of patients per arm.
#' @param n_iterations Integer. Number of Monte Carlo iterations.
#' @param mean_ctrl Numeric vector (length n_times). MVN mean for control arm.
#' @param mean_trt  Numeric vector (length n_times). MVN mean for treatment arm.
#' @param covariance Numeric matrix (n_times x n_times). Shared covariance.
#' @param alpha_coef Numeric. Logistic model intercept (shared).
#' @param beta_ctrl  Numeric. Beta for control arm (typically 0).
#' @param beta_trt   Numeric. Beta for treatment arm.
#' @param gamma      Numeric. Tumour-size effect in logistic model.
#' @param copula_family Integer. Frank copula = 5 in VineCopula.
#' @param n_perm Integer. Number of permutations for the permutation test.
#' @param alpha_level Numeric. Significance level. Default 0.05.
#' @param threshold Numeric. Tumour growth threshold. Default 1.2.
#' @param seed Integer or NULL.
#'
#' @return A list with:
#'   \item{power}{Scalar. Proportion of iterations rejecting H0.}
#'   \item{mean_auc_diff}{Mean AUC difference (trt - ctrl) across iterations.}
#'   \item{p_values}{Numeric vector of per-iteration p-values.}
#'
#' @export
frank_copula_power <- function(n_times,
                               n_patients,
                               n_iterations,
                               mean_ctrl,
                               mean_trt,
                               covariance,
                               alpha_coef,
                               beta_ctrl  = 0,
                               beta_trt,
                               gamma,
                               copula_family = 5,    # 5 = Frank in VineCopula
                               n_perm        = 500,
                               alpha_level   = 0.05,
                               threshold     = 1.2,
                               seed          = NULL) {

  if (!is.null(seed)) set.seed(seed)

  p_values  <- numeric(n_iterations)
  auc_diffs <- numeric(n_iterations)

  observed_auc_diff <- function(pfs_ctrl, pfs_trt) {
    mean(pfs_trt - pfs_ctrl)   # signed AUC difference
  }

  for (i in seq_len(n_iterations)) {

    # --- generate both arms ---
    ctrl <- simulate_arm(n_times, n_patients, mean_ctrl, covariance,
                         alpha_coef, beta_ctrl, gamma, R = 0, threshold)
    trt  <- simulate_arm(n_times, n_patients, mean_trt,  covariance,
                         alpha_coef, beta_trt,  gamma, R = 1, threshold)

    pfs_ctrl <- tryCatch(arm_copula_pfs(ctrl, n_times, copula_family, threshold),
                         error = function(e) NULL)
    pfs_trt  <- tryCatch(arm_copula_pfs(trt,  n_times, copula_family, threshold),
                         error = function(e) NULL)

    if (is.null(pfs_ctrl) || is.null(pfs_trt)) {
      p_values[i]  <- NA
      auc_diffs[i] <- NA
      next
    }

    obs_diff <- observed_auc_diff(pfs_ctrl, pfs_trt)
    auc_diffs[i] <- obs_diff

    # --- permutation test: pool arms, reshuffle, recompute AUC diff ---
    perm_diffs <- numeric(n_perm)

    pool_L <- rbind(ctrl$L, trt$L)
    pool_Z <- rbind(ctrl$Z, trt$Z)
    n_pool <- 2 * n_patients

    for (j in seq_len(n_perm)) {
      idx       <- sample(n_pool)
      perm_ctrl <- list(L = pool_L[idx[1:n_patients], ],
                        Z = pool_Z[idx[1:n_patients], ])
      perm_trt  <- list(L = pool_L[idx[(n_patients+1):n_pool], ],
                        Z = pool_Z[idx[(n_patients+1):n_pool], ])

      pc <- tryCatch(arm_copula_pfs(perm_ctrl, n_times, copula_family, threshold),
                     error = function(e) NULL)
      pt <- tryCatch(arm_copula_pfs(perm_trt,  n_times, copula_family, threshold),
                     error = function(e) NULL)

      perm_diffs[j] <- if (!is.null(pc) && !is.null(pt))
        observed_auc_diff(pc, pt) else NA
    }

    perm_diffs    <- perm_diffs[!is.na(perm_diffs)]
    p_values[i]   <- mean(abs(perm_diffs) >= abs(obs_diff))
  }

  valid <- !is.na(p_values)

  list(
    power         = mean(p_values[valid] < alpha_level),
    mean_auc_diff = mean(auc_diffs[valid]),
    p_values      = p_values
  )
}

#' Power of log-rank test using the same data-generation framework
#'
#' @param n_times,n_patients,n_iterations,mean_ctrl,mean_trt,covariance,
#'   alpha_coef,beta_ctrl,beta_trt,gamma,threshold,seed  same as frank_copula_power
#' @return list: power, p_values
#'
#' @export
logrank_power <- function(n_times,
                          n_patients,
                          n_iterations,
                          mean_ctrl,
                          mean_trt,
                          covariance,
                          alpha_coef,
                          beta_ctrl  = 0,
                          beta_trt,
                          gamma,
                          alpha_level = 0.05,
                          threshold   = 1.2,
                          seed        = NULL) {

  if (!is.null(seed)) set.seed(seed)

  p_values <- numeric(n_iterations)

  for (i in seq_len(n_iterations)) {

    ctrl <- simulate_arm(n_times, n_patients, mean_ctrl, covariance,
                         alpha_coef, beta_ctrl, gamma, R = 0, threshold)
    trt  <- simulate_arm(n_times, n_patients, mean_trt,  covariance,
                         alpha_coef, beta_trt,  gamma, R = 1, threshold)

    # combined PFS event (lesion OR tumour progression)
    events_ctrl <- event_definition(ctrl$L, ctrl$Z, threshold)
    events_trt  <- event_definition(trt$L,  trt$Z,  threshold)

    p_values[i] <- tryCatch(
      km_log_rank_test(events_ctrl, events_trt),
      error = function(e) NA_real_
    )
  }

  valid <- !is.na(p_values)
  list(
    power    = mean(p_values[valid] < alpha_level),
    p_values = p_values
  )
}
