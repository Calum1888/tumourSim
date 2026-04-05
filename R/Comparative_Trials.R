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

#' Power of the Frank copula method using a paired Wilcoxon test
#'
#' Replaces the permutation test with a paired Wilcoxon signed-rank test on
#' the pointwise PFS differences across time points. Much faster than
#' permutation — no inner loop required.
#'
#' @param n_times,n_patients,n_iterations,mean_ctrl,mean_trt,covariance,
#'   alpha_coef,beta_ctrl,beta_trt,gamma,alpha_level,threshold,seed
#'   same as before
#'
#' @return list: power, mean_auc_diff, p_values
#'
#' @importFrom stats wilcox.test
#' @export
frank_copula_power <- function(n_times,
                                      n_patients,
                                      n_iterations,
                                      mean_ctrl,
                                      mean_trt,
                                      covariance,
                                      alpha_coef,
                                      beta_ctrl   = 0,
                                      beta_trt,
                                      gamma,
                                      copula_family = 5,
                                      alpha_level   = 0.05,
                                      threshold     = 1.2,
                                      seed          = NULL) {


  p_values  <- rep(NA_real_, n_iterations)
  auc_diffs <- rep(NA_real_, n_iterations)

  for (i in seq_len(n_iterations)) {

    ctrl <- simulate_arm(n_times, n_patients, mean_ctrl, covariance,
                         alpha_coef, beta_ctrl, gamma, R = 0, threshold)
    trt  <- simulate_arm(n_times, n_patients, mean_trt,  covariance,
                         alpha_coef, beta_trt,  gamma, R = 1, threshold)

    pfs_ctrl <- tryCatch(
      arm_copula_pfs(ctrl, n_times, copula_family, threshold),
      error = function(e) NULL
    )
    pfs_trt <- tryCatch(
      arm_copula_pfs(trt, n_times, copula_family, threshold),
      error = function(e) NULL
    )

    if (is.null(pfs_ctrl) || is.null(pfs_trt)) next

    auc_diffs[i] <- mean(pfs_trt - pfs_ctrl)

    # paired Wilcoxon on pointwise differences across the 5 time points
    p_values[i] <- tryCatch(
      wilcox.test(pfs_trt, pfs_ctrl, paired = TRUE, exact = FALSE)$p.value,
      error   = function(e) NA_real_,
      warning = function(w) wilcox.test(pfs_trt, pfs_ctrl,
                                        paired = TRUE, exact = FALSE)$p.value
    )

    if (i %% 100 == 0)
      cat(sprintf("  iteration %d / %d\n", i, n_iterations))
  }

  valid <- !is.na(p_values)
  list(
    power         = mean(p_values[valid] < alpha_level),
    mean_auc_diff = mean(auc_diffs[valid], na.rm = TRUE),
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
