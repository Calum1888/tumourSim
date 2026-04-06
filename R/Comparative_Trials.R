#' Perform the log rank test for the control and treatment arms
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

  log_rank_test <- survdiff(Surv(time, status) ~ treatment, data = test_data)

  p_val <- 1 - pchisq(log_rank_test$chisq, df = 1)

  return(p_val)

}

#' Simulate one arm's data using the existing framework
#'
#' @inheritParams run_single_simulation
#' @return list with lesion matrix L and tumour size matrix Z
#'
#' @export
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
#'
#' @export
arm_copula_pfs <- function(arm, n_times, copula_family, threshold = 1.2) {
  le <- lesion_event(arm$L)
  te <- tumour_event(arm$Z, threshold = threshold)
  copula_pfs(le, te, n_times, copula_family)
}

#' Power of the Frank copula method using a patient-level permutation test
#'
#' For each iteration, simulates control and treatment arms, computes the
#' observed trapezoidal AUC difference between copula PFS curves, then builds
#' a null distribution by randomly reassigning patients to arms and recomputing
#' the copula PFS curves from scratch each time.
#'
#' Changes vs previous version:
#' \itemize{
#'   \item AUC computed via trapezoidal rule (weighted by time spacing) rather
#'         than a simple unweighted mean across time points.
#'   \item One-sided permutation p-value (treatment better than control) rather
#'         than two-sided, increasing power when direction of effect is known.
#'   \item Per-iteration seeding (\code{seed + i}) so each iteration is
#'         independently reproducible and both methods see the same data when
#'         called with the same seed.
#'   \item Default \code{n_perm} raised to 1000 for stable p-value resolution
#'         at \code{alpha_level = 0.05}.
#' }
#'
#' @param n_times       number of assessment time points
#' @param n_patients    patients per arm
#' @param n_iterations  number of simulated trials
#' @param mean_ctrl     mean tumour trajectory for control arm (length n_times)
#' @param mean_trt      mean tumour trajectory for treatment arm (length n_times)
#' @param covariance    covariance matrix for tumour trajectories (n_times x n_times)
#' @param alpha_coef    baseline hazard intercept
#' @param beta_ctrl     log-hazard coefficient for control arm (default 0)
#' @param beta_trt      log-hazard coefficient for treatment arm
#' @param gamma         tumour burden effect on hazard
#' @param times         numeric vector of assessment times used for trapezoidal
#'                      AUC (default \code{seq_len(n_times)})
#' @param copula_family integer copula family passed to \code{arm_copula_pfs}
#'                      (5 = Frank)
#' @param alpha_level   significance threshold (default 0.05)
#' @param threshold     tumour progression threshold (default 1.2)
#' @param n_perm        permutations per iteration (default 1000)
#' @param seed          integer base seed; iteration i uses \code{seed + i} for
#'                      independent reproducibility. NULL disables seeding.
#'
#' @return list: power, mean_auc_diff, p_values
#'
#' @export
frank_copula_power <- function(n_times,
                               n_patients,
                               n_iterations,
                               mean_ctrl,
                               mean_trt,
                               covariance,
                               alpha_coef,
                               beta_ctrl     = 0,
                               beta_trt,
                               gamma,
                               times         = seq_len(n_times),
                               copula_family = 5,
                               alpha_level   = 0.05,
                               threshold     = 1.2,
                               n_perm        = 1000,
                               seed          = NULL) {

  # Helper: trapezoidal AUC of the difference curve
  auc_diff <- function(pfs_a, pfs_b) {
    diff_curve <- pfs_a - pfs_b
    # trapz: sum of (delta_t * average height of adjacent points)
    dt <- diff(times)
    sum(dt * (head(diff_curve, -1) + tail(diff_curve, -1)) / 2)
  }

  p_values  <- rep(NA_real_, n_iterations)
  auc_diffs <- rep(NA_real_, n_iterations)

  for (i in seq_len(n_iterations)) {

    # FIX 3: per-iteration seed so both methods see the same data
    if (!is.null(seed)) set.seed(seed + i)

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

    # FIX 1: trapezoidal AUC instead of unweighted mean
    obs_diff      <- auc_diff(pfs_trt, pfs_ctrl)
    auc_diffs[i]  <- obs_diff

    # Pool patient-level data
    L_pool  <- rbind(ctrl$L, trt$L)
    Z_pool  <- rbind(ctrl$Z, trt$Z)
    n_total <- nrow(L_pool)

    perm_diffs <- replicate(n_perm, {
      idx    <- sample.int(n_total, n_patients, replace = FALSE)
      perm_a <- list(L = L_pool[ idx, , drop = FALSE],
                     Z = Z_pool[ idx, , drop = FALSE])
      perm_b <- list(L = L_pool[-idx, , drop = FALSE],
                     Z = Z_pool[-idx, , drop = FALSE])

      pfs_a <- tryCatch(
        arm_copula_pfs(perm_a, n_times, copula_family, threshold),
        error = function(e) NULL
      )
      pfs_b <- tryCatch(
        arm_copula_pfs(perm_b, n_times, copula_family, threshold),
        error = function(e) NULL
      )

      if (is.null(pfs_a) || is.null(pfs_b)) return(NA_real_)
      # FIX 1: trapezoidal AUC for permuted differences too
      auc_diff(pfs_a, pfs_b)
    })

    perm_diffs <- perm_diffs[!is.na(perm_diffs)]
    if (length(perm_diffs) == 0) next

    # FIX 2: one-sided p-value (treatment expected to be better, i.e. obs_diff > 0)
    p_values[i] <- mean(perm_diffs >= obs_diff)

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
#' Updated to use per-iteration seeding (\code{seed + i}) to match
#' \code{frank_copula_power}, ensuring both methods simulate from the same
#' datasets when called with the same seed.
#'
#' @param n_times,n_patients,n_iterations,mean_ctrl,mean_trt,covariance,
#'   alpha_coef,beta_ctrl,beta_trt,gamma,threshold  same as frank_copula_power
#' @param alpha_level significance threshold (default 0.05)
#' @param seed integer base seed; iteration i uses \code{seed + i}.
#'             NULL disables seeding.
#'
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
                          beta_ctrl   = 0,
                          beta_trt,
                          gamma,
                          alpha_level = 0.05,
                          threshold   = 1.2,
                          seed        = NULL) {

  p_values <- numeric(n_iterations)

  for (i in seq_len(n_iterations)) {

    # FIX 3: per-iteration seed to match frank_copula_power
    if (!is.null(seed)) set.seed(seed + i)

    ctrl <- simulate_arm(n_times, n_patients, mean_ctrl, covariance,
                         alpha_coef, beta_ctrl, gamma, R = 0, threshold)
    trt  <- simulate_arm(n_times, n_patients, mean_trt,  covariance,
                         alpha_coef, beta_trt,  gamma, R = 1, threshold)

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
