#' Estimate power of the copula-based PFS test via Monte Carlo simulation
#'
#' Compares copula-based PFS curves between a treated (R=1) and control (R=0)
#' arm across n_iterations simulated trials. At each time point, a Z-test is
#' formed using bootstrap standard errors. Power is the proportion of
#' iterations in which H0 is rejected.
#'
#' @param n_times        Integer. Number of follow-up time points.
#' @param n_patients     Integer. Patients per arm.
#' @param n_iterations   Integer. Number of Monte Carlo outer iterations.
#' @param mean           Numeric vector (length n_times). MVN mean for log tumour ratios.
#' @param covariance     Numeric matrix (n_times x n_times). MVN covariance.
#' @param alpha_coef     Numeric. Logistic model intercept.
#' @param beta           Numeric. Treatment effect (set to 0 to estimate type I error).
#' @param gamma          Numeric. Tumour size effect in lesion model.
#' @param copula_family  Integer. Copula family passed to BiCopSelect().
#' @param B              Integer. Bootstrap resamples per arm per iteration.
#' @param alpha_level    Numeric. Significance level (default 0.05).
#' @param threshold      Numeric. Tumour growth threshold (default 1.2).
#' @param seed           Optional integer seed.
#'
#' @return A data frame with columns:
#'   Time       - follow-up time point
#'   Power      - estimated power at that time point
#'   TypeI      - estimated type I error (run separately with beta=0)
#'   MeanDiff   - mean estimated PFS difference (treated - control)
#'
#'   @export
power_copula_pfs <- function(n_times, n_patients, n_iterations,
                             mean, covariance,
                             alpha_coef, beta, gamma,
                             copula_family, B = 200,
                             alpha_level = 0.05,
                             threshold = 1.2,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  reject_matrix <- matrix(NA, nrow = n_iterations, ncol = n_times)
  diff_matrix   <- matrix(NA, nrow = n_iterations, ncol = n_times)

  for (i in seq_len(n_iterations)) {

    # ---- Simulate control arm (R = 0) ----------------------------------------
    res0 <- tryCatch({
      coeffs0 <- generate_coefficients(n_times, n_patients, alpha_coef, beta=0,    gamma, R=0)
      cont0   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z0      <- recover_tumour_sizes(cont0$baseline_tumour_size, cont0$log_tumour_size_ratio)
      p0      <- compute_new_lesion_probability(coeffs0, Z0[, 1:n_times, drop=FALSE])
      L0      <- generate_binary_data(p0)
      bootstrap_copula_pfs(L0, Z0, n_times, copula_family, B=B,
                           alpha_level=alpha_level, true_pfs=rep(0.5, n_times),
                           threshold=threshold)
    }, error = function(e) NULL)

    # ---- Simulate treated arm (R = 1, beta active) ---------------------------
    res1 <- tryCatch({
      coeffs1 <- generate_coefficients(n_times, n_patients, alpha_coef, beta,      gamma, R=1)
      cont1   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z1      <- recover_tumour_sizes(cont1$baseline_tumour_size, cont1$log_tumour_size_ratio)
      p1      <- compute_new_lesion_probability(coeffs1, Z1[, 1:n_times, drop=FALSE])
      L1      <- generate_binary_data(p1)
      bootstrap_copula_pfs(L1, Z1, n_times, copula_family, B=B,
                           alpha_level=alpha_level, true_pfs=rep(0.5, n_times),
                           threshold=threshold)
    }, error = function(e) NULL)

    if (is.null(res0) || is.null(res1)) next

    # ---- Pointwise Z-test at each time t -------------------------------------
    pfs0 <- colMeans(res0$boot_curves)
    pfs1 <- colMeans(res1$boot_curves)

    # Bootstrap SE from CI width (width = 2 * z_{alpha/2} * SE)
    z_crit <- qnorm(1 - alpha_level / 2)
    se0 <- (res0$ci_upper - res0$ci_lower) / (2 * z_crit)
    se1 <- (res1$ci_upper - res1$ci_lower) / (2 * z_crit)

    se_diff   <- sqrt(se0^2 + se1^2)
    Z_stat    <- (pfs1 - pfs0) / se_diff

    reject_matrix[i, ] <- as.integer(abs(Z_stat) > z_crit)
    diff_matrix[i, ]   <- pfs1 - pfs0
  }

  power    <- colMeans(reject_matrix, na.rm = TRUE)
  mean_diff <- colMeans(diff_matrix,  na.rm = TRUE)

  data.frame(
    Time      = seq_len(n_times),
    Power     = round(power,     3),
    MeanDiff  = round(mean_diff, 3)
  )
}
