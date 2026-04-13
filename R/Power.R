#' Estimate power of the copula-based PFS test via Monte Carlo simulation
#'
#' Compares copula-based PFS curves between a treated (R=1) and control (R=0)
#' arm across n_iterations simulated trials. At each time point, a Z-test is
#' formed using bootstrap standard errors. Power is the proportion of
#' iterations in which H0 is rejected. Set beta = 0 to estimate type I error.
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
#'   Time      - follow-up time point
#'   Power     - rejection rate at that time point (= type I error when beta = 0)
#'   MeanDiff  - mean estimated PFS difference (treated - control)
#'
#' @details
#' Parallelisation is controlled by the caller via the \code{future}
#' backend. Call \code{future::plan(multisession)} before running
#' this function to enable parallel computation across iterations.
#' If no plan is set, execution falls back to sequential.
#'
#' @importFrom stats qnorm
#' @importFrom future.apply future_lapply
#' @export
power_copula_pfs <- function(n_times, n_patients, n_iterations,
                             mean, covariance,
                             alpha_coef, beta, gamma,
                             copula_family, B = 100,
                             alpha_level = 0.05,
                             threshold = 1.2,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  true_pfs0 <- get_true_rates(n_times, 100000, mean, covariance,
                              alpha_coef, beta, gamma, R=0, threshold)$surv
  true_pfs1 <- get_true_rates(n_times, 100000, mean, covariance,
                              alpha_coef, beta, gamma, R=1, threshold)$surv

  results <- future_lapply(seq_len(n_iterations), function(i) {

    res0 <- tryCatch({
      coeffs0 <- generate_coefficients(n_times, n_patients, alpha_coef, beta, gamma, R=0)
      cont0   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z0      <- recover_tumour_sizes(cont0$baseline_tumour_size, cont0$log_tumour_size_ratio)
      p0      <- compute_new_lesion_probability(coeffs0, Z0[, 1:n_times, drop=FALSE])
      L0      <- generate_binary_data(p0)
      bootstrap_copula_pfs(L0, Z0, n_times, copula_family, B=B,
                           alpha_level=alpha_level, true_pfs=true_pfs0,
                           threshold=threshold)
    }, error = function(e) NULL)

    res1 <- tryCatch({
      coeffs1 <- generate_coefficients(n_times, n_patients, alpha_coef, beta, gamma, R=1)
      cont1   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z1      <- recover_tumour_sizes(cont1$baseline_tumour_size, cont1$log_tumour_size_ratio)
      p1      <- compute_new_lesion_probability(coeffs1, Z1[, 1:n_times, drop=FALSE])
      L1      <- generate_binary_data(p1)
      bootstrap_copula_pfs(L1, Z1, n_times, copula_family, B=B,
                           alpha_level=alpha_level, true_pfs=true_pfs1,
                           threshold=threshold)
    }, error = function(e) NULL)

    if (is.null(res0) || is.null(res1)) return(NULL)

    pfs0   <- colMeans(res0$boot_curves)
    pfs1   <- colMeans(res1$boot_curves)
    z_crit <- qnorm(1 - alpha_level / 2)
    se0    <- (res0$ci_upper - res0$ci_lower) / (2 * z_crit)
    se1    <- (res1$ci_upper - res1$ci_lower) / (2 * z_crit)
    Z_stat <- (pfs1 - pfs0) / sqrt(se0^2 + se1^2)

    list(
      reject = as.integer(abs(Z_stat) > z_crit),
      diff   = pfs1 - pfs0
    )
  }, future.seed = TRUE)

  results <- Filter(Negate(is.null), results)

  reject_matrix <- do.call(rbind, lapply(results, `[[`, "reject"))
  diff_matrix   <- do.call(rbind, lapply(results, `[[`, "diff"))

  data.frame(
    Time     = seq_len(n_times),
    Power    = round(colMeans(reject_matrix), 3),
    MeanDiff = round(colMeans(diff_matrix),   3)
  )
}

#' Estimate power of the log-rank test via Monte Carlo simulation
#'
#' Compares Kaplan-Meier PFS curves between a treated (R=1) and control (R=0)
#' arm across \code{n_iterations} simulated trials using a log-rank test.
#' Power is the proportion of iterations in which H0 is rejected at
#' \code{alpha_level}. Set \code{beta = 0} to estimate type I error.
#'
#' The output format mirrors \code{power_copula_pfs()} exactly, with a single
#' row (the log-rank test is global, not time-point-specific) reported at
#' \code{Time = NA} and a \code{MeanDiff} column giving the mean difference
#' in restricted mean survival time (RMST) between arms.
#'
#' @param n_times       Integer. Number of follow-up time points.
#' @param n_patients    Integer. Patients per arm.
#' @param n_iterations  Integer. Number of Monte Carlo outer iterations.
#' @param mean          Numeric vector (length \code{n_times}). MVN mean for
#'   log tumour-size ratios.
#' @param covariance    Numeric matrix (\code{n_times x n_times}). MVN
#'   covariance.
#' @param alpha_coef    Numeric. Logistic model intercept.
#' @param beta          Numeric. Treatment effect (set to 0 for type I error).
#' @param gamma         Numeric. Tumour size effect in lesion model.
#' @param alpha_level   Numeric. Significance level. Default \code{0.05}.
#' @param threshold     Numeric. Tumour growth threshold. Default \code{1.2}.
#' @param seed          Optional integer seed.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Power}{Proportion of iterations in which the log-rank test
#'       rejected H0.}
#'     \item{MeanDiff}{Mean difference in restricted mean survival time
#'       (RMST, truncated at \code{n_times}) between treated and control arms
#'       across iterations.}
#'   }
#'
#' @details
#' On each iteration both arms are simulated independently using the shared
#' data-generation pipeline (\code{generate_coefficients},
#' \code{generate_continuous_data}, \code{recover_tumour_sizes},
#' \code{compute_new_lesion_probability}, \code{generate_binary_data},
#' \code{event_definition}). The two event data frames are pooled with an
#' arm indicator and a log-rank test is applied via
#' \code{survival::survdiff()}. The p-value is extracted from the chi-squared
#' statistic and compared against \code{alpha_level}.
#'
#' RMST is computed as the area under each arm's KM curve up to
#' \code{n_times} using the trapezoidal rule, and the difference
#' (treated minus control) is recorded per iteration.
#'
#' Parallelisation is controlled by the caller via the \code{future} backend.
#' Call \code{future::plan(multisession)} before running this function to
#' enable parallel execution. If no plan is set, execution falls back to
#' sequential.
#'
#' @seealso \code{\link{power_copula_pfs}}
#'
#' @importFrom survival survdiff survfit Surv
#' @importFrom stats pchisq
#' @importFrom future.apply future_lapply
#' @importFrom utils head tail read.csv
#' @export
power_logrank_pfs <- function(n_times, n_patients, n_iterations,
                              mean, covariance,
                              alpha_coef, beta, gamma,
                              alpha_level = 0.05,
                              threshold   = 1.2,
                              seed        = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # ---- helper: RMST via trapezoid rule on KM curve ---------------------------
  km_rmst <- function(events, tau) {
    fit     <- survfit(Surv(time, status) ~ 1, data = events)
    summary <- summary(fit, times = seq_len(tau), extend = TRUE)
    times   <- c(0, summary$time)
    survs   <- c(1, summary$surv)
    sum(diff(times) * head(survs, -1))
  }

  # ---- simulation loop -------------------------------------------------------
  results <- future_lapply(seq_len(n_iterations), function(i) {

    tryCatch({

      # --- control arm (R = 0) ---
      coeffs0 <- generate_coefficients(n_times, n_patients, alpha_coef,
                                       beta, gamma, R = 0)
      cont0   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z0      <- recover_tumour_sizes(cont0$baseline_tumour_size,
                                      cont0$log_tumour_size_ratio)
      p0      <- compute_new_lesion_probability(coeffs0,
                                                Z0[, 1:n_times, drop = FALSE])
      L0      <- generate_binary_data(p0)
      events0 <- event_definition(L0, Z0, threshold = threshold)

      # --- treated arm (R = 1) ---
      coeffs1 <- generate_coefficients(n_times, n_patients, alpha_coef,
                                       beta, gamma, R = 1)
      cont1   <- generate_continuous_data(n_times, n_patients, mean, covariance)
      Z1      <- recover_tumour_sizes(cont1$baseline_tumour_size,
                                      cont1$log_tumour_size_ratio)
      p1      <- compute_new_lesion_probability(coeffs1,
                                                Z1[, 1:n_times, drop = FALSE])
      L1      <- generate_binary_data(p1)
      events1 <- event_definition(L1, Z1, threshold = threshold)

      # --- pool arms and run log-rank test ---
      events0$arm <- 0L
      events1$arm <- 1L
      pooled      <- rbind(events0, events1)

      lr      <- survdiff(Surv(time, status) ~ arm, data = pooled)
      p_value <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
      reject  <- as.integer(p_value < alpha_level)

      # --- RMST difference (treated - control) ---
      rmst0 <- km_rmst(events0, tau = n_times)
      rmst1 <- km_rmst(events1, tau = n_times)

      list(reject = reject, diff = rmst1 - rmst0)

    }, error = function(e) NULL)

  }, future.seed = TRUE)

  # ---- aggregate results -----------------------------------------------------
  results <- Filter(Negate(is.null), results)

  if (length(results) == 0) stop("All iterations failed.")

  reject_vec <- vapply(results, `[[`, integer(1), "reject")
  diff_vec   <- vapply(results, `[[`, numeric(1), "diff")

  data.frame(
    Power    = round(mean(reject_vec), 3),
    MeanDiff = round(mean(diff_vec),   3)
  )
}











