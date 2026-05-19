# =============================================================================
# Monte Carlo drivers for coverage and power
# =============================================================================
# This file glues together the data-generating functions, KM baseline,
# AugBin parametric model, and copula reconstruction into Monte Carlo
# studies for:
#
#   * Coverage/width/bias of the three PFS estimators (simulate_metrics_all,
#     evaluate_augbin_parallel).
#   * Power and type-I error of the three tests in a two-arm trial
#     (run_power_study, run_copula_power, run_scenario, run_all_scenarios).
# =============================================================================


#' Single-trial simulator returning a log-rank p-value
#'
#' Simulates a two-arm trial under a single scenario specification and
#' returns the log-rank p-value on the composite endpoint. Used as the
#' inner loop of \code{\link{run_power_study}}.
#'
#' @param s A scenario list. May contain any of \code{mean_vec},
#'   \code{covariance}, \code{alpha}, \code{beta}, \code{gamma},
#'   \code{cens_at}; missing entries fall back to the
#'   \code{defaults} argument.
#' @param n_per_arm Integer. Number of patients per arm.
#' @param defaults A list of defaults used when the scenario does not
#'   override a value. Must contain \code{mean}, \code{covariance},
#'   \code{threshold}, \code{alpha}, \code{gamma}, and optionally
#'   \code{cens_at}.
#'
#' @return Scalar log-rank p-value, or \code{NA_real_} if either arm had
#'   zero events.
#'
#' @importFrom survival survdiff Surv
#' @importFrom stats pchisq
#' @export
simulate_one_trial <- function(s, n_per_arm, defaults) {

  mean_vec <- if (!is.null(s$mean_vec))   s$mean_vec   else defaults$mean
  cov_mat  <- if (!is.null(s$covariance)) s$covariance else defaults$covariance
  alpha    <- if (!is.null(s$alpha))      s$alpha      else defaults$alpha
  gamma    <- if (!is.null(s$gamma))      s$gamma      else defaults$gamma
  cens_at  <- if (!is.null(s$cens_at))    s$cens_at    else (defaults$cens_at %||% Inf)
  beta     <- s$beta
  threshold <- defaults$threshold

  # Control arm
  t0 <- tumour_events(n_per_arm, mean_vec, cov_mat, threshold)
  l0 <- lesion_events(n_per_arm, t0$log_ratios_vs_baseline,
                      alpha = alpha, beta = beta, gamma = gamma,
                      treatment_arm = 0)
  ce0 <- combined_event(l0$events, t0$events); ce0$arm <- 0L

  # Treatment arm
  t1 <- tumour_events(n_per_arm, mean_vec, cov_mat, threshold)
  l1 <- lesion_events(n_per_arm, t1$log_ratios_vs_baseline,
                      alpha = alpha, beta = beta, gamma = gamma,
                      treatment_arm = 1)
  ce1 <- combined_event(l1$events, t1$events); ce1$arm <- 1L

  ce0 <- apply_admin_censoring(ce0, cens_at)
  ce1 <- apply_admin_censoring(ce1, cens_at)

  dat <- rbind(ce0, ce1)

  if (sum(dat$status[dat$arm == 0]) == 0 ||
      sum(dat$status[dat$arm == 1]) == 0) return(NA_real_)

  sd <- try(survdiff(Surv(time, status) ~ arm, data = dat), silent = TRUE)
  if (inherits(sd, "try-error")) return(NA_real_)
  pchisq(sd$chisq, df = 1, lower.tail = FALSE)
}


# Internal null-coalesce
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Power study driver for the KM log-rank test
#'
#' Runs a Monte Carlo power study across one or more scenarios using the
#' Kaplan-Meier log-rank test as the test statistic. For each scenario,
#' \code{M} two-arm trials are simulated via \code{\link{simulate_one_trial}}
#' and the empirical rejection rate is reported alongside a Monte Carlo
#' standard error.
#'
#' @param scenarios A named list of scenario lists. Each scenario may
#'   override \code{mean_vec}, \code{covariance}, \code{alpha},
#'   \code{beta}, \code{gamma}, \code{cens_at}; must contain \code{beta}
#'   and a human-readable \code{label}.
#' @param defaults A list of default DGP parameters (see
#'   \code{\link{simulate_one_trial}}).
#' @param M Integer. Number of replications per scenario. Defaults to 1000.
#' @param n_per_arm Integer. Patients per arm. Defaults to 150.
#' @param alpha_level Numeric in \eqn{(0, 1)}. Rejection threshold.
#'   Defaults to 0.05.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per scenario and columns
#'   \code{scenario}, \code{rejection}, \code{mc_se}, \code{n_valid}.
#'
#' @export
run_power_study <- function(scenarios, defaults, M = 1000, n_per_arm = 150,
                            alpha_level = 0.05, verbose = TRUE) {

  results <- data.frame(
    scenario  = character(),
    rejection = numeric(),
    mc_se     = numeric(),
    n_valid   = integer(),
    stringsAsFactors = FALSE
  )

  for (s_name in names(scenarios)) {
    s <- scenarios[[s_name]]
    if (verbose) cat(sprintf("\n--- %s ---\n", s$label))

    pvals <- numeric(M)
    for (m in seq_len(M)) {
      if (verbose && m %% 200 == 0) cat("  rep", m, "/", M, "\n")
      pvals[m] <- simulate_one_trial(s, n_per_arm, defaults)
    }

    valid     <- !is.na(pvals)
    n_valid   <- sum(valid)
    rejection <- mean(pvals[valid] < alpha_level)
    mc_se     <- sqrt(rejection * (1 - rejection) / n_valid)

    results <- rbind(results, data.frame(
      scenario  = s$label,
      rejection = rejection,
      mc_se     = mc_se,
      n_valid   = n_valid,
      stringsAsFactors = FALSE
    ))
  }

  results
}


#' Single replication of the copula-based two-arm Wald test
#'
#' Inner loop for \code{\link{run_copula_power}}: simulates a control and
#' treatment arm under a scenario, fits per-arm copula PFS curves with
#' bootstrap standard errors, and returns the difference Delta and its
#' standard error at each grid point.
#'
#' @param m Integer. Replication index, used to derive a per-replication
#'   seed.
#' @param s Scenario list (see \code{\link{simulate_one_trial}}).
#' @param n_per_arm Integer.
#' @param base_seed Integer.
#' @param B Integer. Number of bootstrap replicates per arm.
#' @param family Character copula family.
#' @param defaults Defaults list (see \code{\link{simulate_one_trial}}).
#' @param grid Numeric vector of evaluation visits.
#'
#' @return A list with elements \code{delta} and \code{se}, each numeric
#'   vectors of length \code{length(grid)}; both \code{NA_real_} if a
#'   per-arm copula fit failed.
#'
#' @export
one_copula_rep <- function(m, s, n_per_arm, base_seed, B, family,
                           defaults, grid) {

  set.seed(base_seed + m)

  mean_vec <- if (!is.null(s$mean_vec))   s$mean_vec   else defaults$mean
  cov_mat  <- if (!is.null(s$covariance)) s$covariance else defaults$covariance
  alpha    <- if (!is.null(s$alpha))      s$alpha      else defaults$alpha
  gamma    <- if (!is.null(s$gamma))      s$gamma      else defaults$gamma
  cens_at  <- if (!is.null(s$cens_at))    s$cens_at    else Inf
  beta     <- s$beta
  threshold <- defaults$threshold

  # Control arm
  t0 <- tumour_events(n_per_arm, mean_vec, cov_mat, threshold)
  l0 <- lesion_events(n_per_arm, t0$log_ratios_vs_baseline,
                      alpha = alpha, beta = beta, gamma = gamma,
                      treatment_arm = 0)
  fit0 <- copula_arm(t0$events, l0$events, grid, family, B, cens_at)

  # Treatment arm
  t1 <- tumour_events(n_per_arm, mean_vec, cov_mat, threshold)
  l1 <- lesion_events(n_per_arm, t1$log_ratios_vs_baseline,
                      alpha = alpha, beta = beta, gamma = gamma,
                      treatment_arm = 1)
  fit1 <- copula_arm(t1$events, l1$events, grid, family, B, cens_at)

  if (is.null(fit0) || is.null(fit1)) {
    return(list(delta = rep(NA_real_, length(grid)),
                se    = rep(NA_real_, length(grid))))
  }

  list(delta = fit1$pfs - fit0$pfs,
       se    = sqrt(fit1$se^2 + fit0$se^2))
}


#' Power study driver for the copula-based two-arm Wald test
#'
#' Runs a Monte Carlo power study across scenarios using the two-arm
#' Wald test on the difference Delta = S_treat - S_control at a chosen
#' visit, with bootstrap-based per-arm standard errors from the copula
#' model.
#'
#' @param scenarios Named list of scenario lists.
#' @param defaults Defaults list passed through to
#'   \code{\link{one_copula_rep}}.
#' @param M Integer. Number of replications. Defaults to 1000.
#' @param B Integer. Bootstrap replicates per copula fit. Defaults to 200.
#' @param n_per_arm Integer. Patients per arm. Defaults to 150.
#' @param alpha_level Numeric. Rejection threshold. Defaults to 0.05.
#' @param test_visit Integer. Visit at which to perform the Wald test.
#'   Defaults to 5.
#' @param family Character copula family. Defaults to \code{"gumbel"}.
#' @param grid Numeric vector of evaluation visits. Defaults to \code{1:5}.
#' @param base_seed Integer.
#' @param n_cores Integer. Number of parallel workers. \code{NULL} uses
#'   \code{detectCores() - 1}.
#' @param verbose Logical.
#'
#' @return A \code{data.frame} with columns \code{scenario},
#'   \code{rejection}, \code{mc_se}, \code{n_valid}.
#'
#' @details For scenarios with administrative censoring before
#'   \code{test_visit}, the test is performed at \code{min(test_visit, cens_at)}.
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ clusterExport
#' @importFrom pbapply pblapply
#' @importFrom stats qnorm
#' @export
run_copula_power <- function(scenarios, defaults,
                             M = 1000, B = 200,
                             n_per_arm = 150, alpha_level = 0.05,
                             test_visit = 5, family = "gumbel",
                             grid = 1:5, base_seed = 1, n_cores = NULL,
                             verbose = TRUE) {

  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)

  results <- data.frame(
    scenario  = character(),
    rejection = numeric(),
    mc_se     = numeric(),
    n_valid   = integer(),
    stringsAsFactors = FALSE
  )

  for (s_name in names(scenarios)) {
    s <- scenarios[[s_name]]
    if (verbose) cat(sprintf("\n--- %s ---\n", s$label))

    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
      clusterEvalQ(cl, {
        library(MASS); library(mvtnorm); library(survival); library(copula)
        library(tumourSim)
      })
      clusterExport(cl,
                    varlist = c("s", "n_per_arm", "base_seed", "B",
                                "family", "defaults", "grid"),
                    envir = environment())
    }

    out <- pblapply(seq_len(M),
                    function(m) one_copula_rep(m, s, n_per_arm,
                                               base_seed, B, family,
                                               defaults, grid),
                    cl = cl)

    if (!is.null(cl)) { stopCluster(cl); cl <- NULL }

    delta_mat <- do.call(rbind, lapply(out, `[[`, "delta"))
    se_mat    <- do.call(rbind, lapply(out, `[[`, "se"))

    eff_visit <- if (!is.null(s$cens_at) && is.finite(s$cens_at)) {
      min(test_visit, s$cens_at)
    } else test_visit

    z <- qnorm(1 - alpha_level / 2)
    Z_t <- delta_mat[, eff_visit] / se_mat[, eff_visit]
    reject <- abs(Z_t) > z

    valid     <- !is.na(reject)
    n_valid   <- sum(valid)
    rejection <- mean(reject[valid])
    mc_se     <- sqrt(rejection * (1 - rejection) / n_valid)

    results <- rbind(results, data.frame(
      scenario  = s$label,
      rejection = rejection,
      mc_se     = mc_se,
      n_valid   = n_valid,
      stringsAsFactors = FALSE
    ))
  }

  results
}


#' Three-method Monte Carlo study for a single scenario
#'
#' Runs \code{M} two-arm replications under one scenario and applies all
#' three tests (KM log-rank, AugBin pooled-hazard with delta variance,
#' copula L1 permutation) to each replicate. Returns rejection rates and
#' the raw matrix of p-values.
#'
#' @param M Integer. Number of replications.
#' @param B Integer. Bootstrap/permutation replicates for the copula
#'   test.
#' @param n_per_arm Integer. Patients per arm.
#' @param mean_ctrl,mean_trt Numeric vectors of tumour log-ratio means.
#' @param cov_mat Numeric \code{T} x \code{T} covariance matrix.
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#' @param alpha_lesion,beta_lesion,gamma_lesion Lesion-model coefficients.
#' @param gamma_lesion_ctrl,gamma_lesion_trt Optional per-arm overrides for
#'   the tumour-to-lesion coefficient.
#' @param tumour_dist,tumour_df,lesion_link,lesion_nonlinear Optional
#'   model-misspecification arguments passed to
#'   \code{\link{simulate_trial}}.
#' @param copula_family Character copula family. Defaults to \code{"gumbel"}.
#' @param base_seed Integer.
#' @param n_cores Integer or \code{NULL}.
#' @param alpha_level Numeric. Rejection threshold. Defaults to 0.05.
#' @param verbose Logical.
#' @param augbin_test Character. Either \code{"pooled"} (default; uses
#'   \code{\link{augbin_test_pooled}}) or \code{"delta"} (uses
#'   \code{\link{augbin_test_delta}}).
#' @param copula_test_type Character. Either \code{"L1"} (default; uses
#'   \code{\link{copula_L1_permutation_test}}) or \code{"max"} (uses
#'   \code{\link{copula_test}}).
#'
#' @return A list with elements \code{power_km}, \code{power_augbin},
#'   \code{power_copula}, \code{n_failed}, and \code{pvalues}.
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ clusterExport
#' @importFrom pbapply pblapply
#' @export
run_scenario <- function(M, B, n_per_arm, mean_ctrl, mean_trt, cov_mat,
                         threshold = 1.2,
                         alpha_lesion = -2.5, beta_lesion = 0,
                         gamma_lesion = 0.2,
                         gamma_lesion_ctrl = NULL,
                         gamma_lesion_trt  = NULL,
                         tumour_dist = "mvnorm", tumour_df = 4,
                         lesion_link = "logit",
                         lesion_nonlinear = FALSE,
                         copula_family = "gumbel",
                         base_seed = 1, n_cores = NULL,
                         alpha_level = 0.05,
                         verbose = TRUE,
                         augbin_test = c("pooled", "delta"),
                         copula_test_type = c("L1", "max")) {

  augbin_test      <- match.arg(augbin_test)
  copula_test_type <- match.arg(copula_test_type)

  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  T_visits <- length(mean_ctrl)

  one_rep <- function(m) {
    set.seed(base_seed + m)
    trial <- simulate_trial(n_per_arm, mean_ctrl, mean_trt, cov_mat, threshold,
                            alpha_lesion, beta_lesion, gamma_lesion,
                            gamma_lesion_ctrl = gamma_lesion_ctrl,
                            gamma_lesion_trt  = gamma_lesion_trt,
                            tumour_dist = tumour_dist, tumour_df = tumour_df,
                            lesion_link = lesion_link,
                            lesion_nonlinear = lesion_nonlinear)
    p_km <- tryCatch(km_logrank_pvalue(trial), error = function(e) NA_real_)
    p_aug <- if (augbin_test == "pooled") {
      tryCatch(augbin_test_pooled(trial, threshold), error = function(e) NA_real_)
    } else {
      tryCatch(augbin_test_delta(trial, threshold), error = function(e) NA_real_)
    }
    p_cop <- if (copula_test_type == "L1") {
      tryCatch(copula_L1_permutation_test(trial, T_visits, B, copula_family),
               error = function(e) NA_real_)
    } else {
      tryCatch(copula_test(trial, B, grid = 1:T_visits, family = copula_family),
               error = function(e) NA_real_)
    }
    c(km = p_km, aug = p_aug, cop = p_cop)
  }

  if (verbose) {
    cat(sprintf("Scenario: M=%d, B=%d (copula only), n_per_arm=%d, %d cores\n",
                M, B, n_per_arm, n_cores))
    t0 <- Sys.time()
  }

  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
      library(MASS); library(mvtnorm); library(dplyr)
      library(survival); library(copula); library(tumourSim)
    })
    clusterExport(cl,
                  varlist = c("n_per_arm", "mean_ctrl", "mean_trt", "cov_mat",
                              "threshold", "alpha_lesion", "beta_lesion",
                              "gamma_lesion", "gamma_lesion_ctrl",
                              "gamma_lesion_trt",
                              "tumour_dist", "tumour_df",
                              "lesion_link", "lesion_nonlinear",
                              "copula_family", "B",
                              "base_seed", "T_visits",
                              "augbin_test", "copula_test_type"),
                  envir = environment())
    out <- pblapply(seq_len(M), one_rep, cl = cl)
  } else {
    out <- pblapply(seq_len(M), one_rep)
  }

  if (verbose) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  done in %.1f min (%.2f s/rep)\n", el / 60, el / M))
  }

  pmat <- do.call(rbind, out)
  list(
    power_km     = mean(pmat[, "km"]  < alpha_level, na.rm = TRUE),
    power_augbin = mean(pmat[, "aug"] < alpha_level, na.rm = TRUE),
    power_copula = mean(pmat[, "cop"] < alpha_level, na.rm = TRUE),
    n_failed     = colSums(is.na(pmat)),
    pvalues      = pmat
  )
}


#' Three-method power table across multiple scenarios with checkpointing
#'
#' Loops over a list of two-arm scenarios, calls
#' \code{\link{run_scenario}} for each, and assembles a tidy power table.
#' Intermediate results are written to disk as \code{.rds} files after
#' each scenario so long runs can be inspected or resumed.
#'
#' @param scenarios A list of scenario specifications. Each entry must
#'   contain \code{name}, \code{mean_trt}, \code{ctrl_shift} (added to
#'   \code{mean_trt} to obtain \code{mean_ctrl}), and \code{beta}; may
#'   contain optional fields \code{gamma_lesion_ctrl},
#'   \code{gamma_lesion_trt}, \code{tumour_dist}, \code{tumour_df},
#'   \code{lesion_link}, \code{lesion_nonlinear}.
#' @param cov_mat Numeric \code{T} x \code{T} covariance matrix shared by
#'   all scenarios.
#' @param M Integer. Replications per scenario. Defaults to 1000.
#' @param B Integer. Bootstrap/permutation replicates for the copula
#'   test. Defaults to 200.
#' @param n_per_arm Integer. Patients per arm. Defaults to 200.
#' @param n_cores Integer or \code{NULL}.
#' @param base_seed Integer.
#' @param checkpoint_dir Character. Directory for the checkpoint
#'   \code{.rds} files. Defaults to \code{"."}.
#' @param ... Further arguments passed through to
#'   \code{\link{run_scenario}}.
#'
#' @return A list with elements \code{table} (a \code{data.frame} with
#'   columns \code{scenario}, \code{beta}, \code{power_km},
#'   \code{power_augbin}, \code{power_copula}) and \code{raw} (the full
#'   per-scenario result lists).
#'
#' @export
run_all_scenarios <- function(scenarios, cov_mat,
                              M = 1000, B = 200, n_per_arm = 200,
                              n_cores = NULL, base_seed = 1,
                              checkpoint_dir = ".", ...) {

  results <- vector("list", length(scenarios))
  for (s in seq_along(scenarios)) {
    sc <- scenarios[[s]]
    mean_ctrl <- sc$mean_trt + (sc$ctrl_shift %||% 0)
    cat(sprintf("\n=== [%d/%d] %s ===\n", s, length(scenarios), sc$name))
    res <- run_scenario(
      M = M, B = B, n_per_arm = n_per_arm,
      mean_ctrl = mean_ctrl, mean_trt = sc$mean_trt,
      cov_mat = cov_mat,
      beta_lesion = sc$beta,
      gamma_lesion_ctrl = sc$gamma_lesion_ctrl %||% NULL,
      gamma_lesion_trt  = sc$gamma_lesion_trt  %||% NULL,
      tumour_dist = sc$tumour_dist %||% "mvnorm",
      tumour_df   = sc$tumour_df   %||% 4,
      lesion_link = sc$lesion_link %||% "logit",
      lesion_nonlinear = sc$lesion_nonlinear %||% FALSE,
      base_seed = base_seed + 1000 * s,
      n_cores = n_cores,
      ...
    )
    results[[s]] <- list(
      scenario     = sc$name,
      beta         = sc$beta,
      mean_trt     = sc$mean_trt,
      power_km     = res$power_km,
      power_augbin = res$power_augbin,
      power_copula = res$power_copula,
      n_failed     = res$n_failed
    )
    saveRDS(results,
            file.path(checkpoint_dir,
                      sprintf("checkpoint_through_scenario_%d.rds", s)))
    cat(sprintf("  KM = %.3f   AugBin = %.3f   Copula = %.3f\n",
                res$power_km, res$power_augbin, res$power_copula))
  }

  tab <- do.call(rbind, lapply(results, function(r) {
    data.frame(scenario     = r$scenario,
               beta         = r$beta,
               power_km     = r$power_km,
               power_augbin = r$power_augbin,
               power_copula = r$power_copula)
  }))
  list(table = tab, raw = results)
}


#' Monte Carlo coverage and width comparison: KM vs copula PFS
#'
#' Single-arm Monte Carlo study comparing the direct Kaplan-Meier
#' estimator of PFS on the composite endpoint against the copula-based
#' estimator from \code{\link{pfs_copula_boot}} across one or more
#' copula families.
#'
#' @param M Integer. Number of replications. Defaults to 1000.
#' @param B Integer. Bootstrap replicates per copula fit. Defaults to 200.
#' @param n_patients Integer. Patients per simulated trial. Defaults to 150.
#' @param mean,covariance,threshold Tumour-model parameters.
#' @param alpha_lesion,beta_lesion,gamma_lesion Lesion-model parameters.
#' @param treatment_arm Numeric (typically 0/1).
#' @param grid Numeric vector of evaluation visits. Defaults to \code{1:5}.
#' @param families Character vector of copula families. Defaults to
#'   \code{c("clayton", "frank", "gumbel")}.
#' @param true_S Optional numeric vector of true PFS values on
#'   \code{grid}. If \code{NULL}, computed by \code{\link{true_pfs}}.
#' @param verbose Logical. Defaults to \code{TRUE}.
#'
#' @return A list with elements \code{summary} (long-format
#'   \code{data.frame}) and \code{raw} (per-replication matrices and
#'   per-family lists).
#'
#' @importFrom survival survfit Surv
#' @export
simulate_metrics_all <- function(M = 1000, B = 200, n_patients = 150,
                                 mean, covariance, threshold,
                                 alpha_lesion = -2.5, beta_lesion = 0,
                                 gamma_lesion = 0.2, treatment_arm = 0,
                                 grid = 1:5,
                                 families = c("clayton", "frank", "gumbel"),
                                 true_S = NULL, verbose = TRUE) {

  if (is.null(true_S)) {
    if (verbose) cat("Computing true PFS via large simulation...\n")
    true_S <- true_pfs(grid, mean = mean, covariance = covariance,
                       threshold = threshold,
                       alpha = alpha_lesion, beta = beta_lesion,
                       gamma = gamma_lesion, treatment_arm = treatment_arm)
  }

  G <- length(grid)

  pfs_km   <- matrix(NA_real_, M, G)
  width_km <- matrix(NA_real_, M, G)
  cov_km   <- matrix(NA,       M, G)

  pfs_cop   <- stats::setNames(lapply(families, function(f) matrix(NA_real_, M, G)), families)
  width_cop <- stats::setNames(lapply(families, function(f) matrix(NA_real_, M, G)), families)
  cov_cop   <- stats::setNames(lapply(families, function(f) matrix(NA,       M, G)), families)
  theta_cop <- stats::setNames(lapply(families, function(f) numeric(M)),             families)
  aic_cop   <- stats::setNames(lapply(families, function(f) numeric(M)),             families)

  for (m in seq_len(M)) {
    if (verbose && m %% 50 == 0) cat("Replication", m, "/", M, "\n")

    ts <- tumour_events(n_patients, mean, covariance, threshold)
    ls <- lesion_events(n_patients, ts$log_ratios_vs_baseline,
                        alpha = alpha_lesion, beta = beta_lesion,
                        gamma = gamma_lesion, treatment_arm = treatment_arm)

    tt <- find_event_time(ts$events)
    ll <- find_event_time(ls$events)
    ce <- combined_event(ls$events, ts$events)

    kmf  <- survfit(Surv(time, status) ~ 1, data = ce)
    km_s <- summary(kmf, times = grid, extend = TRUE)
    pfs_km[m, ]   <- km_s$surv
    width_km[m, ] <- km_s$upper - km_s$lower
    cov_km[m, ]   <- (km_s$lower <= true_S) & (true_S <= km_s$upper)

    for (fam in families) {
      pt <- try(pfs_copula(ll$time, ll$status, tt$time, tt$status,
                           grid = grid, family = fam), silent = TRUE)
      if (!inherits(pt, "try-error")) {
        theta_cop[[fam]][m] <- pt$theta
        aic_cop[[fam]][m]   <- pt$aic
      }
      rc <- try(pfs_copula_boot(ll$time, ll$status, tt$time, tt$status,
                                grid = grid, family = fam, B = B),
                silent = TRUE)
      if (!inherits(rc, "try-error")) {
        pfs_cop[[fam]][m, ]   <- rc$pfs
        width_cop[[fam]][m, ] <- rc$upper - rc$lower
        cov_cop[[fam]][m, ]   <- (rc$lower <= true_S) & (true_S <= rc$upper)
      }
    }
  }

  km_summary <- data.frame(
    method = "km", family = NA_character_, time = grid, true_S = true_S,
    pfs      = colMeans(pfs_km,   na.rm = TRUE),
    width    = colMeans(width_km, na.rm = TRUE),
    coverage = colMeans(cov_km,   na.rm = TRUE),
    theta    = NA_real_, aic = NA_real_)

  cop_summary <- do.call(rbind, lapply(families, function(fam) {
    data.frame(
      method = "copula", family = fam, time = grid, true_S = true_S,
      pfs      = colMeans(pfs_cop[[fam]],   na.rm = TRUE),
      width    = colMeans(width_cop[[fam]], na.rm = TRUE),
      coverage = colMeans(cov_cop[[fam]],   na.rm = TRUE),
      theta    = mean(theta_cop[[fam]],     na.rm = TRUE),
      aic      = mean(aic_cop[[fam]],       na.rm = TRUE))
  }))

  list(summary = rbind(km_summary, cop_summary),
       raw = list(pfs_km = pfs_km, width_km = width_km, cov_km = cov_km,
                  pfs_cop = pfs_cop, width_cop = width_cop, cov_cop = cov_cop,
                  theta_cop = theta_cop, aic_cop = aic_cop))
}


#' Parallel Monte Carlo evaluator for the AugBin PFS estimator
#'
#' Single-arm Monte Carlo evaluator for the AugBin estimator. Each
#' replication is seeded for reproducibility regardless of platform or
#' core count. Returns mean point estimate, bias, coverage rate, mean CI
#' width, and Monte Carlo standard errors for each.
#'
#' @param M Integer. Number of replications. Defaults to 1000.
#' @param B Integer. Bootstrap replicates per replication. Defaults to 200.
#' @param n_patients Integer.
#' @param mean_vec Numeric vector of length \code{T} (tumour mean).
#' @param cov_mat Numeric \code{T} x \code{T} covariance matrix.
#' @param S_true Numeric vector of true PFS values at visits \code{1:T},
#'   used to evaluate coverage and bias.
#' @param threshold Numeric. Tumour progression threshold. Defaults to 1.2.
#' @param alpha_lesion,beta_lesion,gamma_lesion,treatment_arm Lesion-model
#'   parameters.
#' @param base_seed Integer.
#' @param n_cores Integer or \code{NULL}.
#' @param verbose Logical.
#'
#' @return A list with elements \code{estimate}, \code{bias},
#'   \code{coverage_rate}, \code{avg_ci_width},
#'   \code{mc_se_estimate}, \code{mc_se_coverage}, \code{mc_se_width},
#'   \code{n_failed}.
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ clusterExport
#' @importFrom pbapply pblapply
#' @export
evaluate_augbin_parallel <- function(M = 1000, B = 200, n_patients, mean_vec, cov_mat,
                                     S_true, threshold = 1.2,
                                     alpha_lesion = -2.5, beta_lesion = 0,
                                     gamma_lesion = 0.2, treatment_arm = 0,
                                     base_seed = 1, n_cores = NULL,
                                     verbose = TRUE) {

  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  T_visits  <- length(mean_vec)
  treatment <- rep(treatment_arm, n_patients)

  one_rep <- function(m) {
    set.seed(base_seed + m)

    tumour <- tumour_events(n_patients, mean_vec, cov_mat, threshold)
    lesion <- lesion_events(n_patients, tumour$log_ratios_vs_baseline,
                            alpha = alpha_lesion, beta = beta_lesion,
                            gamma = gamma_lesion,
                            treatment_arm = treatment_arm)

    pe <- tryCatch(
      augbin_estimate(tumour, lesion, treatment, threshold)$S_pfs,
      error = function(e) rep(NA_real_, T_visits)
    )

    bs <- tryCatch(
      augbin_bootstrap(tumour, lesion, treatment, B, threshold),
      error = function(e) list(lower = rep(NA_real_, T_visits),
                               upper = rep(NA_real_, T_visits))
    )

    list(point = pe,
         cover = (S_true >= bs$lower) & (S_true <= bs$upper),
         width = bs$upper - bs$lower)
  }

  if (verbose) {
    cat(sprintf("Running M = %d, B = %d on %d core(s)...\n", M, B, n_cores))
    t0 <- Sys.time()
  }

  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
      library(MASS); library(mvtnorm); library(dplyr); library(tumourSim)
    })
    clusterExport(cl,
                  varlist = c("n_patients", "mean_vec", "cov_mat", "S_true",
                              "threshold", "alpha_lesion", "beta_lesion",
                              "gamma_lesion", "treatment_arm", "treatment",
                              "B", "T_visits", "base_seed"),
                  envir = environment())
  } else {
    cl <- NULL
  }

  out <- pblapply(seq_len(M), one_rep, cl = cl)

  if (verbose) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("Done in %.1f min (%.2f s/replicate).\n", el / 60, el / M))
  }

  point_est <- do.call(rbind, lapply(out, `[[`, "point"))
  cover     <- do.call(rbind, lapply(out, `[[`, "cover"))
  width     <- do.call(rbind, lapply(out, `[[`, "width"))

  est_mean <- colMeans(point_est, na.rm = TRUE)

  list(
    estimate       = est_mean,
    bias           = est_mean - S_true,
    coverage_rate  = colMeans(cover, na.rm = TRUE),
    avg_ci_width   = colMeans(width, na.rm = TRUE),
    mc_se_estimate = apply(point_est, 2, sd, na.rm = TRUE) / sqrt(M),
    mc_se_coverage = apply(cover,     2, sd, na.rm = TRUE) / sqrt(M),
    mc_se_width    = apply(width,     2, sd, na.rm = TRUE) / sqrt(M),
    n_failed       = sum(apply(is.na(point_est), 1, any))
  )
}
