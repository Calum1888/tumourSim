library(tumourSim)
library(future)

plan(multisession)

mean_vec <- c(-0.2, -0.4, -0.56, -0.6, -0.65)

n_times    <- 5
covariance <- matrix(c(
  0.05, 0.05, 0.05, 0.05, 0.05,
  0.05, 0.10, 0.10, 0.10, 0.10,
  0.05, 0.10, 0.14, 0.14, 0.14,
  0.05, 0.10, 0.14, 0.16, 0.16,
  0.05, 0.10, 0.14, 0.16, 0.18
), nrow = 5, byrow = TRUE)


betas <- c(0.0, -0.2, -0.5, -0.8)

power_curve <- lapply(betas, function(b) {
  res <- power_copula_pfs(
    n_times       = n_times,
    n_patients    = 100,
    n_iterations  = 1000,   # increase for stable estimates
    mean          = mean_vec,
    covariance    = covariance,
    alpha_coef    = -1.5,
    beta          = b,
    gamma         = 0.2,
    copula_family = 5,
    B             = 200,
    seed          = 42
  )
  res$beta <- b
  res
})

log_rank_power <- lapply(betas, function(b){
  res_log <- power_logrank_pfs(
    n_times = n_times,
    n_patients = 100,
    n_iterations = 100,
    mean = mean_vec,
    covariance = covariance,
    alpha_coef = -1.5,
    beta = b,
    gamma = 0.2,
    alpha_level = 0.05,
    threshold = 1.2,
    seed = 42
  )
})

power_curve_df <- do.call(rbind, power_curve)

log_power_df <- do.call(rbind, log_rank_power)

plot_power_curve(power_df = power_curve_df,
                 title = "Power Curve for Different Treatment Effects")

# Define your scenarios in a named list
scenarios <- list(
  "Early_Diff" = c(-0.40, -0.60, -0.61, -0.65, -0.70),
  "Increasing" = c(-0.30, -0.50, -0.66, -0.70, -0.75),
  "Crossing"   = c(-0.50, -0.70, -0.63, -0.67, -0.72)
)


fixed_beta <- -0.5

scenario_results <- lapply(names(scenarios), function(s_name) {
  res <- power_copula_pfs(
    n_times       = 5,
    n_patients    = 100,
    n_iterations  = 1000,
    mean          = scenarios[[s_name]], # Pass the scenario vector
    covariance    = covariance,
    alpha_coef    = -1.5,
    beta          = fixed_beta,          # Keep beta fixed
    gamma         = 0.2,
    copula_family = 5,
    B             = 200,
    seed          = 42
  )
  res$scenario <- s_name
  return(res)
})

final_df <- do.call(rbind, scenario_results)

write.csv(power_curve_df, file = "Power_Curve_DF")
write.csv(final_df, file = "cross_diff_increase_scenarios")
