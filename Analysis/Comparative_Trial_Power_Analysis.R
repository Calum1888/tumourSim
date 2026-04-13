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


betas <- c(-0.5)

power_curve <- lapply(betas, function(b) {
  res <- power_copula_pfs(
    n_times       = n_times,
    n_patients    = 100,
    n_iterations  = 100,   # increase for stable estimates
    mean          = mean_vec,
    covariance    = covariance,
    alpha_coef    = -1.5,
    beta          = b,
    gamma         = 0.0,
    copula_family = 5,
    B             = 50,
    seed          = 42
  )
  res$beta <- b
  res
})

power_curve_df <- do.call(rbind, power_curve)


# Define your scenarios in a named list
scenarios <- list(
  "Early_Diff" = c(-0.40, -0.60, -0.61, -0.65, -0.70),
  "Increasing" = c(-0.30, -0.50, -0.66, -0.70, -0.75),
  "Crossing"   = c(-0.50, -0.70, -0.63, -0.67, -0.72)
)

# Fix beta at -0.8 (or your choice)
fixed_beta <- -0.5

scenario_results <- lapply(names(scenarios), function(s_name) {
  res <- power_copula_pfs(
    n_times       = 5,
    n_patients    = 1000,
    n_iterations  = 10,
    mean          = scenarios[[s_name]], # Pass the scenario vector
    covariance    = covariance,
    alpha_coef    = -1.5,
    beta          = fixed_beta,          # Keep beta fixed
    gamma         = 0.2,
    copula_family = 5,
    B             = 100,
    seed          = 42
  )
  res$scenario <- s_name
  return(res)
})

final_df <- do.call(rbind, scenario_results)
