library(tumourSim)

# ---- Simulation parameters ---------------------------------------------------

N_TIMES      <- 5
N_PATIENTS   <- 150
N_ITERATIONS <- 100
B_BOOTSTRAP  <- 50
THRESHOLD    <- 1.2
SEED         <- 42

# Tumour size log-ratio distribution
MEAN_VEC <- mu_tumour_5 <- c(0, 0.036, 0.072, 0.108, 0.144)
COV_MAT  <- matrix(c(0.25, 0.25, 0.25, 0.25, 0.25,
                     0.25, 0.45, 0.45, 0.45, 0.45,
                     0.25, 0.45, 0.50, 0.50, 0.50,
                     0.25, 0.45, 0.50, 0.75, 0.75,
                     0.25, 0.45, 0.50, 0.75, 1.00),
                   nrow=N_TIMES, ncol=N_TIMES)

# Logistic new-lesion model parameters
ALPHA <- -1.5   # baseline intercept
BETA  <- 0   # treatment effect
GAMMA <-  0.2   # tumour size effect
R     <-  0     # control arm (no treatment)

# Copula families to evaluate
copula_families <- list(
  Clayton = 3,
  Gumbel  = 4,
  Frank   = 5
)

# ---- True PFS (large cohort) -------------------------------------------------

cat("Estimating true PFS from 100,000-patient cohort...\n")
set.seed(SEED)
true_pfs <- get_true_rates(
  n_times    = N_TIMES,
  mean       = MEAN_VEC,
  covariance = COV_MAT,
  alpha      = ALPHA,
  beta       = BETA,
  gamma      = GAMMA,
  R          = R,
  threshold  = THRESHOLD
)$surv

# ---- KM estimate -------------------------------------------------------------

cat("Running Kaplan-Meier simulations (", N_ITERATIONS, "iterations)...\n")
set.seed(SEED)
km_results <- run_iterations(
  n_times      = N_TIMES,
  n_patients   = N_PATIENTS,
  n_iterations = N_ITERATIONS,
  mean         = MEAN_VEC,
  covariance   = COV_MAT,
  alpha        = ALPHA,
  beta         = BETA,
  gamma        = GAMMA,
  R            = R,
  threshold    = THRESHOLD
)

cat("\nKaplan-Meier results:\n")
print(km_results)

# ---- Copula estimates --------------------------------------------------------

copula_results <- list()

for (family_name in names(copula_families)) {

  family_id <- copula_families[[family_name]]
  cat("\nRunning", family_name, "copula simulations (", N_ITERATIONS,
      "iterations,", B_BOOTSTRAP, "bootstrap samples each)...\n")

  set.seed(SEED)
  copula_results[[family_name]] <- run_copula_iterations(
    n_times      = N_TIMES,
    n_patients   = N_PATIENTS,
    n_iterations = N_ITERATIONS,
    mean         = MEAN_VEC,
    covariance   = COV_MAT,
    alpha        = ALPHA,
    beta         = BETA,
    gamma        = GAMMA,
    R            = R,
    copula_family = family_id,
    B             = B_BOOTSTRAP,
    threshold     = THRESHOLD
  )

  cat(family_name, "copula results:\n")
  print(copula_results[[family_name]])
}

# ---- Plots -------------------------------------------------------------------

cat("\nGenerating plots...\n")

for (family_name in names(copula_families)) {

  p <- plot_pfs_estimates(
    km_results     = km_results,
    copula_results = copula_results[[family_name]],
    true_pfs       = true_pfs,
    title          = paste("PFS Estimation --", family_name, "Copula")
  )

  filename <- paste0("pfs_plot_", tolower(family_name), ".png")
  ggplot2::ggsave(filename, plot = p, width = 12, height = 5, dpi = 300)
  cat("Saved:", filename, "\n")
}

cat("\nDone.\n")
