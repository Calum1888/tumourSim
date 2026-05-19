library(tumourSim)

# parameters for single arm trial

MEAN5 = c(0.00, 0.036, 0.072, 0.108, 0.144)
COVARIANCE5 <- matrix(
  c(
    0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.45, 0.45, 0.45, 0.45,
    0.25, 0.45, 0.50, 0.50, 0.50,
    0.25, 0.45, 0.50, 0.75, 0.75,
    0.25, 0.45, 0.50, 0.75, 1.00
  ),
  nrow = 5,
  byrow = TRUE
)

MEAN7 <- c(-0.22, -0.29, -0.36, -0.36, -0.43, -0.51, -0.36)

Sigma_7 <- matrix(0.15, nrow = N_TIMES, ncol = N_TIMES)
diag(Sigma_7) <- 0.20
COVARIANCE7 <- Sigma_7

N_SAMPLES_SINGLE_ARM <- 150
THRESHOLD <- 1.2

# parameters for comparative trial

MEAN_COMP  <- c(-0.2, -0.4, -0.56,-0.6,-0.65)
COV_COMP   <-  matrix(c(
  0.05, 0.05, 0.05, 0.05, 0.05,
  0.05, 0.10, 0.10, 0.10, 0.10,
  0.05, 0.10, 0.14, 0.14, 0.14,
  0.05, 0.10, 0.14, 0.16, 0.16,
  0.05, 0.10, 0.14, 0.16, 0.18
), nrow = 5, byrow = TRUE)
ALPHA_COMP <- -2.5
GAMMA_COMP <- 0.2
CENS_COMP  <- Inf # no admin censoring

N_SAMPLES_COMPARATIVE_ARM <- 200

# run simulation for 5 time points
