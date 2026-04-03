# tests/testthat/test-generate-continuous-data.R

# helpers ----
make_cov <- function(n) diag(n)  # valid identity covariance for any n

# happy path ----
test_that("returns a list with the two expected names", {
  out <- generate_continuous_data(3, 5, mean = rep(0, 3), covariance = make_cov(3))
  expect_named(out, c("log_tumour_size_ratio", "baseline_tumour_size"))
})

test_that("log_tumour_size_ratio has dimensions n_patients x n_times", {
  out <- generate_continuous_data(4, 10, mean = rep(0, 4), covariance = make_cov(4))
  expect_equal(dim(out$log_tumour_size_ratio), c(10, 4))
})

test_that("baseline_tumour_size has length n_patients", {
  out <- generate_continuous_data(3, 7, mean = rep(0, 3), covariance = make_cov(3))
  expect_length(out$baseline_tumour_size, 7)
})

test_that("baseline_tumour_size is bounded in [0, 1]", {
  out <- generate_continuous_data(3, 100, mean = rep(0, 3), covariance = make_cov(3))
  expect_true(all(out$baseline_tumour_size >= 0))
  expect_true(all(out$baseline_tumour_size <= 1))
})

test_that("log_tumour_size_ratio contains only finite numeric values", {
  out <- generate_continuous_data(3, 10, mean = rep(0, 3), covariance = make_cov(3))
  expect_true(is.numeric(out$log_tumour_size_ratio))
  expect_true(all(is.finite(out$log_tumour_size_ratio)))
})

test_that("non-identity covariance is accepted and produces correct output shape", {
  cov_mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  out <- generate_continuous_data(2, 5, mean = c(0, 0), covariance = cov_mat)
  expect_equal(dim(out$log_tumour_size_ratio), c(5, 2))
})

test_that("non-zero mean shifts log_tumour_size_ratio distribution", {
  set.seed(1)
  out_zero <- generate_continuous_data(1, 1000, mean = 0,   covariance = matrix(1))
  out_high <- generate_continuous_data(1, 1000, mean = 100, covariance = matrix(1))
  expect_gt(mean(out_high$log_tumour_size_ratio), mean(out_zero$log_tumour_size_ratio))
})

test_that("output is reproducible with set.seed()", {
  set.seed(42)
  out1 <- generate_continuous_data(3, 5, mean = rep(0, 3), covariance = make_cov(3))
  set.seed(42)
  out2 <- generate_continuous_data(3, 5, mean = rep(0, 3), covariance = make_cov(3))
  expect_equal(out1$log_tumour_size_ratio, out2$log_tumour_size_ratio)
  expect_equal(out1$baseline_tumour_size,  out2$baseline_tumour_size)
})

# edge cases ----
test_that("n_patients = 1 returns a 1-row matrix and length-1 vector", {
  out <- generate_continuous_data(3, 1, mean = rep(0, 3), covariance = make_cov(3))
  expect_equal(dim(out$log_tumour_size_ratio), c(1, 3))
  expect_length(out$baseline_tumour_size, 1)
})

test_that("n_times = 1 is accepted and returns correct dimensions", {
  out <- generate_continuous_data(1, 5, mean = 0, covariance = matrix(1))
  expect_equal(dim(out$log_tumour_size_ratio), c(5, 1))
})

# input validation ----
test_that("non-numeric n_times throws an error", {
  expect_error(
    generate_continuous_data("3", 5, rep(0, 3), make_cov(3)),
    "`n_times` must be an integer"
  )
})

test_that("fractional n_times throws an error", {
  expect_error(
    generate_continuous_data(2.5, 5, rep(0, 3), make_cov(3)),
    "`n_times` must be an integer"
  )
})

test_that("non-numeric n_patients throws an error", {
  expect_error(
    generate_continuous_data(3, "5", rep(0, 3), make_cov(3)),
    "`n_patients` must be an integer"
  )
})

test_that("fractional n_patients throws an error", {
  expect_error(
    generate_continuous_data(3, 1.5, rep(0, 3), make_cov(3)),
    "`n_patients` must be an integer"
  )
})

test_that("mean vector of wrong length throws an error", {
  expect_error(
    generate_continuous_data(3, 5, mean = c(0, 0), covariance = make_cov(3)),
    "Length of mean vector must equal the number of time points"
  )
})

test_that("non-vector mean throws an error", {
  expect_error(
    generate_continuous_data(3, 5, mean = matrix(0, 1, 3), covariance = make_cov(3)),
    "Length of mean vector must equal the number of time points"
  )
})

test_that("non-matrix covariance throws an error", {
  expect_error(
    generate_continuous_data(3, 5, mean = rep(0, 3), covariance = rep(1, 9)),
    "Covariance matrix must be a n_times square matrix"
  )
})

test_that("non-square covariance throws an error", {
  expect_error(
    generate_continuous_data(3, 5, mean = rep(0, 3), covariance = matrix(1, nrow = 3, ncol = 2)),
    "Covariance matrix must be a n_times square matrix"
  )
})

test_that("covariance matrix of wrong size throws an error", {
  expect_error(
    generate_continuous_data(3, 5, mean = rep(0, 3), covariance = make_cov(2)),
    "Covariance matrix must be a n_times square matrix"
  )
})
