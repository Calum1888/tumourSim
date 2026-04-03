test_that("returns a list with the four expected names", {
  out <- generate_coefficients(5, 3, 0.1, 0.2, 0.3, 1)
  expect_named(out, c("alpha_r", "beta_r", "gamma_r", "R_indicator"))
})

test_that("alpha, beta, gamma vectors have length n_times", {
  out <- generate_coefficients(n_times = 5, n_patients = 3, 0.1, 0.2, 0.3, 1)
  expect_length(out$alpha_r, 5)
  expect_length(out$beta_r,  5)
  expect_length(out$gamma_r, 5)
})

test_that("R_indicator has length n_patients when R is scalar", {
  out <- generate_coefficients(5, 4, 0.1, 0.2, 0.3, R = 1)
  expect_length(out$R_indicator, 4)
})

test_that("scalar R is broadcast to all patients", {
  out <- generate_coefficients(5, 3, 0.1, 0.2, 0.3, R = 1)
  expect_equal(out$R_indicator, rep(1, 3))
})

test_that("vector R is passed through unchanged", {
  R_in <- c(0, 1, 1)
  out  <- generate_coefficients(5, 3, 0.1, 0.2, 0.3, R = R_in)
  expect_equal(out$R_indicator, R_in)
})

test_that("coefficient vectors are constant-valued", {
  out <- generate_coefficients(4, 2, alpha = 0.5, beta = 1.5, gamma = -0.3, R = 0)
  expect_true(all(out$alpha_r == 0.5))
  expect_true(all(out$beta_r  == 1.5))
  expect_true(all(out$gamma_r == -0.3))
})

test_that("n_times = 1 and n_patients = 1 works", {
  out <- generate_coefficients(1, 1, 0.1, 0.2, 0.3, R = 0)
  expect_length(out$alpha_r,     1)
  expect_length(out$R_indicator, 1)
})

test_that("negative and zero coefficient values are accepted", {
  expect_no_error(generate_coefficients(3, 2, alpha = -1, beta = 0, gamma = -0.5, R = 1))
})

test_that("R can be a vector of zeros and ones", {
  R_in <- c(0, 1, 0, 1)
  out  <- generate_coefficients(3, 4, 0.1, 0.2, 0.3, R = R_in)
  expect_equal(out$R_indicator, R_in)
})

test_that("non-numeric n_times throws an error", {
  expect_error(generate_coefficients("5", 3, 0.1, 0.2, 0.3, 1), "`n_times` must be an integer")
})

test_that("fractional n_times throws an error", {
  expect_error(generate_coefficients(2.5, 3, 0.1, 0.2, 0.3, 1), "`n_times` must be an integer")
})

test_that("non-numeric n_patients throws an error", {
  expect_error(generate_coefficients(5, "3", 0.1, 0.2, 0.3, 1), "`n_patients` must be an integer")
})

test_that("fractional n_patients throws an error", {
  expect_error(generate_coefficients(5, 1.5, 0.1, 0.2, 0.3, 1), "`n_patients` must be an integer")
})

test_that("non-numeric alpha throws an error", {
  expect_error(generate_coefficients(5, 3, "a", 0.2, 0.3, 1), "`alpha` must be numeric")
})

test_that("non-numeric beta throws an error", {
  expect_error(generate_coefficients(5, 3, 0.1, "b", 0.3, 1), "`beta` must be numeric")
})

test_that("non-numeric gamma throws an error", {
  expect_error(generate_coefficients(5, 3, 0.1, 0.2, "c", 1), "`gamma` must be numeric")
})

test_that("R of wrong length throws an error", {
  expect_error(generate_coefficients(5, 3, 0.1, 0.2, 0.3, R = c(1, 0)), "must be length 1 or length n_patients")
})
