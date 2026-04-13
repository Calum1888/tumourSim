# Copula-Based Progression-Free Survival Estimation

An R package implementing a copula-based approach to estimate progression-free survival (PFS) by jointly modelling tumour progression and new lesion appearance, following the RECIST framework.

------------------------------------------------------------------------

## Background

In oncology trials, progression-free survival is defined as the time to the first of two competing events:

-   **Tumour progression** — tumour size exceeds 20% of the subject's rolling minimum (RECIST 1.1 criterion)
-   **New lesion appearance** — a binary indicator switches from 0 to 1

Traditional Kaplan–Meier methods treat PFS as a single combined endpoint. This package instead models the two endpoints separately with marginal Kaplan–Meier estimators and then reconstructs the joint PFS curve via the **survival-copula identity**:

$$S_{\text{PFS}}(t) = S_D(t) + S_Y(t) - 1 + C(1 - S_D(t),\, 1 - S_Y(t);\, \hat{\theta})$$

where $C$ is a bivariate copula fitted to pseudo-observations derived from the marginal event times.

------------------------------------------------------------------------

## Installation

``` r
# Install dependencies
install.packages(c("survival", "VineCopula", "MASS", "ggplot2", "rlang", "parallel", "scales", "patchwork"))

# Install from source
devtools::install_github("Calum1888/tumourSim")
```

------------------------------------------------------------------------

## Data Generation

Two functions simulate clinical trial data under the RECIST model.

### `generate_continuous_data()`

Simulates log tumour size ratios from a multivariate normal distribution and baseline tumour sizes from Uniform(0, 1).

``` r
n_times   <- 5
n_patients <- 150

mu <- c(0, 0.036, 0.072, 0.108, 0.144)

Sigma <- matrix(c(
  0.25, 0.25, 0.25, 0.25, 0.25,
  0.25, 0.45, 0.45, 0.45, 0.45,
  0.25, 0.45, 0.50, 0.50, 0.50,
  0.25, 0.45, 0.50, 0.75, 0.75,
  0.25, 0.45, 0.50, 0.75, 1.00
), nrow = 5)

data <- generate_continuous_data(n_times, n_patients, mean = mu, covariance = Sigma)
# Returns: $log_tumour_size_ratio [n_patients x n_times matrix]
#          $baseline_tumour_size  [n_patients vector]
```

### `generate_coefficients()`

Generates the logistic regression coefficients used to model the probability of a new lesion appearing at each time point.

``` r
coeffs <- generate_coefficients(
  n_times    = 5,
  n_patients = 150,
  alpha      = -2.5,   # baseline log-odds
  beta       = 0.5,    # treatment arm effect (set 0 for single-arm)
  gamma      = 0.3,    # tumour size effect
  R          = 0       # arm indicator: 0 = control, 1 = treatment
)
```

------------------------------------------------------------------------

## Core Functions

### `lesion_event()`

Derives binary event times from a lesion indicator matrix. Each subject's event time is the first column index where the indicator equals 1; subjects with no event are censored at the final time point.

``` r
lesion_data <- matrix(c(0, 0, 1,
                         0, 0, 0), nrow = 2, byrow = TRUE)
lesion_event(lesion_data)
#   time status
# 1    3      1
# 2    3      0
```

### `tumour_event()`

Derives progression event times from continuous tumour size measurements using a rolling-minimum threshold rule. An event is recorded when:

$$z_{it} > 1.2 \times \min(z_{i0}, \ldots, z_{i,t-1})$$

``` r
tumour_event(tumour_size_data, threshold = 1.2)
```

### `copula_pfs()`

The main estimator. Fits a bivariate copula to pseudo-observations from the marginal KM CDFs, then evaluates the PFS curve at each time point using the survival-copula identity.

``` r
# Gaussian copula (family = 1)
pfs <- copula_pfs(
  lesion_events  = lesion_event(lesion_data),
  tumour_events  = tumour_event(tumour_size_data),
  n_times        = 5,
  copula_family  = 1
)
```

Copula family codes follow the `VineCopula` package convention. Passing a vector of families triggers AIC-based selection via `BiCopSelect()`. Recommended families for survival data:

| Family code | Copula   | Tail dependence         |
|-------------|----------|-------------------------|
| 1           | Gaussian | None                    |
| 2           | t        | Upper and lower         |
| 3           | Clayton  | Lower only              |
| 4           | Gumbel   | Upper only              |
| 5           | Frank    | None (symmetric centre) |

### `copula_margin_estimation()`

Returns pointwise Kaplan–Meier survival estimates for each marginal endpoint at times `1:n_times`. Used internally by `copula_pfs()` but exposed for diagnostics.

------------------------------------------------------------------------

## Bootstrap Inference

`bootstrap_copula_pfs()` draws `B` subject-level bootstrap resamples and computes pointwise percentile confidence intervals, coverage against a known true PFS curve, and average CI width.

``` r

result <- bootstrap_copula_pfs(
  lesion_data      = lesion_data,
  tumour_size_data = tumour_size_data,
  n_times          = 5,
  copula_family    = c(1, 3, 4, 5),  # select best-fitting copula
  B                = 500,
  alpha            = 0.05,
  true_pfs         = true_pfs,
  threshold        = 1.2,
  seed             = 42
)

result$ci_lower      # lower CI bound at each time point
result$ci_upper      # upper CI bound at each time point
result$coverage      # 0/1 indicator per time point
result$mean_coverage # proportion of time points where CI contains true value
result$mean_ci_width # average CI width across time points
result$boot_curves   # B x n_times matrix of bootstrap PFS curves
```

Failed bootstrap iterations (e.g. degenerate resamples causing copula fitting to fail) are silently dropped with a warning reporting how many were lost.

------------------------------------------------------------------------

## Power Comparisions

`power_logrank_pfs()` calculates the power of detecting differences between control and treatment arms using the log rank test. This test is used for the Kaplan-Meier estimator. `power_copula_pfs()` calculates the power as well, but uses the bootstrap method to detect the difference between survival curves.

There are 7 different scenarios we wish to test.

| Scenarios | Description |
|----|----|
| 1 | No difference between control and treatment arms ($\beta=0.0$) |
| 2 | Small treatment effect ($\beta=-0.2$) |
| 3 | Medium treatment effect $\beta=-0.5$) |
| 4 | Strong treatment effect $\beta=-0.8$) |
| 5 | Early difference between control and treatment arms for a fixed $\beta=-0.5$ |
| 6 | Late difference between control and treatment arms for a fixed $\beta=-0.5$ |
| 7 | Crossing between control and treatment arms for a fixed $\beta=-0.5$ |

## Testing

Tests are written with `testthat`. Run them with:

``` r
devtools::test()
```

Slow bootstrap tests are marked `skip_on_cran()` and use a small `B` for speed. Coverage can be checked with:

``` r
covr::report()
```
