# Copula-Based Progression-Free Survival Estimation

An R package implementing a copula-based approach to estimate progression-free survival (PFS) by jointly modelling tumour progression and new lesion appearance, following the RECIST framework.

This work accompanies the dissertation *Applications of copula modelling to continuous and binary variables in survival analysis* (Calum Regan, 2025–2026).

---

## Background

In oncology trials, progression-free survival is defined as the time to the first of two competing events:

- **Tumour progression** — tumour size exceeds 20% of the subject's rolling minimum (RECIST 1.1 criterion)
- **New lesion appearance** — a binary indicator switches from 0 to 1

Traditional Kaplan–Meier methods treat PFS as a single combined endpoint. This package instead models the two endpoints separately with marginal Kaplan–Meier estimators and then reconstructs the joint PFS curve via the **survival-copula identity**:

$$S_{\text{PFS}}(t) = S_D(t) + S_Y(t) - 1 + C(1 - S_D(t),\, 1 - S_Y(t);\, \hat{\theta})$$

where $C$ is a bivariate copula fitted to pseudo-observations derived from the marginal event times.

---

## Installation

```r
# Install dependencies
install.packages(c("survival", "VineCopula", "MASS"))

# Install from source
devtools::install_github("Calum1888/tumourSim")
```

---

## Data Generation

Two functions simulate clinical trial data under the RECIST model.

### `generate_continuous_data()`

Simulates log tumour size ratios from a multivariate normal distribution and baseline tumour sizes from Uniform(0, 1).

```r
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

```r
coeffs <- generate_coefficients(
  n_times    = 5,
  n_patients = 150,
  alpha      = -2.5,   # baseline log-odds
  beta       = 0.5,    # treatment arm effect (set 0 for single-arm)
  gamma      = 0.3,    # tumour size effect
  R          = 0       # arm indicator: 0 = control, 1 = treatment
)
```

---

## Core Functions

### `lesion_event()`

Derives binary event times from a lesion indicator matrix. Each subject's event time is the first column index where the indicator equals 1; subjects with no event are censored at the final time point.

```r
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

```r
tumour_event(tumour_size_data, threshold = 1.2)
```

### `copula_pfs()`

The main estimator. Fits a bivariate copula to pseudo-observations from the marginal KM CDFs, then evaluates the PFS curve at each time point using the survival-copula identity.

```r
# Gaussian copula (family = 1)
pfs <- copula_pfs(
  lesion_events  = lesion_event(lesion_data),
  tumour_events  = tumour_event(tumour_size_data),
  n_times        = 5,
  copula_family  = 1
)
```

Copula family codes follow the `VineCopula` package convention. Passing a vector of families triggers AIC-based selection via `BiCopSelect()`. Recommended families for survival data:

| Family code | Copula  | Tail dependence        |
|-------------|---------|------------------------|
| 1           | Gaussian | None                  |
| 2           | t        | Upper and lower        |
| 3           | Clayton  | Lower only             |
| 4           | Gumbel   | Upper only             |
| 5           | Frank    | None (symmetric centre)|

### `copula_margin_estimation()`

Returns pointwise Kaplan–Meier survival estimates for each marginal endpoint at times `1:n_times`. Used internally by `copula_pfs()` but exposed for diagnostics.

---

## Bootstrap Inference

`bootstrap_copula_pfs()` draws `B` subject-level bootstrap resamples and computes pointwise percentile confidence intervals, coverage against a known true PFS curve, and average CI width.

```r
true_pfs <- exp(-0.15 * 1:5)   # known true curve (simulation context)

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

---

## Simulation Results

The table below reproduces the simulation results from the dissertation, comparing the copula estimator to the Kaplan–Meier estimator over 1000 replications with 150 patients.

**5 time points**

| Time | True S(t) | KM estimate | Copula estimate | KM coverage | KM avg CI width |
|------|-----------|-------------|-----------------|-------------|-----------------|
| 1    | 0.561     | 0.562       | —               | 0.955       | 0.159           |
| 2    | 0.304     | 0.307       | —               | 0.955       | 0.148           |
| 3    | 0.193     | 0.195       | —               | 0.961       | 0.128           |
| 4    | 0.101     | 0.103       | —               | 0.960       | 0.100           |
| 5    | 0.055     | 0.056       | —               | 0.954       | 0.078           |

**7 time points**

| Time | True S(t) | KM estimate | Copula estimate | KM coverage | KM avg CI width |
|------|-----------|-------------|-----------------|-------------|-----------------|
| 1    | 0.712     | 0.717       | —               | 0.946       | 0.144           |
| 2    | 0.473     | 0.478       | —               | 0.936       | 0.160           |
| 3    | 0.301     | 0.303       | —               | 0.949       | 0.148           |
| 4    | 0.162     | 0.163       | —               | 0.952       | 0.120           |
| 5    | 0.095     | 0.094       | —               | 0.958       | 0.097           |
| 6    | 0.058     | 0.058       | —               | 0.959       | 0.079           |
| 7    | 0.019     | 0.019       | —               | 0.978       | 0.055           |

---

## Testing

Tests are written with `testthat`. Run them with:

```r
devtools::test()
```

Slow bootstrap tests are marked `skip_on_cran()` and use a small `B` for speed. Coverage can be checked with:

```r
covr::report()
```

---

## Dependencies

| Package      | Role                                        |
|--------------|---------------------------------------------|
| `survival`   | Kaplan–Meier estimation via `survfit()`     |
| `VineCopula` | Copula fitting via `BiCopSelect()`, `BiCopCDF()` |
| `MASS`       | Multivariate normal simulation via `mvrnorm()` |

---

## References

- Eisenhauer et al. (2009). New response evaluation criteria in solid tumours: Revised RECIST guideline (version 1.1). *European Journal of Cancer*, 45(2), 228–247.
- Lin & Wason (2020). Efficient analysis of time-to-event endpoints when the event involves a continuous variable crossing a threshold. *Journal of Statistical Planning and Inference*, 208, 119–129.
- Nelsen, R.B. (2006). *An Introduction to Copulas*. Springer.
- Kaplan, E.L. & Meier, P. (1958). Nonparametric estimation from incomplete observations. *JASA*, 53, 457–481.
- Sklar, M. (1959). Fonctions de répartition à N dimensions et leurs marges. *Annales de l'ISUP*.
