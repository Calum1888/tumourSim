# Copula Modelling for Mixed Continuous and Binary Variables in Survival Analysis

### Applications to Colorectal Cancer Clinical Trials

> A novel statistical framework that improves survival prediction in oncology trials by modelling the *dependence* between tumour growth and the appearance of new lesions — going beyond what Kaplan–Meier and augmented-binary methods can offer.

------------------------------------------------------------------------

## Overview

Clinical trials in oncology typically produce two distinct types of progression signals for each patient:

1.  **A continuous measurement** — the sum of the longest tumour diameters (SLD), tracked longitudinally.
2.  **A binary indicator** — whether a new lesion has appeared at a follow-up visit.

Under the **RECIST** framework, *either* of these crossing a threshold counts as a progression event. Standard methods (Kaplan–Meier, augmented-binary) treat the composite endpoint marginally and **discard the joint structure** between these two signals.

This project develops and evaluates a **copula-based joint survival model** that explicitly captures the dependence between the two endpoints, yielding:

-   More precise survival estimates than Kaplan–Meier in dependent-event settings
-   Increased statistical power to detect treatment effects when tumour growth and new lesions are linked
-   Interpretable measures of clinical dependence (Kendall's τ, tail dependence) inaccessible to traditional methods
-   Robustness to non-Gaussian tumour dynamics where the augmented-binary method breaks down

The framework is validated on a comprehensive simulation study and applied to the **Phase III PRIME colorectal cancer trial** (Amgen, n = 442 patients with full longitudinal data).

------------------------------------------------------------------------

## Key Results

### Headline finding on real trial data

In a stratified analysis of the PRIME colorectal cancer trial, the copula-based permutation test detected a **statistically significant late divergence** between FOLFOX alone and Panitumumab + FOLFOX in the KRAS-Mutant subgroup (**p = 0.008**), where the standard log-rank test found no difference (p = 0.879). This reflects the copula method's sensitivity to non-proportional hazards, which the log-rank test under-weights at late follow-up times.

### Simulation study highlights

| Scenario | Kaplan–Meier | Augmented-Binary | **Copula (Gumbel)** |
|----|:--:|:--:|:--:|
| Type I error (null) | 0.040 | 0.045 | **0.045** ✓ |
| Heavy-tail null | 0.055 | **0.174** ✗ inflated | **0.046** ✓ |
| Small treatment effect | 0.371 | 0.405 | **0.395** |
| Dependent endpoints + small effect | 0.974 | 0.990 | **0.993** ✓ best |
| Dependent endpoints + medium effect | 0.996 | 0.999 | **0.998** |

The copula approach **controls Type I error under heavy-tailed tumour dynamics**, where the augmented-binary method fails (inflated to 17.4%), while delivering competitive or superior power when endpoints are genuinely dependent.

### Survival curve estimation (single-arm trial, 7 time points)

<p align="center">

<img src="KM_GRAPH_7.png" alt="Kaplan-Meier" width="32%"/> <img src="AUGBIN_GRAPH_7.png" alt="Augmented-Binary" width="32%"/> <img src="COPULA_GRAPH_7.png" alt="Gumbel Copula" width="32%"/>

</p>

The Gumbel copula consistently achieves **narrower confidence intervals than Kaplan–Meier** while maintaining ≥95% coverage across all time points.

------------------------------------------------------------------------

## Real-Data Application: PRIME Colorectal Cancer Trial

<p align="center">

<img src="km_trt_by_KRAS.png" alt="Kaplan-Meier by KRAS type" width="48%"/> <img src="cop_trt_strat_by_KRAS.png" alt="Copula by KRAS type" width="48%"/>

</p>

*Left: Kaplan–Meier estimates by treatment, stratified by KRAS type. Right: Gumbel copula estimates of the same. The copula method reveals a sustained late survival advantage for Panitumumab + FOLFOX in the KRAS-Mutant subgroup between years 1 and 2 — a signal the log-rank test misses because of how it weights time points.*

The model also produces **clinically interpretable dependence measures**: - Kendall's τ = 0.256 (FOLFOX alone) vs τ = 0.307 (Panitumumab + FOLFOX) — patients in the combination arm show stronger co-occurrence of tumour growth and new lesion events. - Lower-tail dependence (survival-time scale) → patients who progress early tend to do so on *both* endpoints simultaneously, consistent with shared underlying disease biology.

------------------------------------------------------------------------

## Methodology

### The copula-based survival estimator

For lesion failure time $T_D$ and tumour-progression time $T_Y$, the progression-free survival function is reconstructed via a **survival copula**:

$$
\widehat{S}_{\text{PFS}}(t) = \widehat{S}_D(t) + \widehat{S}_Y(t) - 1 + C\big(1 - \widehat{S}_D(t),\, 1 - \widehat{S}_Y(t);\, \widehat{\theta}\big)
$$

where: - $\widehat{S}_D$ and $\widehat{S}_Y$ are Kaplan–Meier estimates of the marginal survival functions - $C(\cdot, \cdot; \theta)$ is an Archimedean copula (Clayton, Frank, or Gumbel) capturing the dependence - $\widehat{\theta}$ is estimated by maximum likelihood on the bivariate pseudo-observations

Inference uses **bootstrap resampling** (200 samples × 1000 iterations) for confidence intervals, and a **permutation test** on the integrated $L_1$ distance between arm-specific survival curves for hypothesis testing.

### Why Gumbel?

Across simulation scenarios and the real PRIME data, the Gumbel copula consistently produced the lowest AIC and best-calibrated coverage. Its lower-tail dependence (on the survival-time scale) matches the clinical reality that aggressive disease drives both tumour growth and new-lesion formation together.

------------------------------------------------------------------------

## Repository Structure

```         
.
├── DATA_GENERATION.R   # Simulates RECIST-style trial data (single & two-arm)
├── KM.R                # Kaplan-Meier marginals + log-rank baseline
├── AUGBIN.R            # Augmented-binary method (Lin & Wason, 2020)
├── COPULA.R            # Copula fitting, PFS estimator, bootstrap CIs, permutation test
├── HAZARDS.R           # Cox PH and AFT model helpers (illustrative)
├── POWER.R             # Power simulations across 10 scenarios
├── Final_Submission.pdf  # Full dissertation (83 pages)
└── *.png               # Survival curve plots
```

The estimation code is also distributed as an R package, **`tumourSim`**, available at\
👉 [github.com/Calum1888/tumourSim](https://github.com/Calum1888/tumourSim)

------------------------------------------------------------------------

## Technical Stack

-   **Language:** R
-   **Core packages:** `copula`, `survival`, `mvtnorm`, `boot`, `ggplot2`
-   **Statistical techniques:**
    -   Maximum likelihood estimation for Archimedean copulas
    -   Non-parametric bootstrap for confidence intervals
    -   Permutation testing for non-proportional hazards
    -   Two-step estimation framework (Shih & Louis, 1995)
    -   Delta method for variance approximation

------------------------------------------------------------------------

## Reproducing the Results

``` r
# Install the companion package
devtools::install_github("Calum1888/tumourSim")

# Or source the scripts directly
source("DATA_GENERATION.R")
source("KM.R")
source("COPULA.R")
source("POWER.R")

# Run a single-arm simulation
trial  <- simulate_single_arm(n = 150, n_times = 5)
result <- pfs_copula(trial, family = "gumbel", n_boot = 200)

# Run a two-arm power study (1000 iterations)
power_results <- run_power_simulation(scenario = 6, n_sim = 1000)
```

------------------------------------------------------------------------

## What I Built and Learned

This project was my Master of Mathematics (MMATH) final year project. It involved:

-   **Theoretical work:** deriving the survival-copula representation, computing tail dependence and Kendall's τ for Archimedean families, proving convergence of the t-copula to the Gaussian copula as ν → ∞.
-   **Software engineering:** writing a documented, tested R package (`tumourSim`) implementing data generation, three estimation methods, bootstrap inference, and a permutation test.
-   **Simulation design:** ten carefully constructed scenarios isolating treatment effect size, tumour–lesion dependence, and distributional assumptions (Gaussian vs heavy-tailed).
-   **Applied analysis:** end-to-end analysis of a real Phase III oncology trial, including endpoint construction from raw longitudinal data, model selection via AIC, and clinical interpretation of dependence measures.
-   **Communication:** an 83-page dissertation including a self-contained primer on copula theory, survival analysis, and the proposed methodology.

The work demonstrates that **modelling dependence between competing endpoints is not just theoretically interesting — it provides measurable gains in efficiency and interpretability in real clinical trial data**.

------------------------------------------------------------------------

## References

The complete reference list (38 entries) is included in the dissertation. Key methodological foundations:

-   Sklar, M. (1959). *Fonctions de répartition à n dimensions et leurs marges.*
-   Kaplan, E.L. & Meier, P. (1958). *Nonparametric Estimation from Incomplete Observations.* JASA.
-   Lin, C-J. & Wason, J. (2020). *Efficient analysis of time-to-event endpoints when the event involves a continuous variable crossing a threshold.* J. Stat. Plan. Inference.
-   Nelsen, R.B. (2006). *An Introduction to Copulas.* Springer.
-   Eisenhauer, E.A. et al. (2009). *New response evaluation criteria in solid tumours: Revised RECIST guideline (v1.1).* Eur. J. Cancer.

------------------------------------------------------------------------

## Contact

Author: **Calum Regan**\
MMATH Mathematics, 2022-2026\
📄 [Full dissertation (PDF)](Final_Submission.pdf)
