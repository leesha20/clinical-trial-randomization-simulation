# Comparing Clinical Trial Randomization Designs via Simulation

**Author:** Sriya Leesha Gourammagari
**Course:** BIMS/STAT 5304 (Introduction to Human Health Research) — Final Project, UT Dallas

*A 1,000-replicate Monte Carlo on randomization, balance, and what logistic regression quietly hides.*

## Project Overview

This project compares three clinical-trial randomization designs — simple randomization, Pocock–Simon minimization, and randomized play-the-winner — crossed with two analysis strategies (unadjusted vs. covariate-adjusted logistic regression). Because outcomes are simulated from a known data-generating process, every design-by-analysis combination can be scored against ground truth on bias, efficiency, confidence-interval coverage, covariate balance, and power. The central result is that, at realistic trial sizes, *the analysis choice controls whether the estimate is unbiased, while the design choice mostly controls covariate balance and small-sample power.*

## The Question and the Ground Truth

A new drug works; the true treatment effect is fixed by the simulation at **β_drug = +1.00** (odds ratio ≈ 2.72). Outcomes are drawn from a logistic model in which two covariates are strongly prognostic and two are weak noise:

```
logit P(Y = 1) = −0.30 + 1.00·Drug + 3.00·Age65+ − 3.00·Severe − 0.05·BMI − 0.03·SBP0
```

Every estimate that follows is measured against β_drug = 1.00. Age and disease severity are strong prognostic factors (which is exactly why balance and adjustment matter); BMI and baseline SBP are weak nuisance variables. Each simulated trial enrolls N = 200 patients; the study runs 1,000 replicates, with a power sweep from N = 50 to 300.

## Designs (implemented from scratch)

- **Simple** — independent Bernoulli(0.5) allocation; a fair coin per patient.
- **Pocock–Simon minimization** — a biased coin (p = 0.75 toward the less-imbalanced arm) that minimizes current imbalance on age and severity at every assignment. Goal: covariate balance and precision.
- **Randomized play-the-winner (RPW)** — an urn rule (β = 3; success adds to the same arm, failure to the opposite) that drifts allocation toward the winning arm as outcomes accrue. Goal: ethical allocation.

## Estimation and Evaluation

Each design is analyzed two ways: a **naive** model (`Y ~ Drug`) and a **covariate-adjusted** model (`Y ~ Drug + Age + Severity`). Model selection is additionally explored with stepwise AIC and LASSO (`glmnet`). Across 1,000 replicates, the study measures bias relative to β_drug, empirical standard error, 95% CI coverage, the distribution of final allocation, covariate imbalance, and a power curve across sample sizes; a bias²/variance decomposition summarizes mean squared error.

## Key Results

**1. The bias is in the analysis, not the design (non-collapsibility).** The naive model attenuates the treatment effect to a mean β̂ ≈ 0.68 in *every* design — including perfectly-balanced Pocock — a ~32% attenuation. This is not omitted-variable bias from imbalance; it is non-collapsibility of the logistic model: marginalizing over strong prognostic covariates shrinks the log-odds ratio toward zero, so the naive model estimates the (smaller) *marginal* effect rather than the *conditional* truth. Covariate adjustment recovers β̂ ≈ 1.00 with near-nominal coverage across all three designs.

| Design × Model | Mean β̂ | Empirical SE | 95% Coverage |
|---|---|---|---|
| Simple · Naive | 0.68 | 0.299 | 78.5% |
| Simple · Adjusted | 0.99 | 0.380 | 94.2% |
| Pocock · Naive | 0.69 | 0.257 | 86.0% |
| Pocock · Adjusted | 1.03 | 0.381 | 93.6% |
| RPW · Naive | 0.68 | 0.301 | 80.0% |
| RPW · Adjusted | 1.01 | 0.382 | 94.2% |

**2. Each design has a distinct allocation and balance "personality."** Final proportion on drug across replicates: Pocock is a needle on 50% (SD 0.006), Simple is a clean bell around 50% (SD 0.036), and RPW drifts to ≈60% (SD 0.074, ~12× Pocock's spread). Mean total covariate imbalance at N = 200 tells the same story in reverse: Pocock 8.2, Simple 38.8, RPW 86.2 — RPW buys ethical allocation at the cost of the worst balance.

**3. Adjustment cuts MSE ~25%; design barely moves it.** In the bias²/variance decomposition, adjusted models drive bias² to ≈0 at a small variance cost, lowering MSE across all designs; the randomization design has little effect on MSE once the analysis is correct.

**4. Design matters most when under-powered.** On the adjusted model, Pocock reaches ~74% power at N = 150 vs ~64% for Simple/RPW — a ~10-point edge that vanishes by N = 200, where all three designs first cross 80%.

**5. A cautionary result on data-driven variable selection.** Letting a penalty pick covariates can silently destroy the treatment estimate:

| Method | β̂_drug | vs. truth (1.00) |
|---|---|---|
| Full model | 0.66 | −34% |
| Stepwise AIC | 0.68 | −32% |
| LASSO (λ.min) | 0.49 | −51% |
| LASSO (λ.1se) | 0.04 | −96% |

The 1-SE rule is a prediction heuristic; used for inference it shrank a real effect nearly to zero. The penalty should be chosen for the question being asked.

## Repository Structure

```
.
├── Final_project.R              # Single-run: three designs, EDA, naive/adjusted models, model selection
├── Multireplicate_simulation.R  # 1,000-replicate Monte Carlo: bias, SE, coverage, allocation, imbalance
├── Build_tables.R               # Presentation tables
├── Bonus_Graphs.R               # Power curves, imbalance over time, bias²/variance decomposition
├── presentation.pdf             # Slide deck summarizing the study
├── README.md
└── output/                      # Generated tables and figures
```

## Reproducibility

1. Clone this repository.
2. Install the required R packages: `ggplot2`, `dplyr`, `tidyr`, `scales`, `knitr`, `MASS`, `glmnet`, `broom`.
3. Run `Final_project.R` for a single simulated trial, then `Multireplicate_simulation.R` for the 1,000-replicate study; `Build_tables.R` and `Bonus_Graphs.R` produce the tables and figures.
4. All randomness is seeded (`set.seed(5304)`) for reproducibility.

## Notes

All data in this project is synthetic, generated within the scripts from the data-generating process above; no human-subjects data is used.
