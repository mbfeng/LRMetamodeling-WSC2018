# Unbiased Metamodeling via Likelihood Ratios — code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24446361.svg)](https://doi.org/10.5281/zenodo.24446361)

R code for the numerical experiments of

> Dong, J., Feng, M. B., and Nelson, B. L. (2018). Unbiased Metamodeling via Likelihood Ratios.
> In *Proceedings of the 2018 Winter Simulation Conference*, 1778–1789. IEEE.
> https://doi.org/10.1109/WSC.2018.8632506

The paper studies likelihood-ratio (LR) metamodels: simulation outputs at one or several design
points are reweighted by likelihood ratios to predict the mean output at other parameter values.
Section 4 compares the metamodels on a stochastic activity network (SAN) with five exponential
activity times.

## Quick start

Requires R (tested with 4.3.3 and 4.4.0) and the packages `ggplot2` and `lattice`.

```sh
cd Code
Rscript reproduce_paper.R                 # Figures 1 and 2 from the saved summaries (~10 s)
Rscript tests/test_against_original_code.R
```

`reproduce_paper.R` redraws both figures of the paper (Figures 1 and 2) from the summaries in
`results/published_2018/` into `figures/` (PDF and
JPEG, with the file names of the paper's LaTeX source), then runs `verify_against_paper()`.
Set `RERUN_SIMULATIONS <- TRUE` to rerun both experiments first (see Run times).

## Paper result -> code

| Paper | What it shows | Script (simulation) | Function (figure) | Summary |
|---|---|---|---|---|
| Figure 1(a) | true response surface mu(theta1, theta5) | input `data/global_true_values.csv` | `plot_global_metamodels()` | `global_metamodel_summary.csv` |
| Figure 1(b) | variance (MSE) of the SMC estimator, 100 reps | `scripts/run_global_metamodels.R` | `plot_global_metamodels()` | same |
| Figure 1(c)-(e) | MSE of BLR, SLR, RLR relative to Var(SMC) (BLR truncated at 100) | `scripts/run_global_metamodels.R` | `plot_global_metamodels()` | same |
| Figure 1(f)-(h) | MSE of (SMC + LR)/2 relative to Var(SMC) (SMC + BLR truncated at 10) | `scripts/run_global_metamodels.R` | `plot_global_metamodels()` | same |
| Figure 2(a) | MSE vs. iteration t: BLR, SLR, RLR with PW, Crude, Brute | `scripts/run_target_theta.R` | `plot_target_theta()` | `target_theta_summary.csv` |
| Figure 2(b) | RLR with PW, OW, M, Crude, Brute | `scripts/run_target_theta.R` | `plot_target_theta()` | same |
| Figure 2(c) | BLR, SLR, RLR with M, Crude, Brute | `scripts/run_target_theta.R` | `plot_target_theta()` | same |

## How the simulations work

| | Section 4.1 (global, Figure 1) | Section 4.2 (target theta0, Figure 2) |
|---|---|---|
| Model | SAN, Y = max(X1+X4, X1+X3+X5, X2+X5), X_k ~ Exp(rate theta_k) | same |
| Design point(s) | theta* = (1,1,1,1,1) | thetahat_t = 1 / (mean of the first 10t observations), t = 1..100 |
| Prediction point(s) | 20 x 20 grid of (theta1, theta5) in [0.5, 1.5]^2, other rates 1 | theta0 = (1,1,1,1,1) |
| Replications | 100 at theta* per grid point (LR), 100 at each grid point (SMC) | n = 50 at thetahat_t per step; Brute: 50t new reps at thetahat_t |
| Estimators | SMC, BLR, SLR, RLR, (SMC + LR)/2 | Crude, Brute, {BLR, SLR, RLR} x {PW, OW, M} |
| True values | 10^7 replications per grid point (2018 run, `data/global_true_values.csv`) | 83/24 (closed form for theta0 = 1) |
| Macro replications | 10,000 | 10,000 |
| One macro replication | `global_macro_run()` in `R/experiments.R` | `target_macro_run()` in `R/experiments.R` |
| Output | MSE per grid point = mean over macro runs of (estimate - true)^2 | MSE per step t |

As in the 2018 code, each grid point gets its own sample of 100 replications at theta* in
Section 4.1; the MSE at each grid point does not depend on whether samples are shared across grid
points.

## Notation (code -> paper)

| Code | Paper |
|---|---|
| `theta`, `theta_center`, `theta0` | theta, theta* (Section 4.1), theta_0 (Section 4.2) |
| `X`, `Y` | activity times X, time to completion Y = g(X) |
| `density` | f(X given theta) |
| `lr` | likelihood ratio f(X given theta0) / f(X given theta_j) |
| `BLR`, `SLR`, `RLR` | baseline, self-normalized and regression LR metamodels (Section 3.1) |
| `PW`, `OW`, `SM` | proportional weights n_j/n (equal here), optimal weights (estimated inverse variances), mixture LR ("M") (Section 3.2) |
| `n_data`, `n_reps`, `n_steps` | m = 10, n = 50, t = 1, ..., 100 |
| `mse`, `se_mse`, `top1_share` | MSE, its standard error, share of the MSE from the worst 1% of macro runs |

## Folder structure

```
reproduce_paper.R        entry point (figures + verification)
R/san_model.R            SAN simulation, densities, BLR/SLR/RLR
R/experiments.R          one macro replication of each experiment
R/summaries.R            MSE, standard error, tail share over macro replications
R/plots.R                Figures 1 and 2
R/verify.R               verify_against_paper()
scripts/                 the two simulation drivers (resumable, command-line arguments)
data/                    true values of Section 4.1 and the 2018 random-number state
results/published_2018/  summaries extracted from the published 2018 runs
results/rerun*/          summaries and logs of the reruns (November 2026)
figures/                 output of reproduce_paper.R
tests/                   new code vs. the 2018 code on the same random numbers
tools/                   one-off extraction of the published summaries
original_2018/           the 2018 scripts (reference for the tests)
```

## Reproducibility

Both drivers draw their random numbers in the same order as the 2018 code, so they reproduce the
published runs exactly:

* `run_target_theta.R --seed 0` (default): the 2018 run called `set.seed(0)`.
* `run_global_metamodels.R --seed 2018` (default): the 2018 script called `set.seed(1)` and then
  loaded the workspace with the true values, which replaced the random-number state with the one
  saved in that workspace; this state is in `data/global_rng_state_2018.txt`.
* Any other `--seed` value starts a new, independent run (`set.seed(value)`).

Runs are serial (one random-number stream per experiment); each driver saves a checkpoint every
500 macro replications and resumes from it.

## Verification

* **Tests** (`tests/test_against_original_code.R`): the rewritten functions and one macro
  replication of each experiment agree with the 2018 functions on the same random numbers (largest
  relative difference 7e-16).
* **Published runs:** `tools/extract_published_runs.R` reduced the saved 2018 workspaces
  (`SANGlobal_100Samp_10000Macro.RData`, 300 MB; `SAN_10K_Macro.RData`, 123 MB) to the summaries in
  `results/published_2018/`; it checks that they give the same numbers as the 2018 plotting
  scripts and the saved MSE matrix. The redrawn figures agree with the published ones (side-by-side
  comparison; the surfaces are identical, and the MSE curves differ only in font size because
  ggplot2 has changed since 2018).
* **Statements of the paper:** `verify_against_paper()` checks 16 statements that Sections 4.1 and
  4.2 make about the figures (ranges, orderings, where an estimator beats SMC); all hold.
* **Full reruns (2018 random numbers):** both experiments rerun with 10,000 macro replications
  (`results/rerun/`) reproduce the published runs exactly: every MSE, standard error and tail
  share agrees to a relative difference of 6e-14 (Figure 1, 400 grid points x 7 estimators) and
  2e-14 (Figure 2, 100 steps x 12 estimators).
* **Second seed:** independent reruns with `--seed 2026` (`results/rerun_seed2026/`) agree with
  the published runs within Monte Carlo error: of 3,889 compared values, 33 have |z| > 2.5
  (48 expected by chance), the largest |z| is 3.4, and all 16 statements about the figures hold.
  `comparison_with_2018.csv` in each rerun folder lists every z-score.

## Run times

On a shared Linux server (Intel Xeon E5-2667 v4, 3.2 GHz; R 4.3.3), one core per run, four runs
at a time on a busy machine (November 2026):

| Run | 10,000 macro replications |
|---|---|
| `run_global_metamodels.R` (Figure 1) | 13–25 minutes |
| `run_target_theta.R` (Figure 2) | 38–60 minutes |

The lower times are from a less loaded period. On a Windows laptop (R 4.4.0) one macro
replication takes about 0.07 s (global) and 0.3 s (target); the 2018 code took about 0.7 s and
3.4 s. `reproduce_paper.R` without reruns takes about 10 seconds.

## Provenance

* Figure 1: `Codes/SAN_Global/` of the project folder (`SANGlobal_main.R` of 2018-04-08 with
  10^4 macro replications; true values from `SANGlobal_True0.RData`, 10^7 replications per grid
  point). Figure 2: `Codes/SAN_Target/` (`SANEst_helpers.R` and the workspace `SAN_10K_Macro.RData`
  of 2018-03-16; the script in the folder shows 10 macro replications and `set.seed(1)`, the saved
  run has 10^4 and was made with `set.seed(0)`, found by matching the first macro replication).
  The final figures were drawn on 2018-07-22 with `SANGlobal_plot.R` and `SANEst_plot.R`.

## Change log

* 2026-11 (Ben Feng, with Claude Code): code restructured into one entry point, a function
  library and two resumable drivers; plain R with vectorized densities (about 10 times faster than
  the 2018 code, same random-number order); figures drawn from compact summaries; tests against the
  2018 code; reruns and verification. The first server rerun of Section 4.2 differed from the
  2018 run in the PW- and OW-combined RLR estimators (up to 0.6% in MSE): the joint density was
  computed as `X %*% theta`, which an optimized BLAS can round differently for the same row in
  matrices of different sizes, so the likelihood ratios of the current step were 1 +- 1e-16
  instead of exactly 1 and the RLR regression slope of that step became noise. `san_density()`
  now sums the five columns explicitly; the tests then pass with both the reference BLAS and
  OpenBLAS, and the reruns were repeated.
* 2018-03/07: original code.

## Archived code

This code is archived on Zenodo: https://doi.org/10.5281/zenodo.24446361 (all versions;
v1.0.0 is https://doi.org/10.5281/zenodo.24446362).

## License

MIT (see `LICENSE`).
