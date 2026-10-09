BASELINE: PREDICTIVE VALUE OF INDIVIDUAL-LEVEL INFORMATION
=======================================================

This adds the manuscript's M0--M3 comparison to the baseline STAR simulation:
J=40, N_j=150, c=1, Bernoulli(0.3) treatment, and the existing Gaussian DGP.
Existing simulation and application scripts and outputs are unchanged.
Sourcing the new runner loads functions only; it never starts a production run.

MODELS AND TARGET
-----------------
H=(1,H1,H2), where H1=C and H2=M in the existing simulation order.
M0: existing Gaussian microdata hierarchical model.
M1: mean H'a + b*p; common residual variance.
M2: mean H'a + p*H'b; common residual variance.
M3: same mean as M2; residual variance (1,p) Sigma (1,p)'.

M0 and M2/M3 coefficients have matching independent N(0,0.5^2) priors.
M3 and M0 use half-normal(0,0.5) marginal SDs and LKJ(2) correlations.
M1/M2 use a half-normal(0,0.5) common residual SD.
The aggregate Stan program is reused directly from STAR_model_comparison.stan.
M0 uses STAR_superpopulation_model.R in this folder, not the application's
Student-t model. No fold-specific rescaling is performed.

Default target="population": p=pF for ALL training and held-out schools,
and oracle truth is tau_pop_true. This matches the original simulation's
primary target. Set target="finite" to use realized sample free-lunch
proportions throughout and evaluate tau_finite_true. Target-specific fit
folders prevent mixing the two experiments; they share master populations.

For each training school, estimate = mean(y_t)-mean(y_c), with sampling SE
sqrt(var(y_t)/n_t + var(y_c)/n_c). The likelihood is approximately
Normal(estimate | mean, sampling_se^2 + residual_variance).
At least two observations per arm are required; invalid schools cause an error
rather than being dropped or redrawn. Only the held-out covariates enter
prediction. Held-out estimates/SEs never enter a training likelihood.

The default population comparison treats the difference-in-means estimate
as approximately normal about the population ATE under random sampling and
assignment. It does not assume the realized sample composition equals pF.
Plug-in SEs absorb sampling/assignment uncertainty approximately. In the finite
comparison the sampling model instead approximates randomization uncertainty
about the finite-sample ATE. M0/M3 share the mean and residual variance forms,
but differ in likelihood, outcome adjustment, and plug-in SE assumptions.
Interpret predictive differences conditional on these choices.

Unlike the empirical application, simulation truth is observed. Primary scores
are log predictive densities evaluated at that truth, with oracle coverage and
RMSE, rather than scores averaged over a noisy school-only reference posterior.
This avoids reference-estimation noise. Empirical reference-posteriors and
oracle scores are different validation procedures and should not be pooled.

RUN FROM THE CODE DIRECTORY
---------------------------
setwd("SimulationBHM_ModelComparison")
source("STAR_baseline_model_comparison.R")

# Base-R validation, with synthetic checkpoints (no MCMC):
source("test_STAR_baseline_model_comparison.R")

# Small MCMC pilot: four fits, one replication, one held-out school.
run_baseline_model_comparison(list(
  R_rep=1L, replications=1L, schools=1L,
  chains=2L, iter=500L, warmup=250L,
  output_dir="STAR_baseline_comparison_pilot"))

Inspect diagnostics before production. The pilot is not a converged production
result. Use a separate production directory so short pilot chains are not reused.

# Production: 100 replications x 10 schools x 4 models = 4,000 LOSO fits.
run_baseline_model_comparison()

# Score completed checkpoints without fitting or requiring RStan:
run_baseline_model_comparison(list(fit=FALSE))

# Alternative finite-sample estimand:
run_baseline_model_comparison(list(target="finite"))

# Reuse previously generated master populations, without redrawing them:
run_baseline_model_comparison(list(
  population_dir="../STAR_superpopulation_fits",
  output_dir="STAR_baseline_comparison_existing_populations"))
Adjust population_dir to the actual directory containing rep_XXXX folders.
All requested master_population.rds files must exist. The comparison makes its
own M0 LOSO fits, preserving the original checkpoints and recording exact
training-data/model signatures for all four models. No full-sample fits or
nine-scenario estimation are required. Default newly generated masters use the
original population seed, seed=20260828+r-1, and generator.

RESULTS AND RESTARTS
-------------------
STAR_baseline_model_comparison/populations/rep_XXXX/master_population.rds
STAR_baseline_model_comparison/population/rep_XXXX/school_inputs.csv
STAR_baseline_model_comparison/population/rep_XXXX/M0_school_JJ.rds
(similarly M1, M2, M3; finite/ is used for target="finite")
STAR_baseline_model_comparison/population/metrics/
  fold_metrics.csv
  replication.csv
  summary.csv
  paired_replication.csv
  paired_summary.csv
  diagnostics.csv
  metadata.rds

Each fit retains predictive locations/SDs, posterior coefficients/covariance,
HMC diagnostics, seed, sampler configuration, timing, and exact signature.
Completed compatible fits are reused. Changing sampler settings does not refit
completed checkpoints: choose a new output_dir or set overwrite=TRUE explicitly.
Population provenance changes require a new output directory. Changed training
data/model code are rejected unless overwrite=TRUE. Atomic writes reduce the
risk of partial checkpoints. Do not run simultaneous jobs on the same fit paths.

Predictive distributions are posterior normal mixtures for LATENT ATEs,
excluding held-out sampling variance. Density averaging uses log-sum-exp;
interval endpoints use the mixture CDF; predictive variances use the exact
mixture moment formula. No predictive replication noise is introduced.

Metrics first average across the selected schools within each replication,
then across independent replications with equal replication weights. summary.csv
reports mean replication RMSE (and mean MSE separately), coverage, widths,
log score, bias, and predictive SD, together with replication SDs and MCSEs.
paired_summary.csv reports model-minus-comparator log-score differences,
exp(mean difference), and MCSEs computed from paired replication means.
Positive score differences favor the named model. LOSO folds are not treated
as independent observations. With one replication, MCSE is undefined (NA).
Partial runs score ONLY their selected replications/schools; metadata records
that selection. Metrics files are replaced when recalculated for a selection.
Use separate output_dir values if multiple partial summaries must be retained.

VALIDATION STATUS
-----------------
Base-R tests cover mean/variance identities, population versus finite targets,
held-out outcome exclusion, analytic normal scores/quantiles, density scaling,
paired replication MCSEs, scoring-only exports, and stale-checkpoint rejection.
RStan is required only for missing fits. This implementation environment lacks
RStan, so no Stan compilation, MCMC pilot, or production fits were executed.
Review R-hat, ESS, divergences, treedepth hits, and E-BFMI before interpreting
any production rankings.

BUNDLED FILES
-------------
This folder is self-contained: it includes the comparison runner, helpers,
tests, Gaussian microdata model, aggregate Stan model, population generators,
STAR CSV data, and this documentation. All source/model paths resolve within
this folder. R and RStan/C++ remain external software requirements.

Copied from SimulationBHM and the Code-level STAR aggregate-model utilities.
Only dependency paths in the copied files were adapted; source files are intact.
