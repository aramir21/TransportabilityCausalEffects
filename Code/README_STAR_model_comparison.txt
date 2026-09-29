STAR: MICRODATA VERSUS AGGREGATE-DATA LOSO COMPARISON
==================================================

Specification
-------------
M0 reuses the existing 79 microdata LOSO checkpoints in star_output/loso.
As confirmed for this exercise, all models preserve M0's teacher covariates:
C = prop_teacher_eth_cauc and M = prop_degree_master+ (the specific category).

Let h = (1, C, M)' on M0's training-fold standardized scale, and let p denote
the school's unstandardized free-lunch proportion. The models are:

  M1: mean = h'a + b*p;                    residual variance = omega^2
  M2: mean = h'a + p*h'b;                  residual variance = omega^2
  M3: mean = h'a + p*h'b;                  residual variance = (1,p) Sigma (1,p)'

The column order is C then M, matching the stored Gamma arrays. This is just
a permutation of the manuscript's M then C ordering. M2/M3 coefficient order
is vec(Gamma), with independent N(0, 0.5^2) priors, matching M0. M1 uses the
same coefficient prior. Common residual SDs have half-normal(0, 0.5) priors.
M3 uses half-normal(0, 0.5) marginal SDs and LKJ(2) correlation, matching M0.
Cholesky evaluation ensures the quadratic variance remains nonnegative.

For each training school, the aggregate likelihood is

  tau_hat_j ~ Normal(mu_j, s_j^2 + V_res,j).

The second argument here denotes variance. The known sampling SE is the
existing difference-in-means SE: sqrt(var(y_t)/n_t + var(y_c)/n_c).
This retains the normal approximation and sampling-variance convention of
the existing application, rather than introducing a new variance estimator.
Latent school effects are integrated out during fitting. Predictive densities
for latent ATEs include V_res only, not the held-out sampling variance.

School summaries are recovered in raw outcome units from each M0 checkpoint's
direct_ATE and outcome_scale. H values are recovered by inverting its stored
transformation. For fold j, all training summaries and covariates use fold j's
M0 outcome scale and H scaling. The held-out estimate and SE never enter the
training likelihood. All folds retain the same treatment and sample as M0.

Validation target
-----------------
All models are evaluated against the SAME school-only reference:

  tau_j | D_j ~ Normal(tau_hat_j, s_j^2).

This follows from the approximate normal sampling model and an improper flat
prior. No full-data hierarchical posterior is used in this reference.

The log score is E_reference[log(mean_r Normal_density(tau; mu_r, sd_r))].
It is computed by adaptive one-dimensional numerical integration with
log-sum-exp for the inner mixture. Integration error estimates are saved.
Predictive interval endpoints use the normal-mixture CDF; posterior coverage
uses the reference normal CDF. No predictive replications are needed.

Scores use the same full-fit control SD solely as a reporting unit, recovered
from the existing latent school metrics CSV. Scaling the outcome by a common
positive constant shifts all log scores equally and leaves paired differences
unchanged. The school-only reference is different from the previous full-fit
reference; the old -0.676602 score is not a comparable M0 baseline here.

Reported measures include log score, RMSE of predictive versus reference means,
posterior-expected RMSE, posterior coverage, binary reference-mean coverage,
and interval widths. Model summaries give equal weight to schools. Pairwise
log-score differences are model minus comparator, so positive favors model.
These are predictive score differences, not Bayes factors or posterior model
probabilities. No independent-fold significance test is supplied: LOSO folds
share training observations.

Matching M0/M3 mean and residual variance functions does not make their
likelihoods identical. The comparison retains the microdata model's outcome
adjustment and Student-t likelihood versus approximate normal aggregate
estimates with plug-in SEs. Interpret the information comparison conditional
on these modeling and sampling-variance assumptions. M3 covariance learning
may be weak if free-lunch composition varies little; inspect diagnostics.

Run from the Code directory
---------------------------
Current status: ONE fit exists for each of M1, M2 and M3, each using 78 schools
and omitting school 1. Thus three aggregate fits exist; 234 remain for the
other omitted schools. N = 78 in a checkpoint is the number of training
schools in that one fit. M0's 79 checkpoints already exist separately.

To run the COMPLETE workflow for all 79 schools:

  options(star.comparison.run = NULL)
  source("ApplicationSTARV3_model_comparison_run.R")

The runner explicitly selects all 79 schools regardless of earlier
star.comparison or star.comparison.metrics settings. Step 1 fits missing
M1--M3 folds and reuses compatible checkpoints. Step 2 calculates all metrics
for M0--M3, producing 316 model-school rows and overall summaries under
star_output/model_comparison/. No M0 refitting is needed. If interrupted,
rerun the same command; completed fits are reused and metrics recalculated.
RStan is loaded and compiled only when a fit is missing.

To run the complete workflow for school 1 only:

  options(star.comparison.run = list(schools = 1L))
  source("ApplicationSTARV3_model_comparison_run.R")
  options(star.comparison.run = NULL)

This test has completed: ALL school-1 metrics for M0--M3 are saved under
star_output/model_comparison/metrics_schools_1/. Its log scores agree with
the earlier pilot results. The remaining 234 fits were not run in this test.

The individual scripts can also be used manually. For school-1 fitting only:

  options(star.comparison = list(schools = 1L))
  source("ApplicationSTARV3_model_comparison_fit.R")

For all remaining fits, followed by all-school scoring:

  options(star.comparison = NULL)
  source("ApplicationSTARV3_model_comparison_fit.R")
  source("ApplicationSTARV3_model_comparison_metrics.R")

This would complete 79 fits per aggregate model (237 total), reusing the
three already saved. The metrics script calculates scores; it does not fit
models. Its default requires all 237 checkpoints; a subset requires only
that subset's checkpoints.

Fitting requires RStan and a working C++ toolchain. Scoring needs only base R.
No M0 refitting is required. Sampler settings and output paths can be overridden
through options(star.comparison = list(...)); see defaults in the fit script.
The metrics script uses options(star.comparison.metrics = list(...)).
It supports a schools subset and a separate fit_dir. For school-1 scoring only:

  options(star.comparison.metrics = list(schools = 1L,
    fit_dir = "star_output/model_comparison",
    output_dir = "star_output/model_comparison/metrics_schools_1"))
  source("ApplicationSTARV3_model_comparison_metrics.R")
  options(star.comparison.metrics = NULL)

To evaluate just the new M0 baseline in a separate output directory:

  options(star.comparison.metrics = list(models = "M0",
    output_dir = "star_output/model_comparison_M0"))
  source("ApplicationSTARV3_model_comparison_metrics.R")
  options(star.comparison.metrics = NULL)

Outputs of a completed all-school run under star_output/model_comparison:
  school_inputs.csv
  M1_school_01.rds, ... (similarly for M2 and M3)
  comparison_school_metrics.csv
  comparison_summary.csv
  comparison_pairwise_log_scores.csv
  comparison_coverage.csv
  comparison_diagnostics.csv
  comparison_metadata.rds

Checkpoints contain predictive component means/SDs, configuration, exact Stan
data/model signatures, and diagnostics. Changes to the model or training data
are rejected on resume; use a new output directory or explicit overwrite.
Existing M0 outputs and the earlier latent validation files are not modified.
Review R-hat, ESS, divergences, treedepth hits and E-BFMI before interpreting
rankings. Overall 79-school comparisons require all fits to complete.
Complete school-1 metrics are already available in metrics_schools_1/.

Checks:
  source("test_STAR_model_comparison.R")

These cover analytic normal log scores, density scaling, interval quantiles,
M0/M3 mean and variance identities, and exclusion of held-out outcomes.

School-1 heterogeneity comparison
--------------------------------
ApplicationSTARV3_school1_variance_comparison.R reports posterior covariance
parameters and the decomposition E[V_res | training] + Var(mu | training)
for the school-1 ATE. Results are written in metrics_schools_1/ as
posterior_variance_parameters.csv and predictive_variance_decomposition.csv.
All parameters are converted to the common full-fit control-SD scale.
These are LOSO parameter posteriors with school 1 excluded.

The original aggregate pilot checkpoints retained only predictive parameters.
For M1/M2, sds are the common residual SD draws. The full M3 covariance cannot
be inferred from its school-1 predictive SD alone. The script
recover_STAR_M3_school1_covariance.R reruns that fit with the saved specification
and seed, saving Sigma and L_Sigma draws in M3_school_01_covariance.rds without
replacing the original checkpoint. It also compares the sorted predictive
draws with the original results and saves HMC diagnostics.
New fits now retain coefficients and covariance parameter draws in posterior;
existing checkpoints are reused unchanged. L_Sigma is lower triangular with
Sigma = L_Sigma L_Sigma'. A posterior mean Cholesky factor multiplied by its
transpose is generally not the posterior mean covariance matrix.
