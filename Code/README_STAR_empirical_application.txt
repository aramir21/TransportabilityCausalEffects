PROJECT STAR EMPIRICAL APPLICATION
==================================

Purpose
-------
This application estimates heterogeneous treatment effects across schools in
the Tennessee Project STAR experiment. It also conducts leave-one-school-out
(LOSO) validation of the latent school average treatment effects (ATEs).

For the separate M0--M3 microdata-versus-aggregate model comparison using a
common school-only reference posterior, see README_STAR_model_comparison.txt.
Its fitting and metrics scripts are ApplicationSTARV3_model_comparison_fit.R
and ApplicationSTARV3_model_comparison_metrics.R. It preserves the existing
M0 covariates and fits, and writes to separate model_comparison output folders.
To fit missing aggregate folds and score all four models across all 79 schools:
  options(star.comparison.run = NULL)
  source("ApplicationSTARV3_model_comparison_run.R")

The maintained scripts are:

1. ApplicationSTARV3_clean.R
   Prepares the STAR data, estimates the full hierarchical model, reports
   posterior summaries, and constructs school-level ATE distributions.

2. ApplicationSTARV3_LOSO.R
   Estimates one model for each omitted school. With 79 observed schools, this
   produces 79 models based on 78 training schools each.

3. ApplicationSTARV3_LOSO_metrics.R
   Compares the latent ATE predicted before observing school j with the latent
   ATE posterior obtained after including school j in the full model.

4. PlotSTAR_school_LOSO_vs_full.R
   Produces an overlaid histogram for a selected school's LOSO latent ATE
   predictive distribution and full-data latent ATE posterior.

5. STAR_reference_categories.txt
   Documents the reference categories and interpretation of X, Z, W, and H.

6. ApplicationSTARV3_full_diagnostics.R
   Loads the saved full-sample Stan fit and produces parameter-level,
   chain-level, and overall convergence and HMC diagnostic tables.


Data and required software
--------------------------
The data are AER::STAR. The main R packages are:

  AER, dplyr, forcats, ggplot2, MASS, rlang, rstan, stringr, and tidyr.

The analysis uses RStan and therefore requires a working C++ toolchain.


Configuration
-------------
The main configuration block is at the beginning of
ApplicationSTARV3_clean.R. The default application uses:

  grade   = kindergarten
  subject = reading

Available grades are kindergarten, first, second, and third. Available subjects
are reading and math. The full fit and all LOSO fits must use the same grade,
subject, treatment definition, covariates, and model specification.


Treatment and estimand
----------------------
The treatment indicator is:

  Tr = 1: small class
  Tr = 0: regular class or regular class with an aide

The treatment contrast is therefore small class versus the combined non-small
class group.

For student i in school j, the heterogeneous treatment effect is

  tau_ij = Z_i' delta_j.

The average treatment effect for the observed composition of school j is

  tau_j = Zbar_j' delta_j,

where Zbar_j is the within-school average of the individual treatment-effect
moderators.


Regressors and reference categories
-----------------------------------
X contains individual-level outcome regressors:

  intercept, female, Caucasian, African American, standardized age,
  Q1, Q2, Q3, and free lunch.

The individual-level reference categories are male, an ethnicity other than
Caucasian or African American, Q4 birth, and non-free lunch. Age zero denotes
the estimation-sample mean age.

Z contains an intercept and free lunch. Thus, the first random treatment
coefficient is the treatment effect for a non-free-lunch student, and the
second is the free-lunch treatment-effect difference.

W contains standardized school-level predictors of the outcome. The omitted
compositional categories are non-Caucasian teacher ethnicity, apprentice career
ladder, bachelor's degree, inner-city location, and non-free lunch.

H contains an intercept, the standardized proportion of Caucasian teachers,
and the standardized proportion of teachers with a master's-plus degree.

Teacher-characteristic proportions are student weighted because AER::STAR does
not provide a teacher identifier.


Hierarchical model
------------------
The school-specific treatment coefficients satisfy

  delta_j = H_j Gamma + u_j,
  u_j ~ MVN(0, Sigma_delta).

The outcome model is

  y_ij ~ Student-t(5, mu_ij, sigma),

where

  mu_ij = X_ij' beta + W_j' gamma + Tr_ij Z_ij' delta_j.

The outcome is standardized using the control-group mean and standard deviation
in the final estimation sample. Treatment effects are therefore expressed in
control-group standard-deviation units.

The Stan model uses a non-centred Cholesky parameterization for the multivariate
school effects. In the current code, the standard deviations in Sigma_delta
have half-normal N+(0, 0.5^2) priors. Change this explicitly if the intended
specification is a half-Cauchy prior.


Full-sample analysis
--------------------
Run from the project directory:

  source("ApplicationSTARV3_clean.R")

Important outputs are written under star_output/, including:

  STAR_kindergarten_reading_fit.RData
  beta_summary.csv
  gamma_summary.csv
  Gamma_summary.csv
  Sigma_delta_summary.csv
  school_treatment_effects.csv
  external_validity_summary.csv
  external_validity_R2.pdf
  observed_school_mean_ATEs.pdf

Change the grade and subject in the configuration block when estimating a
different application. The saved-fit filename changes accordingly.


Leave-one-school-out estimation
-------------------------------
Run:

  source("ApplicationSTARV3_LOSO.R")

The Stan model is compiled once. Each fold removes every observation from one
school and recomputes all training-dependent transformations using only the
remaining 78 schools. These include:

  - the control-group outcome mean and standard deviation;
  - age standardization;
  - W standardization; and
  - H standardization.

Each completed fold is checkpointed under star_output/loso/:

  loso_school_01.rds
  ...
  loso_school_79.rds

The filenames use consecutive Stan school indices, not original STAR school
IDs. The original ID is stored in result$metadata$original_schoolid. This is
important because original school ID 77 is absent. No model is estimated for an
unobserved school; the 79 folds correspond to the 79 schools in the estimation
sample.

The LOSO script is resumable. Existing checkpoints are skipped unless
overwrite = TRUE. Set schools = 1L for a one-fold test and schools = NULL for all
unfinished folds.


Latent-ATE predictive validation
--------------------------------
Run:

  source("ApplicationSTARV3_LOSO_metrics.R")

For school j, the validation compares:

  p(tau_j | D_{-j}, H_j, Zbar_j)

with

  p(tau_j | D_{1:79}).

The first distribution is the latent ATE prediction before observing school j.
The second is the latent ATE posterior after including that school. The direct
difference-in-means estimate and its standard error are not added to the LOSO
distribution, because the validation target is tau_j rather than a noisy
estimate of tau_j.

Because every LOSO fold uses its own control-group standard deviation, the
script converts LOSO draws to the full-fit control-SD scale before comparison.

The reported metrics are:

  1. RMSE between the LOSO and full-data posterior means.

  2. Posterior-expected RMSE:

       sqrt(mean((m_j_LOSO - m_j_full)^2 + Var(tau_j | full data))).

  3. Posterior or soft coverage at 50%, 80%, 90%, and 95%:

       P(tau_j lies in the LOSO interval | full data).

  4. Binary coverage of the full-data posterior mean.

  5. Average and median lengths of the latent LOSO predictive intervals.

  6. Posterior-expected log predictive score: for each full-data ATE draw,
     evaluate the normal mixture density implied by all LOSO Gamma and
     Sigma_delta draws, take its natural logarithm, then average these logs.
     Both component means and SDs are converted to the full-fit control-SD
     scale. Computation uses log-sum-exp for numerical stability. The overall
     score weights schools equally; higher is better for a fixed full-data
     reference posterior. Scores depend on the units of the ATE. The KL
     comparison holds with that reference posterior fixed, not generally
     across schools or models with different full-data reference posteriors.

Outputs are written under star_output/loso_latent_metrics/:

  LOSO_latent_RMSE_summary.csv
  LOSO_latent_coverage_summary.csv
  LOSO_latent_school_metrics.csv
  LOSO_vs_full_latent_ATE_draws.rds
  LOSO_latent_log_score_summary.csv

If the latent metrics and comparison draws have already been saved, compute
only the log score without refitting or loading RStan:

  source("ApplicationSTARV3_LOSO_log_score.R")

This base-R script reads the existing comparison draws, school metrics, and
all 79 LOSO checkpoints. Use files from the same analysis run. It writes
LOSO_latent_log_score_by_school.csv and LOSO_latent_log_score_summary.csv.
Both scoring workflows use STAR_latent_predictive_helpers.R; numerical checks
can be run with source("test_STAR_log_score.R").


School-specific comparison figure
---------------------------------
After running the latent metrics script, run:

  source("PlotSTAR_school_LOSO_vs_full.R")

Set school_index in that file to select a school. For school 1, the output is:

  star_output/loso_latent_metrics/
  school_01_LOSO_vs_full_latent_ATE.pdf

The light-grey histogram is the latent LOSO predictive distribution before the
school is observed. The dark-grey histogram is the latent ATE posterior from the
full model after the school is included.


Full-fit convergence diagnostics
--------------------------------
After estimating or loading the full-sample model, run:

  source("ApplicationSTARV3_full_diagnostics.R")

The script reads:

  star_output/STAR_kindergarten_reading_fit.RData

It reports parameter-level R-hat and effective sample sizes, divergences,
maximum-treedepth hits, E-BFMI by chain, acceptance statistics, leapfrog-step
summaries, and step sizes. If the posterior package is installed, the parameter
table uses rank-normalized split-R-hat and bulk/tail effective sample sizes;
otherwise, it uses the diagnostics reported by rstan.

Outputs are written under star_output/full_fit_diagnostics/:

  full_fit_overall_diagnostics.csv
  full_fit_parameter_diagnostics.csv
  full_fit_chain_diagnostics.csv

The concise parameter table includes beta, gamma, Gamma, sd_delta, Omega,
Sigma_delta, and sigma. Standard rstan HMC diagnostic messages are also printed
in the R Console.


Reproducibility and diagnostics
-------------------------------
All scripts use explicit seeds. Before interpreting results, inspect split
R-hat, effective sample sizes, divergences, maximum treedepth, and E-BFMI for the
full model and every LOSO fit. The LOSO checkpoints contain a diagnostics list.

Do not compare results from fits that use different grades, subjects, treatment
definitions, covariate sets, priors, or standardization conventions.
