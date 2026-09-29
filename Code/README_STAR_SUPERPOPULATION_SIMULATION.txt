STAR SUPER-POPULATION MONTE CARLO SIMULATION
============================================

PURPOSE
-------

This directory contains a three-stage Monte Carlo pipeline for the
STAR-calibrated hierarchical treatment-effect model:

1. generate a new STAR-like super-population in every replication and fit the
   nine full-sample scenarios;
2. perform leave-one-school-out (LOSO) estimation for ten schools in every
   scenario and replication; and
3. calculate oracle predictive-performance measures and average them over
   Monte Carlo replications within each scenario.

The observation model uses a normal likelihood, not a Student-t likelihood.


REQUIRED FILES
--------------

STAR_Data_Simulation.csv
    Empirical STAR data used to construct the joint 79-school distribution of
    (H1, H2, W, pF).

star_baseline_simulation.R
    Reads and validates the empirical STAR school population.

run_star_simulation_exercises.R
    Generates the master super-population and constructs the nine nested
    scenarios. Sourcing this file only loads functions; it does not begin the
    simulation automatically.

STAR_superpopulation_model.R
    Contains the Gaussian hierarchical Stan model, constructs Stan inputs,
    calculates HMC diagnostics, and provides safe checkpoint-saving utilities.
    This is a supporting file and is not normally run by itself.

STAR_superpopulation_estimation.R
    Stage 1. Generates each replication, estimates all nine full-sample
    scenarios, and saves the master population, scenario inputs, posterior
    draws, and diagnostics.

STAR_superpopulation_LOSO.R
    Stage 2. Performs ten separate one-school-out fits for each scenario and
    replication.

STAR_superpopulation_LOSO_metrics.R
    Stage 3. Evaluates LOSO predictions against the true population school ATE
    and averages performance measures over replications within each scenario.


SOFTWARE REQUIREMENTS
---------------------

The scripts require R and the rstan package. Install and check rstan before
starting the simulation:

    install.packages("rstan")
    library(rstan)
    rstan::stan_version()

Run the scripts from this directory:

    setwd("C:/ANDRES/HierarchicalMetaAnalysis/Code")


SUPER-POPULATION DESIGN
-----------------------

Every Monte Carlo replication generates a new master population containing:

    J_max = 80 schools
    N_max = 450 students per school
    total = 36,000 students

School covariates (H1, H2, W, pF) are sampled jointly, with replacement, from
the empirical distribution of the 79 STAR schools.

Within a replication, the same school covariates, student covariates,
treatment assignments, outcome errors, and standardized school heterogeneity
are shared across scenarios. Across replications, all these quantities are
regenerated.

Treatment assignment is:

    T_ij ~ Bernoulli(0.3).

The nine unique scenarios are:

    Scenario    J     N_j    c
    --------------------------------
    J_10        10    150    1
    J_20        20    150    1
    baseline    40    150    1
    J_80        80    150    1
    N_50        40     50    1
    N_450       40    450    1
    c_0         40    150    0
    c_0.5       40    150    0.5
    c_2         40    150    2

The baseline (J=40, N_j=150, c=1) is stored only once.

School samples are nested:

    J=10: schools 1,...,10
    J=20: schools 1,...,20
    J=40: schools 1,...,40
    J=80: schools 1,...,80

Student samples are also nested within school:

    N_j=50:  students 1,...,50
    N_j=150: students 1,...,150
    N_j=450: students 1,...,450

For residual heterogeneity,

    v_j(c) = sqrt(c) Sigma_0^(1/2) epsilon_j,
    epsilon_j ~ N(0,I_2).

The same epsilon_j is used for every value of c within a replication.


ESTIMATION MODEL
----------------

The Gaussian hierarchical model is:

    y_ij = X_ij' beta + W_j' gamma
           + T_ij Z_ij' theta_j + u_ij,

    u_ij ~ Normal(0, sigma^2),

    theta_j = Gamma' H_j + v_j,

    v_j ~ MVN(0, Sigma).

The regressors match the theoretical simulation DGP:

    X_ij = (1, female_ij, age_ij)
    Z_ij = (1, free_lunch_ij)
    W_j  = avg_experience_j
    H_j  = (1, H1_j, H2_j)

The model uses a non-centred parameterization for the school random effects, a
Cholesky representation of Sigma, half-normal priors on standard deviations,
and an LKJ prior on the random-effect correlation matrix.

The simulated outcome is already expressed in control-group standard-deviation
units. The estimation and LOSO scripts therefore do not rescale the outcome or
the STAR-calibrated covariates within a fold.


RECOMMENDED SMALL TEST
----------------------

Do not begin with all 100 replications. First edit the configuration near the
top of STAR_superpopulation_estimation.R to use:

    R_rep = 1L
    replications = 1L
    chains = 2L
    iter = 500L
    warmup = 250L
    refresh = 50L

Then run Stage 1:

    source("STAR_superpopulation_estimation.R")

For a small LOSO test, edit STAR_superpopulation_LOSO.R to use:

    replications = 1L
    scenarios = "baseline"
    chains = 2L
    iter = 500L
    warmup = 250L
    refresh = 50L

Then run Stage 2:

    source("STAR_superpopulation_LOSO.R")

Finally, edit STAR_superpopulation_LOSO_metrics.R consistently:

    replications = 1L
    scenarios = "baseline"

and run Stage 3:

    source("STAR_superpopulation_LOSO_metrics.R")

Inspect the saved HMC diagnostics before proceeding to production. In
particular, check divergences, maximum R-hat, effective sample sizes,
treedepth hits, and E-BFMI.


PRODUCTION EXECUTION
--------------------

After the small test succeeds, restore the production settings in
STAR_superpopulation_estimation.R:

    R_rep = 100L
    replications = NULL
    chains = 2L
    iter = 2000L
    warmup = 1000L

Run the scripts in this exact order:

    source("STAR_superpopulation_estimation.R")
    source("STAR_superpopulation_LOSO.R")
    source("STAR_superpopulation_LOSO_metrics.R")

Setting replications=NULL means that all available replications are processed.
Setting scenarios=NULL means that all nine scenarios are processed.

The full production design contains:

    100 x 9 = 900 full-sample fits,

and

    100 x 9 x 10 = 9,000 LOSO fits.

Consequently, the full experiment can require substantial computation time.
The estimation scripts save each completed fit separately and skip completed
files when restarted, unless overwrite=TRUE.

Independent jobs can process different replications by assigning disjoint
values to the replications configuration. Two jobs must never write the same
replication/scenario/fold simultaneously.


STAGE 1 OUTPUTS
---------------

By default, Stage 1 writes to:

    STAR_superpopulation_fits/

The directory structure is:

    STAR_superpopulation_fits/
      rep_0001/
        master_population.rds
        scenario_manifest.rds
        baseline/
          scenario_input.rds
          full_fit.rds
        J_10/
          scenario_input.rds
          full_fit.rds
        ...

master_population.rds
    Contains the single 80-school x 450-student realization for the
    replication. Scenario data are reconstructed from this object rather than
    stored repeatedly.

scenario_manifest.rds
    Contains the scenario table, replication information, population seed, and
    the ten common LOSO school indices.

scenario_input.rds
    Contains the scenario settings, true school effects, LOSO schools, master
    population path, and expected Stan dimensions.

full_fit.rds
    Contains compact posterior draws for beta, gamma, Gamma, Sigma, random-
    effect scales and correlations, and sigma, together with HMC diagnostics.

Set save_full_stanfit=TRUE only if complete stanfit objects are required. These
objects can consume substantial disk space.


STAGE 2: TEN-SCHOOL LOSO
-----------------------

The ten LOSO targets are the common nested schools 1,...,10. Each school is
omitted separately: this is ten leave-one-school-out fits, not one joint
ten-school holdout.

Using the same ten schools across scenarios permits direct counterfactual
comparisons. In the J=10 scenario, all ten schools are evaluated one at a time.

Each checkpoint is stored as:

    rep_XXXX/scenario/loso/loso_school_JJ.rds

The checkpoint contains:

    posterior draws of Gamma and Sigma;
    held-out H_j and population zbar_j;
    the true population and finite-sample school ATEs;
    HMC diagnostics;
    scenario and replication identifiers; and
    timing and seed information.

The primary predictive distribution integrates over a new realization of the
unobserved school random effect. It does not condition on the held-out school's
outcomes.


STAGE 3: PERFORMANCE MEASURES
-----------------------------

The primary oracle target is the population school ATE:

    tau_j,pop^true = (1,pF_j)' theta_j.

The finite-sample ATE based on the realized free-lunch proportion is also
stored, but it is not the primary target.

For every LOSO fold, the metrics script calculates:

    prediction error;
    absolute error;
    squared error;
    predictive standard deviation;
    50%, 80%, 90%, and 95% predictive intervals;
    interval coverage of the true population ATE; and
    predictive interval length.

Aggregation proceeds in two steps:

1. average over the ten LOSO schools within each replication and scenario;
2. average the replication-level measures over replications within each
   scenario.

This gives equal weight to every Monte Carlo replication.

The final files are:

    STAR_superpopulation_metrics/LOSO_fold_metrics.csv
    STAR_superpopulation_metrics/LOSO_replication_metrics.csv
    STAR_superpopulation_metrics/LOSO_scenario_average_metrics.csv

LOSO_fold_metrics.csv
    One row per replication, scenario, and held-out school.

LOSO_replication_metrics.csv
    Performance averaged over the ten held-out schools within each replication.

LOSO_scenario_average_metrics.csv
    Final scenario-level averages, standard deviations across replications,
    and Monte Carlo standard errors.


RESTARTING AND PARTIAL RUNS
---------------------------

Both estimation stages use checkpoint files. Existing completed fits are
skipped when overwrite=FALSE.

To rerun selected replications, set, for example:

    replications = 21:30

To run selected scenarios during LOSO, set, for example:

    scenarios = c("baseline", "J_80", "c_2")

Restore replications=NULL and scenarios=NULL when calculating final metrics
over the complete experiment.

Do not set overwrite=TRUE unless the intention is to replace valid completed
fits.


REPRODUCIBILITY
---------------

Every replication and model fit has a deterministic seed. A saved master
population is reused if Stage 1 is restarted, so completed and subsequent
scenario fits remain based on the same realized replication.

Within a replication:

    H, W, pF, X, Z, T, u, and epsilon are common across scenarios;
    school and student samples are nested;
    changing c rescales the same epsilon_j realization.

Across replications, a new independent master population is generated.


IMPORTANT INTERPRETATION
------------------------

Increasing N_j improves information about treatment effects in already
observed schools. Increasing J improves information about Gamma and Sigma.
Residual heterogeneity creates predictive uncertainty about the treatment
effect in a new school, even when its observed modifiers are known.

The common-random-numbers and nested-sample design ensures that scenario
comparisons are not dominated by unrelated population, treatment, or outcome
draws.
