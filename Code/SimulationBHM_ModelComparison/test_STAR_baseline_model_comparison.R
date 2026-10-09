# Base-R checks; no RStan or sampling required. Run from SimulationBHM_ModelComparison.
source("STAR_baseline_model_comparison.R")
set.seed(81)
master <- generate_star_master(read_star_population())
scenario <- construct_star_scenario(master, 40L, 150L, 1)
pop <- baseline_school_summaries(scenario$data, "population")
finite <- baseline_school_summaries(scenario$data, "finite")
stopifnot(nrow(pop) == 40L, identical(pop$p, master$schools$pF[1:40]),
  isTRUE(all.equal(finite$p, scenario$school_truth$sample_pF)))
for (model in c("M1", "M2", "M3")) {
  a <- baseline_aggregate_input(pop, 1L, model)
  altered <- pop
  altered$estimate[1L] <- 1e6
  altered$sampling_se[1L] <- 1e6
  stopifnot(identical(a, baseline_aggregate_input(altered, 1L, model)), a$N == 39L)
}
changed <- scenario$data
changed$y[changed$sim_school_id == 1L] <- 1e6
stopifnot(identical(
  make_star_stan_input(scenario$data, scenario$school_truth, 1L)$stan_data,
  make_star_stan_input(changed, scenario$school_truth, 1L)$stan_data))

# M2/M3 reproduce the DGP mean; coefficient order is vec(Gamma).
a <- baseline_aggregate_input(pop, 1L, "M3")
Gamma <- t(master$parameters$Gamma_0)
stopifnot(isTRUE(all.equal(as.numeric(a$X %*% as.vector(Gamma)),
  scenario$school_truth$mu_systematic[-1L], tolerance = 1e-12)))
Sigma <- matrix(c(.2, -.07, -.07, .3), 2L)
L <- t(chol(Sigma))
p <- pop$p
quadratic <- Sigma[1,1] + 2*p*Sigma[1,2] + p^2*Sigma[2,2]
stable <- (L[1,1] + p*L[2,1])^2 + (p*L[2,2])^2
stopifnot(isTRUE(all.equal(quadratic, stable)))
posterior <- list(Gamma = array(rep(Gamma, 3L), c(3L, 3L, 2L)),
                  Sigma_delta = array(0, c(3L, 2L, 2L)))
# Assign explicitly to preserve draw/row/column ordering.
for (s in 1:3) { posterior$Gamma[s,,] <- Gamma; posterior$Sigma_delta[s,,] <- Sigma }
moments <- baseline_micro_prediction(posterior, c(1, pop$H1[1L], pop$H2[1L]), pop$p[1L])
stopifnot(max(abs(moments$conditional_locations - scenario$school_truth$mu_systematic[1L])) < 1e-12,
          max(abs(moments$conditional_sds^2 - quadratic[1L])) < 1e-12)

row <- baseline_oracle_score(rep(.3, 4), rep(.7, 4), .1, c(.5, .95))
stopifnot(abs(row$log_predictive_score - dnorm(.1, .3, .7, log = TRUE)) < 1e-12,
          abs(row$predictive_sd - .7) < 1e-12,
          abs(row$lower_95 - qnorm(.025, .3, .7)) < 1e-9)
scaled <- baseline_oracle_score(rep(.6, 4), rep(1.4, 4), .2, c(.5, .95))
stopifnot(abs(scaled$log_predictive_score - row$log_predictive_score + log(2)) < 1e-12)

# Paired MCSE is computed from replication means, with matched folds.
folds <- expand.grid(replication = 1:3, school = 1:2, model = c("M0", "M1", "M2", "M3"))
folds$error <- .1
folds$squared_error <- .01
folds$predictive_sd <- .7
folds$covered_95 <- 1
folds$length_95 <- 2
folds$log_predictive_score <- ifelse(folds$model == "M0", folds$replication, 0)
tables <- baseline_comparison_tables(folds, .95)
pair <- tables$paired_summary[tables$paired_summary$model == "M0" &
  tables$paired_summary$comparator == "M1", ]
stopifnot(nrow(pair) == 1L, abs(pair$mean_log_score_difference - 2) < 1e-12,
          abs(pair$monte_carlo_se - 1/sqrt(3)) < 1e-12)

# End-to-end scoring-only run using synthetic posterior checkpoints.
temporary <- tempfile("baseline_comparison_")
dir.create(temporary)
tryCatch({
  output <- file.path(temporary, "comparison")
  master_path <- file.path(output, "populations", "rep_0001", "master_population.rds")
  atomic_save_rds(master, master_path)
  provenance <- list(origin = NULL, seed = 20260828L,
    generator = paste(readLines("run_star_simulation_exercises.R"), collapse = "\n"),
    empirical_schools = read_star_population())
  atomic_save_rds(provenance, paste0(master_path, ".provenance.rds"))
  code <- list(M0 = star_normal_stan_code,
    aggregate = paste(readLines("STAR_model_comparison.stan"), collapse = "\n"))
  micro <- make_star_stan_input(scenario$data, scenario$school_truth, 1L)
  for (model in c("M0", "M1", "M2", "M3")) {
    data <- if (model == "M0") micro$stan_data else baseline_aggregate_input(pop, 1L, model)
    signature <- list(version = 1L, replication = 1L, school = 1L, model = model,
      target = "population", stan_data = data, h = micro$metadata$h_heldout,
      p = pop$p[1L], model_code = code[[if (model == "M0") "M0" else "aggregate"]])
    checkpoint <- list(signature = signature, locations = rep(.2, 4), sds = rep(.5, 4),
      diagnostics = list(max_rhat = 1, min_n_eff = 100, divergences = 0,
        treedepth_hits = 0, min_ebfmi = 1))
    atomic_save_rds(checkpoint, file.path(output, "population", "rep_0001",
      sprintf("%s_school_01.rds", model)))
  }
  result <- run_baseline_model_comparison(list(R_rep = 1L, replications = 1L,
    schools = 1L, output_dir = output, fit = FALSE))
  stopifnot(nrow(result$fold_metrics) == 4L, nrow(result$paired_summary) == 6L)
  # A stale training likelihood must be rejected on resume.
  path <- file.path(output, "population", "rep_0001", "M1_school_01.rds")
  stale <- readRDS(path)
  stale$signature$stan_data$estimate[1L] <- 999
  saveRDS(stale, path)
  failure <- try(run_baseline_model_comparison(list(R_rep = 1L, replications = 1L,
    schools = 1L, output_dir = output, fit = FALSE)), silent = TRUE)
  stopifnot(inherits(failure, "try-error"))
}, finally = unlink(temporary, recursive = TRUE))
cat("Baseline comparison checks passed.\n")
