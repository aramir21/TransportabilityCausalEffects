# Run from SimulationBHM, or source this file with chdir=TRUE from Code.
source("STAR_baseline_model_comparison_helpers.R")

run_baseline_model_comparison <- function(config = list()) {
  cfg <- modifyList(list(
    R_rep = 100L, replications = NULL, schools = 1:10,
    target = "population", csv_file = "STAR_Data_Simulation.csv",
    # Optional directory containing existing rep_XXXX/master_population.rds.
    population_dir = NULL, output_dir = "STAR_baseline_model_comparison",
    population_seed = 20260828L, fit_seed = 20261010L,
    chains = 4L, iter = 2000L, warmup = 1000L, cores = 1L,
    adapt_delta = .95, max_treedepth = 12L, refresh = 100L,
    interval_levels = c(.5, .8, .9, .95), overwrite = FALSE,
    fit = TRUE, score = TRUE
  ), config)
  cfg$target <- match.arg(cfg$target, c("population", "finite"))
  reps <- if (is.null(cfg$replications)) seq_len(cfg$R_rep) else cfg$replications
  stopifnot(length(reps) > 0L, all(reps == as.integer(reps)),
    all(reps >= 1L & reps <= cfg$R_rep), !anyDuplicated(reps),
    length(cfg$schools) > 0L, !anyDuplicated(cfg$schools),
    all(cfg$schools %in% 1:40), cfg$warmup >= 0L, cfg$warmup < cfg$iter,
    cfg$chains >= 1L, cfg$cores >= 1L,
    all(cfg$interval_levels > 0 & cfg$interval_levels < 1),
    !anyDuplicated(round(100 * cfg$interval_levels)))
  models <- c("M0", "M1", "M2", "M3")
  model_code <- list(M0 = star_normal_stan_code,
    aggregate = paste(readLines(file.path("..", "STAR_model_comparison.stan")), collapse = "\n"))
  compiled <- list()
  folds <- diagnostics <- list()
  settings <- star_scenario_table()
  settings <- settings[settings$scenario == "baseline", , drop = FALSE]
  destination <- file.path(cfg$output_dir, cfg$target)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  for (r in reps) {
    rep_name <- sprintf("rep_%04d", r)
    rep_dir <- file.path(destination, rep_name)
    # Cache a common master outside the target-specific output directory.
    master_path <- file.path(cfg$output_dir, "populations", rep_name, "master_population.rds")
    origin <- if (is.null(cfg$population_dir)) NULL else
      file.path(cfg$population_dir, rep_name, "master_population.rds")
    provenance <- list(origin = if (is.null(origin)) NULL else
      normalizePath(origin, mustWork = FALSE), seed = cfg$population_seed + r - 1L,
      generator = paste(readLines("run_star_simulation_exercises.R"), collapse = "\n"),
      empirical_schools = if (is.null(origin)) read_star_population(cfg$csv_file) else NULL)
    provenance_path <- paste0(master_path, ".provenance.rds")
    if (file.exists(master_path)) {
      if (!file.exists(provenance_path) || !identical(readRDS(provenance_path), provenance))
        stop("Population provenance changed. Use a new output_dir.")
      master <- readRDS(master_path)
      if (!is.null(origin) && !identical(master, readRDS(origin)))
        stop("Source master population changed: ", origin)
    } else {
      if (!cfg$fit) stop("Missing population cache: ", master_path)
      if (is.null(origin)) {
        set.seed(provenance$seed)
        master <- generate_star_master(provenance$empirical_schools)
      } else {
        if (!file.exists(origin)) stop("Missing existing master: ", origin)
        master <- readRDS(origin)
      }
      atomic_save_rds(master, master_path)
      atomic_save_rds(provenance, provenance_path)
    }
    scenario <- construct_star_scenario(master, 40L, 150L, 1)
    summaries <- baseline_school_summaries(scenario$data, cfg$target)
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summaries, file.path(rep_dir, "school_inputs.csv"), row.names = FALSE)
    for (j in cfg$schools) {
      # Held-out outcome summaries are used only for optional validation,
      # never in a training likelihood or preprocessing scale.
      micro <- make_star_stan_input(scenario$data, scenario$school_truth, j)
      h <- micro$metadata$h_heldout
      p <- summaries$p[summaries$school == j]
      truth_row <- scenario$school_truth[scenario$school_truth$sim_school_id == j, ]
      truth <- if (cfg$target == "population") truth_row$tau_pop_true else truth_row$tau_finite_true
      for (model in models) {
        path <- file.path(rep_dir, sprintf("%s_school_%02d.rds", model, j))
        stan_data <- if (model == "M0") micro$stan_data else
          baseline_aggregate_input(summaries, j, model)
        code_key <- if (model == "M0") "M0" else "aggregate"
        signature <- list(version = 1L, replication = r, school = j,
          model = model, target = cfg$target, stan_data = stan_data,
          h = h, p = p, model_code = model_code[[code_key]])
        if (file.exists(path) && !cfg$overwrite) {
          result <- readRDS(path)
          if (!identical(result$signature, signature))
            stop("Checkpoint specification changed: ", path, ". Use a new output_dir or overwrite.")
        } else {
          if (!cfg$fit) stop("Missing checkpoint: ", path)
          if (!requireNamespace("rstan", quietly = TRUE))
            stop("Install rstan and a working C++ toolchain to fit models.")
          if (is.null(compiled[[code_key]]))
            compiled[[code_key]] <- rstan::stan_model(model_code = model_code[[code_key]])
          seed <- cfg$fit_seed + 10000L * r + 100L * match(model, models) + j
          message("Baseline LOSO: rep ", r, ", ", model, ", school ", j)
          elapsed <- system.time(fit <- rstan::sampling(compiled[[code_key]],
            data = stan_data, chains = cfg$chains, iter = cfg$iter,
            warmup = cfg$warmup, cores = cfg$cores, seed = seed,
            refresh = cfg$refresh, control = list(adapt_delta = cfg$adapt_delta,
              max_treedepth = cfg$max_treedepth)))
          if (model == "M0") {
            posterior <- rstan::extract(fit,
              pars = c("Gamma", "Sigma_delta", "sd_delta", "Omega", "sigma"))
            prediction <- baseline_micro_prediction(posterior, h, p)
            locations <- prediction$conditional_locations
            sds <- prediction$conditional_sds
            diagnostic <- summarise_star_fit_diagnostics(fit, cfg$max_treedepth)
          } else {
            prediction <- rstan::extract(fit, pars = c("predictive_location", "predictive_sd"))
            locations <- prediction$predictive_location
            sds <- prediction$predictive_sd
            posterior <- rstan::extract(fit, pars = c("coefficients", "Sigma",
              if (model == "M3") c("sd_delta", "L_Omega") else "common_sd"))
            diagnostic <- comparison_diagnostics(fit, cfg$max_treedepth, model)
          }
          result <- list(signature = signature, config = cfg, seed = seed,
            settings = settings, posterior = posterior, locations = locations, sds = sds,
            diagnostics = diagnostic, elapsed_seconds = unname(elapsed["elapsed"]))
          atomic_save_rds(result, path)
          rm(fit)
          invisible(gc())
        }
        if (cfg$score) {
          row <- baseline_oracle_score(result$locations, result$sds, truth, cfg$interval_levels)
          folds[[length(folds) + 1L]] <- cbind(
            data.frame(replication = r, school = j, model = model, target = cfg$target), row)
          d <- result$diagnostics
          diagnostics[[length(diagnostics) + 1L]] <- data.frame(replication = r,
            school = j, model = model, max_rhat = d$max_rhat, min_n_eff = d$min_n_eff,
            divergences = d$divergences, treedepth_hits = d$treedepth_hits,
            min_ebfmi = d$min_ebfmi)
        }
      }
    }
  }
  if (!cfg$score) return(invisible(list(config = cfg, output_dir = destination)))
  folds <- do.call(rbind, folds)
  tables <- baseline_comparison_tables(folds, cfg$interval_levels)
  exports <- c(list(fold_metrics = folds, diagnostics = do.call(rbind, diagnostics)), tables)
  metrics_dir <- file.path(destination, "metrics")
  dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(exports))
    write.csv(exports[[name]], file.path(metrics_dir, paste0(name, ".csv")), row.names = FALSE)
  atomic_save_rds(list(config = cfg, settings = settings, replications = reps,
    schools = cfg$schools, model_code = model_code), file.path(metrics_dir, "metadata.rds"))
  message("Comparison metrics saved to ", normalizePath(metrics_dir))
  invisible(exports)
}

# Sourcing defines functions only. Call run_baseline_model_comparison() to run.
