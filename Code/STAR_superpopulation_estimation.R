# Stage 1: generate every super-population replication and estimate all nine
# full-sample scenarios under the Gaussian hierarchical model.

source("run_star_simulation_exercises.R")
source("STAR_superpopulation_model.R")

estimation_config <- list(
  R_rep = 100L,
  replications = NULL,       # NULL means 1:R_rep; e.g. 1:2 for a test
  seed = 20260828L,
  chains = 2L,
  iter = 2000L,
  warmup = 1000L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  refresh = 100L,
  cores = 4L,
  overwrite = FALSE,
  save_full_stanfit = FALSE,
  output_dir = "STAR_superpopulation_fits"
)

if (estimation_config$warmup >= estimation_config$iter) {
  stop("warmup must be smaller than iter.")
}
if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan first.")
dir.create(estimation_config$output_dir, recursive = TRUE, showWarnings = FALSE)
replications <- estimation_config$replications
if (is.null(replications)) replications <- seq_len(estimation_config$R_rep)
if (any(replications < 1L | replications > estimation_config$R_rep)) {
  stop("Invalid replication index.")
}

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, min(estimation_config$cores,
                               parallel::detectCores())))
message("Compiling the Gaussian hierarchical model once...")
compiled_model <- rstan::stan_model(model_code = star_normal_stan_code)
star_population <- read_star_population("STAR_Data_Simulation.csv")
scenario_table <- star_scenario_table()

for (r in replications) {
  rep_dir <- file.path(estimation_config$output_dir,
                       sprintf("rep_%04d", r))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  master_path <- file.path(rep_dir, "master_population.rds")

  if (file.exists(master_path)) {
    master <- readRDS(master_path)
  } else {
    set.seed(estimation_config$seed + r - 1L)
    master <- generate_star_master(star_population)
    atomic_save_rds(master, master_path)
  }

  # The ten common nested schools are the LOSO targets in every scenario.
  # J=10 contains exactly schools 1,...,10; randomization changes their order.
  set.seed(estimation_config$seed + 1000000L + r)
  common_loso_schools <- sample.int(10L, 10L, replace = FALSE)
  manifest <- list(replication = r, master_path = master_path,
                   scenario_table = scenario_table,
                   loso_schools = common_loso_schools,
                   population_seed = estimation_config$seed + r - 1L)
  atomic_save_rds(manifest, file.path(rep_dir, "scenario_manifest.rds"))

  for (s in seq_len(nrow(scenario_table))) {
    q <- scenario_table[s, ]
    scenario_dir <- file.path(rep_dir, q$scenario)
    dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(scenario_dir, "full_fit.rds")
    if (file.exists(result_path) && !estimation_config$overwrite) {
      message("Skipping completed full fit: replication ", r,
              ", scenario ", q$scenario)
      next
    }

    scenario <- construct_star_scenario(master, q$J, q$N_j, q$c_scale)
    model_input <- make_star_stan_input(scenario$data, scenario$school_truth)
    scenario_input <- list(
      replication = r, settings = q,
      school_truth = scenario$school_truth,
      loso_schools = common_loso_schools,
      # LOSO reconstructs data from master_path, avoiding duplicated storage.
      master_path = normalizePath(master_path, mustWork = TRUE),
      stan_dimensions = lapply(model_input$stan_data[c(
        "N", "J", "Kx", "Kz", "Kw", "Kh")], identity)
    )
    atomic_save_rds(scenario_input,
                    file.path(scenario_dir, "scenario_input.rds"))

    fit_seed <- estimation_config$seed + 100000L * r + s
    message("Full fit: replication ", r, ", scenario ", q$scenario)
    elapsed <- system.time({
      fit <- rstan::sampling(
        compiled_model, data = model_input$stan_data,
        chains = estimation_config$chains, iter = estimation_config$iter,
        warmup = estimation_config$warmup, seed = fit_seed,
        refresh = estimation_config$refresh,
        control = list(adapt_delta = estimation_config$adapt_delta,
                       max_treedepth = estimation_config$max_treedepth)
      )
    })
    posterior <- rstan::extract(
      fit, pars = c("beta", "gamma", "Gamma", "sd_delta", "Omega",
                    "Sigma_delta", "sigma"))
    fit_result <- list(
      replication = r, settings = q, posterior = posterior,
      diagnostics = summarise_star_fit_diagnostics(
        fit, estimation_config$max_treedepth),
      elapsed_seconds = unname(elapsed[["elapsed"]]), seed = fit_seed,
      model_config = estimation_config,
      column_names = model_input$metadata$column_names
    )
    atomic_save_rds(fit_result, result_path)
    if (estimation_config$save_full_stanfit) {
      atomic_save_rds(fit, file.path(scenario_dir, "full_stanfit.rds"))
    }
    rm(fit, posterior, fit_result, scenario, model_input)
    invisible(gc())
  }
}

message("Full-scenario estimation complete: ",
        normalizePath(estimation_config$output_dir))
