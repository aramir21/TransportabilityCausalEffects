# Stage 2: ten leave-one-school-out fits for each scenario and replication.
# Run STAR_superpopulation_estimation.R first to create the manifests/inputs.

source("run_star_simulation_exercises.R")
source("STAR_superpopulation_model.R")

loso_config <- list(
  input_dir = "STAR_superpopulation_fits",
  replications = NULL,       # NULL discovers every completed replication
  scenarios = NULL,          # NULL uses all nine scenarios
  chains = 2L,
  iter = 2000L,
  warmup = 1000L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  refresh = 100L,
  seed = 20260829L,
  cores = 4L,
  overwrite = FALSE,
  save_full_stanfit = FALSE
)

if (loso_config$warmup >= loso_config$iter) stop("warmup must be below iter.")
if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan first.")

rep_dirs <- list.dirs(loso_config$input_dir, recursive = FALSE,
                      full.names = TRUE)
rep_dirs <- rep_dirs[grepl("^rep_[0-9]{4}$", basename(rep_dirs))]
if (!length(rep_dirs)) stop("No replication directories found.")
rep_number <- as.integer(sub("^rep_", "", basename(rep_dirs)))
if (!is.null(loso_config$replications)) {
  keep <- rep_number %in% as.integer(loso_config$replications)
  rep_dirs <- rep_dirs[keep]
  rep_number <- rep_number[keep]
}
ordering <- order(rep_number)
rep_dirs <- rep_dirs[ordering]
rep_number <- rep_number[ordering]

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, min(loso_config$cores, parallel::detectCores())))
message("Compiling the Gaussian hierarchical model once for all LOSO fits...")
compiled_model <- rstan::stan_model(model_code = star_normal_stan_code)

for (rr in seq_along(rep_dirs)) {
  r <- rep_number[rr]
  rep_dir <- rep_dirs[rr]
  manifest_path <- file.path(rep_dir, "scenario_manifest.rds")
  if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)
  manifest <- readRDS(manifest_path)
  master <- readRDS(manifest$master_path)
  settings <- manifest$scenario_table
  if (!is.null(loso_config$scenarios)) {
    settings <- settings[settings$scenario %in% loso_config$scenarios, ]
  }
  if (!nrow(settings)) stop("No requested scenarios found.")

  for (s in seq_len(nrow(settings))) {
    q <- settings[s, ]
    scenario_dir <- file.path(rep_dir, q$scenario)
    input_path <- file.path(scenario_dir, "scenario_input.rds")
    if (!file.exists(input_path)) stop("Missing scenario input: ", input_path)
    scenario_input <- readRDS(input_path)
    scenario <- construct_star_scenario(master, q$J, q$N_j, q$c_scale)

    loso_dir <- file.path(scenario_dir, "loso")
    dir.create(loso_dir, recursive = TRUE, showWarnings = FALSE)
    for (omitted_school in scenario_input$loso_schools) {
      result_path <- file.path(
        loso_dir, sprintf("loso_school_%02d.rds", omitted_school))
      if (file.exists(result_path) && !loso_config$overwrite) {
        message("Skipping completed LOSO: rep ", r, ", ", q$scenario,
                ", school ", omitted_school)
        next
      }

      fold <- make_star_stan_input(
        scenario$data, scenario$school_truth, omitted_school)
      fit_seed <- loso_config$seed + 100000L * r +
        100L * match(q$scenario, manifest$scenario_table$scenario) +
        omitted_school
      message("LOSO fit: rep ", r, ", scenario ", q$scenario,
              ", school ", omitted_school)
      elapsed <- system.time({
        fit <- rstan::sampling(
          compiled_model, data = fold$stan_data,
          chains = loso_config$chains, iter = loso_config$iter,
          warmup = loso_config$warmup, seed = fit_seed,
          refresh = loso_config$refresh,
          control = list(adapt_delta = loso_config$adapt_delta,
                         max_treedepth = loso_config$max_treedepth)
        )
      })
      posterior <- rstan::extract(
        fit, pars = c("Gamma", "Sigma_delta", "sd_delta", "Omega", "sigma"))
      fold_result <- list(
        replication = r, scenario = q$scenario, settings = q,
        metadata = fold$metadata, posterior = posterior,
        diagnostics = summarise_star_fit_diagnostics(
          fit, loso_config$max_treedepth),
        elapsed_seconds = unname(elapsed[["elapsed"]]), seed = fit_seed,
        model_config = loso_config
      )
      atomic_save_rds(fold_result, result_path)
      if (loso_config$save_full_stanfit) {
        atomic_save_rds(
          fit, file.path(loso_dir,
                         sprintf("loso_school_%02d_stanfit.rds",
                                 omitted_school)))
      }
      rm(fit, posterior, fold_result, fold)
      invisible(gc())
    }
  }
  rm(master)
  invisible(gc())
}

message("All requested super-population LOSO fits are complete.")
