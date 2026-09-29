# Fit M1--M3; reuse existing M0 LOSO fits. No microdata packages are needed.
source("STAR_model_comparison_helpers.R")
comparison_config <- modifyList(list(
  loso_dir = "star_output/loso", output_dir = "star_output/model_comparison",
  schools = seq_len(79L), models = c("M1", "M2", "M3"),
  # These 78-row models are small; sequential chains avoid Windows worker
  # startup costs while retaining four independent chains per fit.
  chains = 4L, iter = 2000L, warmup = 1000L, cores = 1L,
  seed = 20260928L, adapt_delta = .95, max_treedepth = 12L,
  refresh = 100L, overwrite = FALSE
), getOption("star.comparison", list()))
stopifnot(comparison_config$warmup < comparison_config$iter,
          all(comparison_config$schools %in% seq_len(79L)),
          all(comparison_config$models %in% c("M1", "M2", "M3")))
inputs <- read_comparison_inputs(comparison_config$loso_dir)
dir.create(comparison_config$output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(inputs$schools, file.path(comparison_config$output_dir, "school_inputs.csv"),
          row.names = FALSE)
compiled <- NULL
for (model in comparison_config$models) {
  for (j in comparison_config$schools) {
    data_j <- aggregate_fold_data(inputs, j, model)
    path <- file.path(comparison_config$output_dir, sprintf("%s_school_%02d.rds", model, j))
    signature <- list(model = model, stan_data = data_j,
                      model_code = readLines("STAR_model_comparison.stan"),
                      original_schoolid = inputs$schools$original_schoolid[j])
    if (file.exists(path) && !comparison_config$overwrite) {
      old <- readRDS(path)
      if (!identical(old$signature, signature))
        stop("Checkpoint specification changed: ", path, ". Use a new output_dir or overwrite.")
      message("Skipping ", model, " school ", j)
      next
    }
    if (is.null(compiled)) {
      if (!requireNamespace("rstan", quietly = TRUE))
        stop("Install rstan and its C++ toolchain to fit aggregate models.")
      compiled <- rstan::stan_model(file = "STAR_model_comparison.stan")
    }
    message("Fitting ", model, " school ", j)
    fit <- rstan::sampling(compiled, data = data_j,
      chains = comparison_config$chains, iter = comparison_config$iter,
      warmup = comparison_config$warmup, cores = comparison_config$cores,
      seed = comparison_config$seed + 1000L * match(model, c("M1", "M2", "M3")) + j,
      refresh = comparison_config$refresh,
      control = list(adapt_delta = comparison_config$adapt_delta,
                     max_treedepth = comparison_config$max_treedepth))
    draws <- rstan::extract(fit, pars = c("predictive_location", "predictive_sd"))
    parameter_draws <- rstan::extract(fit, pars = c("coefficients", "Sigma",
      if (model == "M3") c("sd_delta", "L_Omega") else "common_sd"))
    result <- list(signature = signature, config = comparison_config,
      fold = j, model = model, outcome_scale = inputs$metadata[[j]]$outcome_scale,
      locations = draws$predictive_location, sds = draws$predictive_sd,
      posterior = parameter_draws,
      diagnostics = comparison_diagnostics(fit, comparison_config$max_treedepth, model))
    temporary <- tempfile(tmpdir = comparison_config$output_dir, fileext = ".rds")
    saveRDS(result, temporary)
    if (file.exists(path) && !file.remove(path)) stop("Cannot replace ", path)
    if (!file.rename(temporary, path)) stop("Cannot save ", path)
    message("Saved ", model, " school ", j, "; divergences = ",
            result$diagnostics$divergences, "; max R-hat = ", result$diagnostics$max_rhat)
    rm(fit)
    invisible(gc())
  }
}
