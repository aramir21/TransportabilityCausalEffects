# One entry point: complete aggregate fits, then score M0--M3.
# Default: all 79 schools. Existing compatible checkpoints are reused.
# Run from the Code directory:
# source("ApplicationSTARV3_model_comparison_run.R")

run_star_model_comparison <- function(
    schools = seq_len(79L),
    loso_dir = "star_output/loso",
    fit_dir = "star_output/model_comparison",
    metrics_dir = NULL,
    full_metrics_path = "star_output/loso_latent_metrics/LOSO_latent_school_metrics.csv",
    chains = 4L, iter = 2000L, warmup = 1000L, cores = 1L,
    seed = 20260928L, adapt_delta = .95, max_treedepth = 12L,
    refresh = 0L) {
  stopifnot(length(schools) > 0L, !anyDuplicated(schools),
            all(schools %in% seq_len(79L)))
  schools <- sort(as.integer(schools))
  if (is.null(metrics_dir)) {
    metrics_dir <- if (identical(schools, seq_len(79L))) fit_dir else
      file.path(fit_dir, paste0("metrics_schools_", paste(schools, collapse = "_")))
  }
  if (!file.exists(full_metrics_path)) stop("Missing reporting-scale file: ", full_metrics_path)
  old_options <- options(c("star.comparison", "star.comparison.metrics"))
  on.exit(options(old_options), add = TRUE)
  options(star.comparison = list(
    schools = schools, models = c("M1", "M2", "M3"),
    loso_dir = loso_dir, output_dir = fit_dir, overwrite = FALSE,
    chains = chains, iter = iter, warmup = warmup, cores = cores,
    seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    refresh = refresh
  ))
  message("Step 1/2: fit M1--M3 for ", length(schools),
          " omitted schools; reuse existing checkpoints.")
  source("ApplicationSTARV3_model_comparison_fit.R", local = new.env(parent = globalenv()))

  options(star.comparison.metrics = list(
    schools = schools, models = c("M0", "M1", "M2", "M3"),
    loso_dir = loso_dir, fit_dir = fit_dir, output_dir = metrics_dir,
    full_metrics_path = full_metrics_path
  ))
  message("Step 2/2: calculate all four models' metrics for ", length(schools), " schools.")
  source("ApplicationSTARV3_model_comparison_metrics.R", local = new.env(parent = globalenv()))
  results <- read.csv(file.path(metrics_dir, "comparison_school_metrics.csv"))
  counts <- table(factor(results$model, levels = c("M0", "M1", "M2", "M3")))
  stopifnot(nrow(results) == 4L * length(schools), all(counts == length(schools)),
            !anyDuplicated(results[c("model", "fold")]),
            all(is.finite(results$log_predictive_score)))
  message("Complete: ", nrow(results), " model-school metric rows saved in ",
          normalizePath(metrics_dir))
  invisible(results)
}

# Old star.comparison school-1 settings do not restrict this runner.
# To request a subset explicitly, set options(star.comparison.run = list(schools = 1L)).
do.call(run_star_model_comparison, getOption("star.comparison.run", list()))
