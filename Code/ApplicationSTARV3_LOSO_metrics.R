# Latent-ATE predictive validation for the 79 leave-one-school-out STAR fits
#
# This compares p(tau_j | D_{-j}, H_j, Zbar_j), what was predicted before
# observing school j, with p(tau_j | D_{1:J}), what the full model learned after
# observing it. The direct ATE estimate and its SE are deliberately not used.

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

metrics_config <- list(
  loso_dir = "star_output/loso",
  full_fit_path = "star_output/STAR_kindergarten_reading_fit.RData",
  output_dir = "star_output/loso_latent_metrics",
  expected_folds = 79L,
  interval_levels = c(0.50, 0.80, 0.90, 0.95),
  predictive_replicates_per_posterior_draw = 4L,
  seed = 20260823L,
  save_predictive_draws = TRUE
)

if (metrics_config$predictive_replicates_per_posterior_draw < 1L) {
  stop("predictive_replicates_per_posterior_draw must be positive.")
}
if (any(metrics_config$interval_levels <= 0 |
        metrics_config$interval_levels >= 1)) {
  stop("All interval levels must lie strictly between zero and one.")
}
if (!file.exists(metrics_config$full_fit_path)) {
  stop("Full-sample fit not found: ", metrics_config$full_fit_path)
}

dir.create(metrics_config$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(metrics_config$seed)

# Recreate the estimation-sample ordering without fitting another model.
old_prepare_option <- getOption("star.prepare_only")
options(star.prepare_only = TRUE)
source("ApplicationSTARV3_clean.R", local = FALSE)
options(star.prepare_only = old_prepare_option)

# Load the full fit separately so it does not replace preprocessing objects.
full_fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(metrics_config$full_fit_path, envir = full_fit_environment)
if (!"fit" %in% loaded_names) {
  stop("The full-fit file does not contain an object named 'fit'.")
}
fit_full <- full_fit_environment$fit
full_delta <- rstan::extract(fit_full, pars = "delta")$delta

if (length(dim(full_delta)) != 3L || dim(full_delta)[2L] != J ||
    dim(full_delta)[3L] != ncol(Z)) {
  stop("The full posterior delta dimensions do not match the prepared sample.")
}

# Prefer the outcome scale saved with the fitted object and verify that it
# agrees with the current preprocessing.
if ("sample_description" %in% loaded_names &&
    !is.null(full_fit_environment$sample_description$control_sd)) {
  full_outcome_scale <- full_fit_environment$sample_description$control_sd
} else {
  full_outcome_scale <- control_sd
}
if (!is.finite(full_outcome_scale) || full_outcome_scale <= 0) {
  stop("The full-sample outcome scale is invalid.")
}
if (!isTRUE(all.equal(full_outcome_scale, control_sd, tolerance = 1e-10))) {
  stop(
    "The saved fit and current preprocessing use different outcome scales. ",
    "Use the clean-script configuration that generated the full fit."
  )
}

# -----------------------------------------------------------------------------
# Checkpoint discovery and utility functions
# -----------------------------------------------------------------------------

fold_files <- list.files(
  metrics_config$loso_dir,
  pattern = "^loso_school_[0-9]+\\.rds$",
  full.names = TRUE
)

extract_fold_number <- function(path) {
  as.integer(sub("^loso_school_([0-9]+)\\.rds$", "\\1", basename(path)))
}

fold_numbers <- unname(vapply(fold_files, extract_fold_number, integer(1)))
fold_files <- fold_files[order(fold_numbers)]
fold_numbers <- sort(fold_numbers)

if (length(fold_files) != metrics_config$expected_folds) {
  stop(
    "Expected ", metrics_config$expected_folds,
    " LOSO checkpoints but found ", length(fold_files), "."
  )
}
if (!identical(fold_numbers, seq_len(metrics_config$expected_folds))) {
  stop("LOSO checkpoint indices are duplicated or nonconsecutive.")
}
if (!identical(fold_numbers, sort(unique(analysis_data$stan_school)))) {
  stop("LOSO fold indices do not match the full estimation sample.")
}

interval_bounds <- function(draws, level) {
  alpha <- 1 - level
  bounds <- stats::quantile(
    draws,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    type = 8
  )
  c(lower = bounds[1L], upper = bounds[2L],
    length = bounds[2L] - bounds[1L])
}

source("STAR_latent_predictive_helpers.R")

full_school_ate_draws <- function(school_index, zbar) {
  S <- dim(full_delta)[1L]
  Kz <- dim(full_delta)[3L]
  delta_j <- matrix(
    full_delta[, school_index, , drop = FALSE],
    nrow = S,
    ncol = Kz
  )
  drop(delta_j %*% as.numeric(zbar))
}

# -----------------------------------------------------------------------------
# Compare pre-observation predictions with full-data latent-ATE posteriors
# -----------------------------------------------------------------------------

school_results <- vector("list", length(fold_files))
comparison_draws <- vector("list", length(fold_files))

for (j in seq_along(fold_files)) {
  fold <- readRDS(fold_files[j])
  school_index <- fold_numbers[j]
  if (!identical(as.integer(fold$metadata$omitted_school), school_index)) {
    stop("Checkpoint metadata does not match filename: ", fold_files[j])
  }

  loso_prediction <- latent_loso_ate_moments(fold)
  # Simulation is needed only for coverage and predictive intervals.
  predictive_locations <- rep(
    loso_prediction$conditional_locations,
    each = metrics_config$predictive_replicates_per_posterior_draw
  )
  predictive_sds <- rep(
    loso_prediction$conditional_sds,
    each = metrics_config$predictive_replicates_per_posterior_draw
  )
  loso_prediction$latent_draws <- stats::rnorm(
    length(predictive_locations), mean = predictive_locations, sd = predictive_sds
  )

  # LOSO outcomes use the 78-school control SD. Convert all LOSO effects to the
  # full-fit control-SD scale before making comparisons.
  fold_outcome_scale <- as.numeric(fold$metadata$outcome_scale)
  if (!is.finite(fold_outcome_scale) || fold_outcome_scale <= 0) {
    stop("Invalid fold outcome scale in school ", school_index, ".")
  }
  scale_conversion <- fold_outcome_scale / full_outcome_scale
  loso_draws <- loso_prediction$latent_draws * scale_conversion
  loso_conditional_locations <-
    loso_prediction$conditional_locations * scale_conversion

  zbar <- fold$metadata$zbar_heldout
  full_draws <- full_school_ate_draws(school_index, zbar)
  # Scaling both the component mean and SD includes the density Jacobian.
  log_predictive_density <- latent_log_predictive_density(
    full_draws,
    loso_conditional_locations,
    loso_prediction$conditional_sds * scale_conversion
  )

  loso_mean <- mean(loso_conditional_locations)
  full_mean <- mean(full_draws)
  full_variance <- stats::var(full_draws)
  mean_error <- loso_mean - full_mean

  row <- data.frame(
    fold = school_index,
    original_schoolid = as.character(fold$metadata$original_schoolid)[1L],
    heldout_N = fold$metadata$heldout_N,
    fold_outcome_scale = fold_outcome_scale,
    full_outcome_scale = full_outcome_scale,
    scale_conversion = scale_conversion,
    loso_predictive_mean = loso_mean,
    full_posterior_mean = full_mean,
    full_posterior_sd = stats::sd(full_draws),
    posterior_mean_error = mean_error,
    squared_posterior_mean_error = mean_error^2,
    posterior_expected_squared_error = mean_error^2 + full_variance,
    log_predictive_score = mean(log_predictive_density),
    log_score_full_draws = length(full_draws),
    log_score_loso_draws = length(loso_conditional_locations),
    stringsAsFactors = FALSE
  )

  for (level in metrics_config$interval_levels) {
    label <- sprintf("%02d", round(100 * level))
    loso_interval <- interval_bounds(loso_draws, level)

    row[[paste0("loso_lower_", label)]] <- loso_interval[["lower"]]
    row[[paste0("loso_upper_", label)]] <- loso_interval[["upper"]]
    row[[paste0("loso_length_", label)]] <- loso_interval[["length"]]

    # Soft coverage is the full posterior probability that the latent ATE lies
    # in the interval predicted before observing the school.
    row[[paste0("posterior_coverage_", label)]] <- mean(
      full_draws >= loso_interval[["lower"]] &
        full_draws <= loso_interval[["upper"]]
    )

    # Secondary binary summary based on the full posterior mean.
    row[[paste0("full_mean_covered_", label)]] <- as.integer(
      full_mean >= loso_interval[["lower"]] &&
        full_mean <= loso_interval[["upper"]]
    )
  }

  school_results[[j]] <- row
  comparison_draws[[j]] <- list(
    fold = school_index,
    original_schoolid = row$original_schoolid,
    loso_latent_ATE = loso_draws,
    full_posterior_latent_ATE = full_draws,
    loso_log_density_at_full_draws = log_predictive_density
  )
  message("Compared fold ", school_index, " of ", length(fold_files))
}

school_metrics <- do.call(rbind, school_results)
rownames(school_metrics) <- NULL

log_score_summary <- data.frame(
  metric = "Mean posterior-expected latent ATE log predictive score",
  value = mean(school_metrics$log_predictive_score),
  schools = nrow(school_metrics),
  weighting = "Equal weight per school",
  scale = "Full-fit control-group standard deviation",
  logarithm = "Natural log",
  interpretation = "Higher is better for a fixed full-data reference posterior"
)
utils::write.csv(
  log_score_summary,
  file.path(metrics_config$output_dir, "LOSO_latent_log_score_summary.csv"),
  row.names = FALSE
)

# Point-prediction RMSE relative to what was learned after observing each school.
rmse_posterior_means <- sqrt(mean(
  school_metrics$squared_posterior_mean_error
))

# Posterior expected RMSE integrates uncertainty remaining about each latent ATE:
# E[(m_j^- - tau_j)^2 | D] = (m_j^- - m_j^+)^2 + Var(tau_j | D).
rmse_posterior_expected <- sqrt(mean(
  school_metrics$posterior_expected_squared_error
))

coverage_summary <- do.call(rbind, lapply(
  metrics_config$interval_levels,
  function(level) {
    label <- sprintf("%02d", round(100 * level))
    posterior_coverage <- school_metrics[[paste0("posterior_coverage_", label)]]
    mean_covered <- school_metrics[[paste0("full_mean_covered_", label)]]
    interval_length <- school_metrics[[paste0("loso_length_", label)]]

    data.frame(
      nominal_coverage = level,
      average_posterior_coverage = mean(posterior_coverage),
      median_posterior_coverage = stats::median(posterior_coverage),
      full_posterior_mean_coverage = mean(mean_covered),
      average_interval_length = mean(interval_length),
      median_interval_length = stats::median(interval_length),
      schools = nrow(school_metrics)
    )
  }
))

rmse_summary <- data.frame(
  metric = c(
    "RMSE between LOSO and full posterior means",
    "Posterior expected RMSE for latent ATE"
  ),
  value = c(rmse_posterior_means, rmse_posterior_expected),
  schools = nrow(school_metrics),
  interpretation = c(
    "Point-prediction disagreement after observing each school",
    "Includes uncertainty remaining in the full posterior latent ATE"
  )
)

utils::write.csv(
  school_metrics,
  file.path(metrics_config$output_dir, "LOSO_latent_school_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  coverage_summary,
  file.path(metrics_config$output_dir, "LOSO_latent_coverage_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  rmse_summary,
  file.path(metrics_config$output_dir, "LOSO_latent_RMSE_summary.csv"),
  row.names = FALSE
)

if (isTRUE(metrics_config$save_predictive_draws)) {
  saveRDS(
    comparison_draws,
    file.path(metrics_config$output_dir, "LOSO_vs_full_latent_ATE_draws.rds"),
    compress = "xz"
  )
}

cat("\nLatent-ATE RMSE\n")
print(rmse_summary, row.names = FALSE)
cat("\nLatent-ATE log predictive score\n")
print(log_score_summary, row.names = FALSE)
cat("\nLatent-ATE coverage and LOSO interval lengths\n")
print(coverage_summary, row.names = FALSE)
message("Latent metrics written to: ", normalizePath(metrics_config$output_dir))
