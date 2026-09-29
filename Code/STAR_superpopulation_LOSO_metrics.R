# Stage 3: oracle LOSO metrics, averaged over Monte Carlo replications within
# each scenario. Predictions target tau_j,pop^true using zbar=(1,pF_j).

metrics_config <- list(
  input_dir = "STAR_superpopulation_fits",
  output_dir = "STAR_superpopulation_metrics",
  replications = NULL,
  scenarios = NULL,
  expected_loso_schools = 10L,
  interval_levels = c(0.50, 0.80, 0.90, 0.95),
  predictive_replicates_per_draw = 4L,
  seed = 20260830L,
  save_predictive_draws = FALSE
)

if (any(metrics_config$interval_levels <= 0 |
        metrics_config$interval_levels >= 1)) {
  stop("Interval levels must lie strictly between zero and one.")
}
if (metrics_config$predictive_replicates_per_draw < 1L) {
  stop("predictive_replicates_per_draw must be positive.")
}
dir.create(metrics_config$output_dir, recursive = TRUE, showWarnings = FALSE)

interval_bounds <- function(draws, level) {
  alpha <- 1 - level
  q <- stats::quantile(draws, c(alpha / 2, 1 - alpha / 2),
                       names = FALSE, type = 8)
  c(lower = q[1L], upper = q[2L], length = q[2L] - q[1L])
}

simulate_loso_population_ate <- function(fold, replicates_per_draw) {
  Gamma <- fold$posterior$Gamma
  Sigma <- fold$posterior$Sigma_delta
  h <- as.numeric(fold$metadata$h_heldout)
  zbar <- as.numeric(fold$metadata$zbar_heldout_pop)
  dG <- dim(Gamma)
  dS <- dim(Sigma)
  if (length(dG) != 3L || length(dS) != 3L) {
    stop("Invalid Gamma or Sigma posterior dimensions.")
  }
  S <- dG[1L]
  if (length(h) != dG[2L] || length(zbar) != dG[3L] ||
      !identical(dS, c(S, dG[3L], dG[3L]))) {
    stop("Posterior dimensions disagree with held-out covariates.")
  }

  location <- variance <- numeric(S)
  for (s in seq_len(S)) {
    Gamma_s <- matrix(Gamma[s, , ], nrow = dG[2L], ncol = dG[3L])
    Sigma_s <- matrix(Sigma[s, , ], nrow = dG[3L], ncol = dG[3L])
    location[s] <- drop(crossprod(zbar, t(Gamma_s) %*% h))
    variance[s] <- drop(crossprod(zbar, Sigma_s %*% zbar))
  }
  if (any(!is.finite(variance)) || any(variance < -1e-10)) {
    stop("Negative or nonfinite predictive variance.")
  }
  variance <- pmax(variance, 0)
  draws <- stats::rnorm(
    S * replicates_per_draw,
    mean = rep(location, each = replicates_per_draw),
    sd = sqrt(rep(variance, each = replicates_per_draw)))
  list(location = location, draws = draws)
}

rep_dirs <- list.dirs(metrics_config$input_dir, recursive = FALSE,
                      full.names = TRUE)
rep_dirs <- rep_dirs[grepl("^rep_[0-9]{4}$", basename(rep_dirs))]
rep_number <- as.integer(sub("^rep_", "", basename(rep_dirs)))
if (!is.null(metrics_config$replications)) {
  keep <- rep_number %in% as.integer(metrics_config$replications)
  rep_dirs <- rep_dirs[keep]
  rep_number <- rep_number[keep]
}
ordering <- order(rep_number)
rep_dirs <- rep_dirs[ordering]
rep_number <- rep_number[ordering]
if (!length(rep_dirs)) stop("No replication directories found.")

fold_rows <- list()
saved_draws <- list()
counter <- 0L
for (rr in seq_along(rep_dirs)) {
  r <- rep_number[rr]
  manifest <- readRDS(file.path(rep_dirs[rr], "scenario_manifest.rds"))
  settings <- manifest$scenario_table
  if (!is.null(metrics_config$scenarios)) {
    settings <- settings[settings$scenario %in% metrics_config$scenarios, ]
  }

  for (s in seq_len(nrow(settings))) {
    scenario_name <- settings$scenario[s]
    loso_dir <- file.path(rep_dirs[rr], scenario_name, "loso")
    files <- list.files(loso_dir, "^loso_school_[0-9]+\\.rds$",
                        full.names = TRUE)
    school_number <- as.integer(sub(
      "^loso_school_([0-9]+)\\.rds$", "\\1", basename(files)))
    files <- files[order(school_number)]
    school_number <- sort(school_number)
    expected <- sort(as.integer(manifest$loso_schools))
    if (length(files) != metrics_config$expected_loso_schools ||
        !identical(school_number, expected)) {
      stop("Incomplete LOSO checkpoints for replication ", r,
           ", scenario ", scenario_name, ".")
    }

    for (f in seq_along(files)) {
      fold <- readRDS(files[f])
      set.seed(metrics_config$seed + 100000L * r + 100L * s + f)
      prediction <- simulate_loso_population_ate(
        fold, metrics_config$predictive_replicates_per_draw)
      truth <- as.numeric(fold$metadata$tau_pop_true)
      prediction_mean <- mean(prediction$location)
      error <- prediction_mean - truth
      row <- data.frame(
        replication = r, scenario = scenario_name,
        experiment = settings$experiment[s], J = settings$J[s],
        N_j = settings$N_j[s], c_scale = settings$c_scale[s],
        school = school_number[f], truth = truth,
        predictive_mean = prediction_mean, error = error,
        absolute_error = abs(error), squared_error = error^2,
        predictive_sd = stats::sd(prediction$draws),
        stringsAsFactors = FALSE
      )
      for (level in metrics_config$interval_levels) {
        label <- sprintf("%02d", round(100 * level))
        interval <- interval_bounds(prediction$draws, level)
        row[[paste0("lower_", label)]] <- interval[["lower"]]
        row[[paste0("upper_", label)]] <- interval[["upper"]]
        row[[paste0("length_", label)]] <- interval[["length"]]
        row[[paste0("covered_", label)]] <- as.integer(
          truth >= interval[["lower"]] && truth <= interval[["upper"]])
      }
      counter <- counter + 1L
      fold_rows[[counter]] <- row
      if (metrics_config$save_predictive_draws) {
        saved_draws[[counter]] <- list(
          replication = r, scenario = scenario_name,
          school = school_number[f], truth = truth,
          predictive_draws = prediction$draws)
      }
    }
    message("Metrics: replication ", r, ", scenario ", scenario_name)
  }
}

fold_metrics <- do.call(rbind, fold_rows)
rownames(fold_metrics) <- NULL

# Equal-weight aggregation: first across ten LOSO schools within a replication,
# then across replications within a scenario.
group_keys <- unique(fold_metrics[c("replication", "scenario", "experiment",
                                    "J", "N_j", "c_scale")])
rep_rows <- vector("list", nrow(group_keys))
for (g in seq_len(nrow(group_keys))) {
  key <- group_keys[g, ]
  x <- fold_metrics[
    fold_metrics$replication == key$replication &
      fold_metrics$scenario == key$scenario, ]
  base <- data.frame(key, schools = nrow(x), bias = mean(x$error),
                     MAE = mean(x$absolute_error),
                     RMSE = sqrt(mean(x$squared_error)),
                     mean_predictive_sd = mean(x$predictive_sd))
  for (level in metrics_config$interval_levels) {
    label <- sprintf("%02d", round(100 * level))
    base[[paste0("coverage_", label)]] <- mean(x[[paste0("covered_", label)]])
    base[[paste0("mean_length_", label)]] <- mean(x[[paste0("length_", label)]])
  }
  rep_rows[[g]] <- base
}
replication_metrics <- do.call(rbind, rep_rows)
rownames(replication_metrics) <- NULL

metric_names <- setdiff(names(replication_metrics),
                        c("replication", "scenario", "experiment", "J",
                          "N_j", "c_scale", "schools"))
scenario_rows <- list()
counter <- 0L
for (scenario_name in unique(replication_metrics$scenario)) {
  x <- replication_metrics[replication_metrics$scenario == scenario_name, ]
  for (metric in metric_names) {
    values <- x[[metric]]
    counter <- counter + 1L
    scenario_rows[[counter]] <- data.frame(
      scenario = scenario_name, experiment = x$experiment[1L],
      J = x$J[1L], N_j = x$N_j[1L], c_scale = x$c_scale[1L],
      metric = metric, mean_over_replications = mean(values),
      sd_across_replications = stats::sd(values),
      monte_carlo_se = stats::sd(values) / sqrt(length(values)),
      replications = length(values), stringsAsFactors = FALSE)
  }
}
scenario_metrics <- do.call(rbind, scenario_rows)
rownames(scenario_metrics) <- NULL

utils::write.csv(fold_metrics,
                 file.path(metrics_config$output_dir, "LOSO_fold_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(
  replication_metrics,
  file.path(metrics_config$output_dir, "LOSO_replication_metrics.csv"),
  row.names = FALSE)
utils::write.csv(
  scenario_metrics,
  file.path(metrics_config$output_dir, "LOSO_scenario_average_metrics.csv"),
  row.names = FALSE)
if (metrics_config$save_predictive_draws) {
  saveRDS(saved_draws,
          file.path(metrics_config$output_dir, "LOSO_predictive_draws.rds"),
          compress = "xz")
}

print(scenario_metrics, row.names = FALSE)
message("Metrics written to: ", normalizePath(metrics_config$output_dir))
