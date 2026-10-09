# Baseline M0--M3 LOSO comparison. Source from SimulationBHM.
source("run_star_simulation_exercises.R")
source("STAR_superpopulation_model.R")
source(file.path("..", "STAR_model_comparison_helpers.R"), chdir = TRUE)

baseline_school_summaries <- function(data, target = c("population", "finite")) {
  target <- match.arg(target)
  rows <- lapply(sort(unique(data$sim_school_id)), function(j) {
    d <- data[data$sim_school_id == j, , drop = FALSE]
    treated <- d$y[d$treatment == 1L]
    control <- d$y[d$treatment == 0L]
    if (min(length(treated), length(control)) < 2L)
      stop("Need at least two observations in each treatment arm: school ", j)
    if (any(!is.finite(as.matrix(d[c("y", "H1", "H2", "pF")]))))
      stop("Nonfinite school data.")
    for (v in c("H1", "H2", "pF"))
      if (length(unique(d[[v]])) != 1L) stop("School covariate varies: ", v)
    se <- sqrt(var(treated) / length(treated) + var(control) / length(control))
    if (!is.finite(se) || se <= 0) stop("Invalid sampling SE.")
    data.frame(school = j, H1 = d$H1[1L], H2 = d$H2[1L],
      p = if (target == "population") d$pF[1L] else mean(d$free_lunch),
      estimate = mean(treated) - mean(control), sampling_se = se)
  })
  do.call(rbind, rows)
}

baseline_aggregate_input <- function(summaries, j, model) {
  stopifnot(model %in% c("M1", "M2", "M3"),
            sum(summaries$school == j) == 1L,
            all(is.finite(as.matrix(summaries))),
            all(summaries$p >= 0 & summaries$p <= 1))
  h <- cbind(intercept = 1, H1 = summaries$H1, H2 = summaries$H2)
  x <- if (model == "M1") cbind(h, p = summaries$p) else
    cbind(h, h * summaries$p)
  train <- summaries$school != j
  list(N = sum(train), K = ncol(x), composition_variance = as.integer(model == "M3"),
    X = x[train, , drop = FALSE], estimate = summaries$estimate[train],
    sampling_se = summaries$sampling_se[train], p = summaries$p[train],
    x_new = as.numeric(x[!train, ]), p_new = summaries$p[!train])
}

baseline_micro_prediction <- function(posterior, h, p) {
  latent_loso_ate_moments(list(posterior = posterior,
    metadata = list(h_heldout = h, zbar_heldout = c(1, p))))
}

baseline_oracle_score <- function(locations, sds, truth, levels) {
  stopifnot(length(truth) == 1L, is.finite(truth),
    length(locations) > 0L, length(locations) == length(sds),
    all(is.finite(c(locations, sds))), all(sds > 0),
    all(levels > 0 & levels < 1))
  error <- mean(locations) - truth
  row <- data.frame(truth = truth, predictive_mean = mean(locations),
    error = error, squared_error = error^2,
    predictive_sd = sqrt(mean(sds^2) + mean((locations - mean(locations))^2)),
    log_predictive_score = latent_log_predictive_density(truth, locations, sds))
  for (level in levels) {
    label <- sprintf("%02d", round(100 * level))
    lower <- mixture_quantile((1 - level) / 2, locations, sds)
    upper <- mixture_quantile((1 + level) / 2, locations, sds)
    row[[paste0("lower_", label)]] <- lower
    row[[paste0("upper_", label)]] <- upper
    row[[paste0("length_", label)]] <- upper - lower
    row[[paste0("covered_", label)]] <- as.integer(truth >= lower && truth <= upper)
  }
  row
}

baseline_comparison_tables <- function(folds, levels) {
  # Independent Monte Carlo replications, not overlapping LOSO folds,
  # are the units used for Monte Carlo standard errors.
  rows <- lapply(split(folds, interaction(folds$replication, folds$model, drop = TRUE)),
    function(x) {
      row <- data.frame(replication = x$replication[1L], model = x$model[1L],
        schools = nrow(x), bias = mean(x$error),
        RMSE = sqrt(mean(x$squared_error)), MSE = mean(x$squared_error),
        log_predictive_score = mean(x$log_predictive_score),
        predictive_sd = mean(x$predictive_sd))
      for (level in levels) {
        label <- sprintf("%02d", round(100 * level))
        row[[paste0("coverage_", label)]] <- mean(x[[paste0("covered_", label)]])
        row[[paste0("length_", label)]] <- mean(x[[paste0("length_", label)]])
      }
      row
    })
  reps <- do.call(rbind, rows)
  metrics <- setdiff(names(reps), c("replication", "model", "schools"))
  summary <- do.call(rbind, lapply(split(reps, reps$model), function(x)
    do.call(rbind, lapply(metrics, function(metric) {
      values <- x[[metric]]
      data.frame(model = x$model[1L], metric = metric, mean = mean(values),
        sd = sd(values), monte_carlo_se = sd(values) / sqrt(length(values)),
        replications = length(values))
    }))))
  pairs <- combn(unique(folds$model), 2L, simplify = FALSE)
  paired_reps <- do.call(rbind, lapply(pairs, function(pair) {
    a <- folds[folds$model == pair[1L], c("replication", "school", "log_predictive_score")]
    b <- folds[folds$model == pair[2L], c("replication", "school", "log_predictive_score")]
    x <- merge(a, b, by = c("replication", "school"), suffixes = c("_a", "_b"))
    if (nrow(x) != nrow(a) || nrow(x) != nrow(b)) stop("Unpaired folds.")
    differences <- x$log_predictive_score_a - x$log_predictive_score_b
    do.call(rbind, lapply(split(differences, x$replication), function(d)
      data.frame(model = pair[1L], comparator = pair[2L],
        mean_log_score_difference = mean(d)))) |>
      transform(replication = as.integer(names(split(differences, x$replication))))
  }))
  paired <- do.call(rbind, lapply(split(paired_reps,
    interaction(paired_reps$model, paired_reps$comparator, drop = TRUE)), function(x) {
      d <- x$mean_log_score_difference
      data.frame(model = x$model[1L], comparator = x$comparator[1L],
        mean_log_score_difference = mean(d), geometric_density_ratio = exp(mean(d)),
        monte_carlo_se = sd(d) / sqrt(length(d)), replications = length(d))
    }))
  list(replication = reps, summary = summary,
       paired_replication = paired_reps, paired_summary = paired)
}
