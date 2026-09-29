# All four predictions evaluated against ONE common, school-only reference.
source("STAR_model_comparison_helpers.R")
comparison_metrics_config <- modifyList(list(
  loso_dir = "star_output/loso", output_dir = "star_output/model_comparison",
  full_metrics_path = "star_output/loso_latent_metrics/LOSO_latent_school_metrics.csv",
  schools = seq_len(79L), fit_dir = NULL,
  models = c("M0", "M1", "M2", "M3"), levels = c(.5, .8, .9, .95)
), getOption("star.comparison.metrics", list()))
cfg <- comparison_metrics_config
if (is.null(cfg$fit_dir)) cfg$fit_dir <- cfg$output_dir
stopifnot(all(cfg$models %in% c("M0", "M1", "M2", "M3")),
          length(cfg$schools) > 0L, !anyDuplicated(cfg$schools),
          all(cfg$schools %in% seq_len(79L)),
          "M0" %in% cfg$models, !anyDuplicated(cfg$models),
          all(cfg$levels > 0 & cfg$levels < 1))
inputs <- read_comparison_inputs(cfg$loso_dir)
# The common scale is only a reporting unit; no full-data posterior is used.
previous <- read.csv(cfg$full_metrics_path, colClasses = c(original_schoolid = "character"))
stopifnot(identical(as.integer(previous$fold), inputs$schools$fold),
          identical(previous$original_schoolid, inputs$schools$original_schoolid),
          isTRUE(all.equal(previous$fold_outcome_scale,
            vapply(inputs$metadata, function(m) m$outcome_scale, numeric(1)))))
reporting_sd <- unique(previous$full_outcome_scale)
stopifnot(length(reporting_sd) == 1L, is.finite(reporting_sd), reporting_sd > 0)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
# Fail before scoring if any requested aggregate fit is missing.
for (model in setdiff(cfg$models, "M0")) {
  paths <- file.path(cfg$fit_dir, sprintf("%s_school_%02d.rds", model, cfg$schools))
  if (!all(file.exists(paths))) stop("Missing ", model, " checkpoints; run the fitting script first.")
}
rows <- list()
diagnostics <- list()
for (j in cfg$schools) {
  reference_mean <- inputs$schools$estimate_raw[j] / reporting_sd
  reference_se <- inputs$schools$se_raw[j] / reporting_sd
  for (model in cfg$models) {
    message("Scoring ", model, " school ", j, "...")
    if (model == "M0") {
      fold <- readRDS(inputs$files[j])
      prediction <- latent_loso_ate_moments(fold)
      locations <- prediction$conditional_locations
      sds <- prediction$conditional_sds
      scale <- fold$metadata$outcome_scale
      diagnostic <- fold$diagnostics
    } else {
      fold <- readRDS(file.path(cfg$fit_dir, sprintf("%s_school_%02d.rds", model, j)))
      stopifnot(fold$fold == j, fold$model == model,
                identical(fold$signature$stan_data, aggregate_fold_data(inputs, j, model)),
                identical(fold$signature$model_code, readLines("STAR_model_comparison.stan")),
                fold$signature$original_schoolid == inputs$schools$original_schoolid[j],
                isTRUE(all.equal(fold$outcome_scale, inputs$metadata[[j]]$outcome_scale)))
      locations <- fold$locations
      sds <- fold$sds
      scale <- fold$outcome_scale
      diagnostic <- fold$diagnostics
    }
    score <- score_school_reference(locations * scale / reporting_sd,
                                    sds * scale / reporting_sd,
                                    reference_mean, reference_se, cfg$levels)
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(fold = j, original_schoolid = inputs$schools$original_schoolid[j], model = model), score)
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      fold = j, model = model, max_rhat = diagnostic$max_rhat,
      min_n_eff = diagnostic$min_n_eff, divergences = diagnostic$divergences,
      treedepth_hits = diagnostic$treedepth_hits, min_ebfmi = diagnostic$min_ebfmi)
    message("Scored ", model, " school ", j)
  }
}
school_metrics <- do.call(rbind, rows)
baseline <- school_metrics[school_metrics$model == "M0", ]
school_metrics$log_score_difference_vs_M0 <- school_metrics$log_predictive_score -
  baseline$log_predictive_score[match(school_metrics$fold, baseline$fold)]
summary <- do.call(rbind, lapply(cfg$models, function(model) {
  d <- school_metrics[school_metrics$model == model, ]
  data.frame(model = model, schools = nrow(d), mean_log_score = mean(d$log_predictive_score),
    mean_log_score_difference_vs_M0 = mean(d$log_score_difference_vs_M0),
    schools_with_higher_log_score_than_M0 = sum(d$log_score_difference_vs_M0 > 0),
    RMSE_reference_means = sqrt(mean(d$squared_mean_error)),
    posterior_expected_RMSE = sqrt(mean(d$posterior_expected_squared_error)))
}))
coverage <- do.call(rbind, lapply(cfg$models, function(model) {
  d <- school_metrics[school_metrics$model == model, ]
  do.call(rbind, lapply(cfg$levels, function(level) {
    label <- sprintf("%02d", round(100 * level))
    data.frame(model = model, level = level,
      posterior_coverage = mean(d[[paste0("posterior_coverage_", label)]]),
      reference_mean_coverage = mean(d[[paste0("reference_mean_covered_", label)]]),
      mean_interval_length = mean(d[[paste0("length_", label)]]),
      median_interval_length = median(d[[paste0("length_", label)]]))
  }))
}))
# Include each interval level in the one-row-per-model summary as well.
for (level in cfg$levels) {
  label <- sprintf("%02d", round(100 * level))
  level_coverage <- coverage[coverage$level == level, ]
  level_coverage <- level_coverage[match(summary$model, level_coverage$model), ]
  for (measure in c("posterior_coverage", "reference_mean_coverage",
                    "mean_interval_length", "median_interval_length")) {
    summary[[paste0(measure, "_", label)]] <- level_coverage[[measure]]
  }
}
if (length(cfg$models) > 1L) {
  pairs <- combn(cfg$models, 2, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    a <- school_metrics[school_metrics$model == pair[1], ]
    b <- school_metrics[school_metrics$model == pair[2], ]
    difference <- b$log_predictive_score[match(a$fold, b$fold)] - a$log_predictive_score
    data.frame(model = pair[2], comparator = pair[1], schools = length(difference),
               mean_log_score_difference = mean(difference),
               median_log_score_difference = median(difference),
               schools_with_higher_log_score = sum(difference > 0))
  }))
  write.csv(pairwise, file.path(cfg$output_dir, "comparison_pairwise_log_scores.csv"), row.names = FALSE)
}
write.csv(school_metrics, file.path(cfg$output_dir, "comparison_school_metrics.csv"), row.names = FALSE)
write.csv(summary, file.path(cfg$output_dir, "comparison_summary.csv"), row.names = FALSE)
write.csv(coverage, file.path(cfg$output_dir, "comparison_coverage.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(cfg$output_dir, "comparison_diagnostics.csv"), row.names = FALSE)
saveRDS(list(config = cfg, reporting_sd = reporting_sd,
             reference = "School-only normal likelihood with flat prior",
             h_names = inputs$h_names), file.path(cfg$output_dir, "comparison_metadata.rds"))
print(summary, row.names = FALSE)
