# Compute the manuscript's log score from previously saved latent metrics.
# Requires only base R; no model estimation or new predictive simulation.
source("STAR_latent_predictive_helpers.R")
output_dir <- "star_output/loso_latent_metrics"
comparison <- readRDS(file.path(output_dir, "LOSO_vs_full_latent_ATE_draws.rds"))
metrics <- read.csv(file.path(output_dir, "LOSO_latent_school_metrics.csv"),
                    colClasses = c(original_schoolid = "character"))
stopifnot(length(comparison) == 79L, nrow(metrics) == 79L,
          identical(as.integer(metrics$fold), seq_len(79L)))
results <- vector("list", length(comparison))
for (j in seq_along(comparison)) {
  saved <- comparison[[j]]
  fold <- readRDS(sprintf("star_output/loso/loso_school_%02d.rds", j))
  stopifnot(saved$fold == j, fold$metadata$omitted_school == j,
            as.character(saved$original_schoolid) == metrics$original_schoolid[j],
            as.character(fold$metadata$original_schoolid) == metrics$original_schoolid[j],
            isTRUE(all.equal(as.numeric(fold$metadata$outcome_scale),
                             metrics$fold_outcome_scale[j])))
  prediction <- latent_loso_ate_moments(fold)
  conversion <- metrics$fold_outcome_scale[j] / metrics$full_outcome_scale[j]
  stopifnot(is.finite(conversion), conversion > 0)
  log_density <- latent_log_predictive_density(
    saved$full_posterior_latent_ATE,
    prediction$conditional_locations * conversion,
    prediction$conditional_sds * conversion
  )
  results[[j]] <- data.frame(
    fold = j, original_schoolid = metrics$original_schoolid[j],
    log_predictive_score = mean(log_density),
    log_score_full_draws = length(log_density),
    log_score_loso_draws = length(prediction$conditional_locations)
  )
  message("Scored school ", j, " of ", length(comparison))
}
scores <- do.call(rbind, results)
summary <- data.frame(
  metric = "Mean posterior-expected latent ATE log predictive score",
  value = mean(scores$log_predictive_score), schools = nrow(scores),
  weighting = "Equal weight per school",
  scale = "Full-fit control-group standard deviation",
  logarithm = "Natural log",
  interpretation = "Higher is better for a fixed full-data reference posterior"
)
write.csv(scores, file.path(output_dir, "LOSO_latent_log_score_by_school.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_dir, "LOSO_latent_log_score_summary.csv"),
          row.names = FALSE)
print(summary, row.names = FALSE)
