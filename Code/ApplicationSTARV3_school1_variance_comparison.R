# All quantities on the common full-sample control-SD scale; all fits omit school 1.
source("STAR_latent_predictive_helpers.R")
output_dir <- "star_output/model_comparison/metrics_schools_1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
full_sd <- read.csv("star_output/loso_latent_metrics/LOSO_latent_school_metrics.csv")$full_outcome_scale[1]
m0 <- readRDS("star_output/loso/loso_school_01.rds")
a0 <- m0$metadata$outcome_scale / full_sd
p <- as.numeric(m0$metadata$zbar_heldout[2])
summaries <- list()
add_summary <- function(model, parameter, x) {
  summaries[[length(summaries) + 1L]] <<- data.frame(model = model, parameter = parameter,
    mean = mean(x), median = median(x), sd = sd(x),
    lower_95 = unname(quantile(x, .025)), upper_95 = unname(quantile(x, .975)))
}
summarise_matrix <- function(model, name, x) {
  for (i in 1:2) for (j in 1:i)
    add_summary(model, sprintf("%s[%d,%d]", name, i, j), x[, i, j])
}
Sigma0 <- m0$posterior$Sigma_delta * a0^2
summarise_matrix("M0", "Sigma", Sigma0)
moments0 <- latent_loso_ate_moments(m0)
predictions <- list(M0 = list(locations = moments0$conditional_locations * a0,
                              sds = moments0$conditional_sds * a0))
for (model in c("M1", "M2", "M3")) {
  fit <- readRDS(sprintf("star_output/model_comparison/%s_school_01.rds", model))
  a <- fit$outcome_scale / full_sd
  predictions[[model]] <- list(locations = fit$locations * a, sds = fit$sds * a)
  if (model != "M3") {
    add_summary(model, "common_sd", fit$sds * a)
    add_summary(model, "common_variance", (fit$sds * a)^2)
  }
}
recovered_path <- "star_output/model_comparison/M3_school_01_covariance.rds"
if (file.exists(recovered_path)) {
  recovered <- readRDS(recovered_path)
  a <- recovered$outcome_scale / full_sd
  Sigma3 <- recovered$posterior$Sigma * a^2
  L3 <- recovered$posterior$L_Sigma * a
  summarise_matrix("M3", "Sigma", Sigma3)
  summarise_matrix("M3", "L_Sigma", L3)
  print(recovered$agreement)
}
reference_mean <- unname(m0$metadata$direct_ATE["estimate"] * a0)
reference_sd <- unname(m0$metadata$direct_ATE["se"] * a0)
decomposition <- do.call(rbind, lapply(names(predictions), function(model) {
  pred <- predictions[[model]]
  residual_variance <- mean(pred$sds^2)
  mean_uncertainty <- mean((pred$locations - mean(pred$locations))^2)
  data.frame(model = model, predictive_mean = mean(pred$locations),
    reference_mean = reference_mean, reference_sd = reference_sd,
    mean_conditional_sd = mean(pred$sds),
    expected_residual_variance = residual_variance,
    variance_of_conditional_mean = mean_uncertainty,
    total_predictive_variance = residual_variance + mean_uncertainty,
    predictive_sd = sqrt(residual_variance + mean_uncertainty))
}))
parameters <- do.call(rbind, summaries)
write.csv(parameters, file.path(output_dir, "posterior_variance_parameters.csv"), row.names = FALSE)
write.csv(decomposition, file.path(output_dir, "predictive_variance_decomposition.csv"), row.names = FALSE)
cat("Free-lunch proportion:", p, "\nCommon reporting SD:", full_sd, "\n")
print(parameters, row.names = FALSE)
print(decomposition, row.names = FALSE)
