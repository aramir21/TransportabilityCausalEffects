# Compile and fit one fold of each model using the production sampler settings.
stopifnot(requireNamespace("rstan", quietly = TRUE))
translation <- rstan::stanc(file = "STAR_model_comparison.stan")
stopifnot(isTRUE(translation$status))
options(star.comparison = list(schools = 1L))
source("ApplicationSTARV3_model_comparison_fit.R")
