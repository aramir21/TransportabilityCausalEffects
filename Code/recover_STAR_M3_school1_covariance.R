# Recover parameters omitted from the original pilot checkpoint, keeping
# that original checkpoint intact. Same data, model, sampler settings and seed.
source("STAR_model_comparison_helpers.R")
old <- readRDS("star_output/model_comparison/M3_school_01.rds")
cfg <- old$config
model <- rstan::stan_model(model_code = paste(old$signature$model_code, collapse = "\n"))
fit <- rstan::sampling(model, data = old$signature$stan_data,
  chains = cfg$chains, iter = cfg$iter, warmup = cfg$warmup, cores = cfg$cores,
  seed = cfg$seed + 3000L + 1L, refresh = 0L,
  control = list(adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_treedepth))
posterior <- rstan::extract(fit, pars = c("Sigma", "sd_delta", "L_Omega",
                                       "predictive_location", "predictive_sd"))
L_Sigma <- array(0, dim = dim(posterior$Sigma))
for (r in seq_len(dim(L_Sigma)[1])) {
  L_Sigma[r, , ] <- diag(posterior$sd_delta[r, ]) %*% posterior$L_Omega[r, , ]
  stopifnot(max(abs(tcrossprod(L_Sigma[r, , ]) - posterior$Sigma[r, , ])) < 1e-10)
}
posterior$L_Sigma <- L_Sigma
agreement <- c(
  location_max_sorted_difference = max(abs(sort(old$locations) - sort(posterior$predictive_location))),
  sd_max_sorted_difference = max(abs(sort(old$sds) - sort(posterior$predictive_sd))))
result <- list(posterior = posterior, original_signature = old$signature,
  outcome_scale = old$outcome_scale, agreement = agreement,
  diagnostics = comparison_diagnostics(fit, cfg$max_treedepth, "M3"))
saveRDS(result, "star_output/model_comparison/M3_school_01_covariance.rds")
print(agreement)
print(result$diagnostics)
