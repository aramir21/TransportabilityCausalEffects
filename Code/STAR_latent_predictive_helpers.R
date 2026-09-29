# Normal component parameters for each LOSO hyperposterior draw.
# Gamma is stored as Kh by Kz, so the mean is zbar' Gamma' h.
latent_loso_ate_moments <- function(fold) {
  Gamma_draws <- fold$posterior$Gamma
  Sigma_draws <- fold$posterior$Sigma_delta
  h <- as.numeric(fold$metadata$h_heldout)
  zbar <- as.numeric(fold$metadata$zbar_heldout)

  dimensions_Gamma <- dim(Gamma_draws)
  dimensions_Sigma <- dim(Sigma_draws)
  if (length(dimensions_Gamma) != 3L || length(dimensions_Sigma) != 3L) {
    stop("Gamma or Sigma_delta has invalid dimensions.")
  }

  S <- dimensions_Gamma[1L]
  Kh <- dimensions_Gamma[2L]
  Kz <- dimensions_Gamma[3L]
  if (length(h) != Kh || length(zbar) != Kz ||
      !identical(dimensions_Sigma, c(S, Kz, Kz))) {
    stop("Posterior and held-out covariate dimensions do not agree.")
  }

  # Conditional predictive mean and variance of the latent school ATE under
  # each hyperposterior draw.
  tau_location <- numeric(S)
  tau_variance <- numeric(S)
  for (s in seq_len(S)) {
    Gamma_s <- matrix(Gamma_draws[s, , ], nrow = Kh, ncol = Kz)
    Sigma_s <- matrix(Sigma_draws[s, , ], nrow = Kz, ncol = Kz)
    mean_delta_s <- drop(t(Gamma_s) %*% h)
    tau_location[s] <- drop(crossprod(zbar, mean_delta_s))
    tau_variance[s] <- drop(crossprod(zbar, Sigma_s %*% zbar))
  }

  tolerance <- 1e-10
  if (any(tau_variance < -tolerance) || any(!is.finite(tau_variance))) {
    stop("A latent predictive ATE variance is negative or nonfinite.")
  }
  tau_variance <- pmax(tau_variance, 0)

  list(
    conditional_locations = tau_location,
    conditional_sds = sqrt(tau_variance)
  )
}

# Integrate the normal density over LOSO hyperposterior draws BEFORE taking
# logs. The caller averages these logs over full-data ATE draws.
# Log-sum-exp avoids underflow;
# processing one evaluation point at a time avoids an S by R allocation.
latent_log_predictive_density <- function(full_draws, locations, sds) {
  if (!length(full_draws) || !length(locations) ||
      length(locations) != length(sds) ||
      any(!is.finite(c(full_draws, locations, sds))) || any(sds <= 0)) {
    stop("Log scoring requires finite draws and strictly positive component SDs.")
  }
  vapply(full_draws, function(tau) {
    log_components <- stats::dnorm(tau, mean = locations, sd = sds, log = TRUE)
    maximum <- max(log_components)
    if (!is.finite(maximum)) return(maximum)
    maximum + log(mean(exp(log_components - maximum)))
  }, numeric(1))
}
