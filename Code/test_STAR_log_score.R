source("STAR_latent_predictive_helpers.R")
x <- c(-2, 0, 3)
stopifnot(isTRUE(all.equal(
  latent_log_predictive_density(x, c(0, 0), c(1, 1)),
  dnorm(x, log = TRUE)
)))
stopifnot(isTRUE(all.equal(
  latent_log_predictive_density(x, c(-1, 2), c(0.5, 1.5)),
  log((dnorm(x, -1, 0.5) + dnorm(x, 2, 1.5)) / 2)
)))
stopifnot(isTRUE(all.equal(
  latent_log_predictive_density(2 * x, c(0, 0), c(2, 2)),
  dnorm(x, log = TRUE) - log(2)
)))
stopifnot(is.finite(latent_log_predictive_density(1000, c(0, 1), c(1, 2))))
stopifnot(inherits(try(
  latent_log_predictive_density(x, 0, 0), silent = TRUE
), "try-error"))
cat("Log-score numerical checks passed\n")
