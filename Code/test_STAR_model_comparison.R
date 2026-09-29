source("STAR_model_comparison_helpers.R")

# Closed-form cross entropy of two normals, including the scaling Jacobian.
score <- score_school_reference(c(.3, .3), c(.7, .7), -.1, .2)
expected <- -.5 * log(2 * pi * .7^2) - (.2^2 + (-.1 - .3)^2) / (2 * .7^2)
stopifnot(abs(score$log_predictive_score - expected) < 1e-7,
          abs(score$lower_95 - qnorm(.025, .3, .7)) < 1e-8)
rescaled <- score_school_reference(c(.6, .6), c(1.4, 1.4), -.2, .4)
stopifnot(abs(rescaled$log_predictive_score - score$log_predictive_score + log(2)) < 1e-7,
          abs(rescaled$posterior_coverage_95 - score$posterior_coverage_95) < 1e-8)

inputs <- read_comparison_inputs("star_output/loso")
for (j in c(1L, 40L, 79L)) {
  m1 <- aggregate_fold_data(inputs, j, "M1")
  m2 <- aggregate_fold_data(inputs, j, "M2")
  m3 <- aggregate_fold_data(inputs, j, "M3")
  stopifnot(m1$K == 4L, m2$K == 6L, m3$N == 78L,
            identical(m2$X, m3$X))
  h <- as.numeric(inputs$metadata[[j]]$h_heldout)
  p <- inputs$schools$p[j]
  Gamma <- matrix(seq_len(6) / 10, 3, 2)
  stopifnot(abs(sum(m3$x_new * as.vector(Gamma)) -
                drop(t(c(1, p)) %*% t(Gamma) %*% h)) < 1e-12)
  Sigma <- matrix(c(.2, -.03, -.03, .1), 2)
  stopifnot(abs(drop(t(c(1, p)) %*% Sigma %*% c(1, p)) -
                (Sigma[1, 1] + 2 * p * Sigma[1, 2] + p^2 * Sigma[2, 2])) < 1e-12)
  changed <- inputs
  changed$schools$estimate_raw[j] <- 1e9
  changed$schools$se_raw[j] <- 1e8
  stopifnot(identical(m3, aggregate_fold_data(changed, j, "M3")))
}
for (path in c("ApplicationSTARV3_model_comparison_fit.R",
               "ApplicationSTARV3_model_comparison_metrics.R")) invisible(parse(path))
cat("Model-comparison analytic, scaling, design, and held-out exclusion checks passed.\n")
