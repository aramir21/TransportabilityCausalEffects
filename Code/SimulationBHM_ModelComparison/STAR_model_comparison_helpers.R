source("STAR_latent_predictive_helpers.R")

# Recover school summaries in raw outcome units from existing M0 checkpoints.
# Invert each school's saved H transformation, then apply the TRAINING fold's
# transformation to all schools when fitting aggregate models.
read_comparison_inputs <- function(loso_dir, expected_schools = 79L) {
  files <- file.path(loso_dir, sprintf("loso_school_%02d.rds", seq_len(expected_schools)))
  if (!all(file.exists(files))) stop("Missing M0 LOSO checkpoints.")
  metadata <- lapply(files, function(path) readRDS(path)$metadata)
  h_names <- names(metadata[[1]]$H_scaling$center)
  if (!identical(h_names, c("prop_teacher_eth_cauc", "prop_degree_master+")))
    stop("This comparison preserves the existing M0 teacher covariates.")
  rows <- lapply(seq_along(metadata), function(j) {
    m <- metadata[[j]]
    stopifnot(m$omitted_school == j,
              is.finite(m$outcome_scale), m$outcome_scale > 0,
              all(is.finite(m$H_scaling$scale)), all(m$H_scaling$scale > 0),
              identical(names(m$H_scaling$center), h_names),
              identical(names(m$zbar_heldout), c("intercept", "free_lunch")))
    raw_h <- m$h_heldout[h_names] * m$H_scaling$scale + m$H_scaling$center
    values <- c(raw_h, p = unname(m$zbar_heldout[2]),
                estimate_raw = unname(m$direct_ATE["estimate"] * m$outcome_scale),
                se_raw = unname(m$direct_ATE["se"] * m$outcome_scale))
    if (any(!is.finite(values)) || values["se_raw"] <= 0 ||
        values["p"] < 0 || values["p"] > 1) stop("Invalid school summaries: ", j)
    data.frame(fold = j, original_schoolid = as.character(m$original_schoolid),
               as.list(values), check.names = FALSE)
  })
  schools <- do.call(rbind, rows)
  if (anyDuplicated(schools$original_schoolid)) stop("Duplicate school IDs.")
  list(schools = schools, metadata = metadata, h_names = h_names, files = files)
}

aggregate_fold_data <- function(inputs, j, model) {
  stopifnot(model %in% c("M1", "M2", "M3"))
  d <- inputs$schools
  m <- inputs$metadata[[j]]
  h <- sweep(as.matrix(d[inputs$h_names]), 2, m$H_scaling$center, "-")
  h <- sweep(h, 2, m$H_scaling$scale, "/")
  h <- cbind(intercept = 1, h)
  # M2/M3 columns have the same order and N(0, .5^2) priors as vec(Gamma).
  x <- if (model == "M1") cbind(h, p = d$p) else cbind(h, h * d$p)
  colnames(x) <- if (model == "M1") c(colnames(h), "p") else
    c(colnames(h), paste0("p:", colnames(h)))
  train <- seq_len(nrow(d)) != j
  list(N = sum(train), K = ncol(x), composition_variance = as.integer(model == "M3"),
       X = x[train, , drop = FALSE], estimate = d$estimate_raw[train] / m$outcome_scale,
       sampling_se = d$se_raw[train] / m$outcome_scale, p = d$p[train],
       x_new = as.numeric(x[j, ]), p_new = d$p[j])
}

comparison_diagnostics <- function(fit, max_treedepth, model) {
  parameters <- c("coefficients", "predictive_location", "predictive_sd",
                  if (model == "M3") c("sd_delta", "Sigma") else "common_sd")
  s <- rstan::summary(fit, pars = parameters)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  list(max_rhat = max(s[, "Rhat"]), min_n_eff = min(s[, "n_eff"]),
       divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
       treedepth_hits = sum(vapply(sampler, function(x)
         sum(x[, "treedepth__"] >= max_treedepth), numeric(1))),
       min_ebfmi = min(vapply(sampler, function(x)
         mean(diff(x[, "energy__"])^2) / var(x[, "energy__"]), numeric(1))))
}

# Deterministic normal-mixture quantiles: no predictive replication required.
mixture_quantile <- function(probability, locations, sds) {
  endpoints <- locations + qnorm(probability) * sds
  bounds <- range(endpoints)
  if (diff(bounds) < 1e-12) return(mean(bounds))
  uniroot(function(x) mean(pnorm(x, locations, sds)) - probability,
          interval = bounds, tol = 1e-9)$root
}

# School-only reference: N(estimate, se^2), from normal sampling likelihood
# and an improper flat prior. The latent predictive SD excludes sampling SE.
score_school_reference <- function(locations, sds, estimate, se,
                                   levels = c(.5, .8, .9, .95)) {
  stopifnot(length(locations) == length(sds), length(locations) > 0,
            all(is.finite(c(locations, sds, estimate, se))), all(sds > 0), se > 0)
  log_score <- integrate(function(z) {
    latent_log_predictive_density(estimate + se * z, locations, sds) * dnorm(z)
  }, lower = -Inf, upper = Inf, rel.tol = 1e-6, abs.tol = 1e-7,
     subdivisions = 300L, stop.on.error = TRUE)
  error <- mean(locations) - estimate
  result <- data.frame(predictive_mean = mean(locations), reference_mean = estimate,
                       reference_sd = se, squared_mean_error = error^2,
                       posterior_expected_squared_error = error^2 + se^2,
                       log_predictive_score = log_score$value,
                       integration_abs_error = log_score$abs.error)
  for (level in levels) {
    label <- sprintf("%02d", round(100 * level))
    lower <- mixture_quantile((1 - level) / 2, locations, sds)
    upper <- mixture_quantile((1 + level) / 2, locations, sds)
    result[[paste0("lower_", label)]] <- lower
    result[[paste0("upper_", label)]] <- upper
    result[[paste0("length_", label)]] <- upper - lower
    result[[paste0("posterior_coverage_", label)]] <-
      pnorm(upper, estimate, se) - pnorm(lower, estimate, se)
    result[[paste0("reference_mean_covered_", label)]] <-
      as.integer(estimate >= lower && estimate <= upper)
  }
  result
}
