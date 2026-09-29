# Convergence and HMC diagnostics for the full Project STAR fit

diagnostics_config <- list(
  fit_path = file.path(
    "star_output",
    "STAR_kindergarten_reading_fit.RData"
  ),
  output_dir = file.path("star_output", "full_fit_diagnostics"),
  max_treedepth = 10L,
  rhat_warning = 1.01,
  ebfmi_warning = 0.30,
  # Exclude the large school-level delta and delta_raw blocks from the concise
  # parameter table. Add them here if their individual diagnostics are needed.
  parameters = c(
    "beta", "gamma", "Gamma", "sd_delta",
    "Omega", "Sigma_delta", "sigma"
  )
)

if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("Install the rstan package before running this script.")
}
if (!file.exists(diagnostics_config$fit_path)) {
  stop("Saved fit not found: ", diagnostics_config$fit_path)
}
dir.create(diagnostics_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load into an isolated environment so existing workspace objects are unchanged.
fit_environment <- new.env(parent = emptyenv())
loaded_objects <- load(diagnostics_config$fit_path, envir = fit_environment)
if (!"fit" %in% loaded_objects) {
  stop("The RData file does not contain an object named 'fit'.")
}
fit <- fit_environment$fit
if (!methods::is(fit, "stanfit")) {
  stop("The object named 'fit' is not a stanfit object.")
}

# -----------------------------------------------------------------------------
# Parameter-level convergence diagnostics
# -----------------------------------------------------------------------------

available_parameters <- dimnames(as.array(fit))[[3L]]
requested_regex <- paste0(
  "^(",
  paste(diagnostics_config$parameters, collapse = "|"),
  ")($|\\[)"
)
selected_parameter_names <- grep(
  requested_regex,
  available_parameters,
  value = TRUE
)
if (length(selected_parameter_names) == 0L) {
  stop("None of the requested parameters are present in the saved fit.")
}

# Use rank-normalized split-R-hat and bulk/tail ESS when posterior is available.
# Otherwise fall back to the R-hat and n_eff reported by rstan.
if (requireNamespace("posterior", quietly = TRUE)) {
  selected_array <- as.array(fit, pars = selected_parameter_names)
  selected_draws <- posterior::as_draws_array(selected_array)
  parameter_diagnostics <- as.data.frame(
    posterior::summarise_draws(
      selected_draws,
      "mean", "sd", "rhat", "ess_bulk", "ess_tail"
    )
  )
  names(parameter_diagnostics)[names(parameter_diagnostics) == "variable"] <-
    "parameter"
  ess_columns <- c("ess_bulk", "ess_tail")
  ess_method <- "Rank-normalized bulk and tail ESS from posterior"
} else {
  rstan_summary <- rstan::summary(
    fit,
    pars = diagnostics_config$parameters
  )$summary
  parameter_diagnostics <- data.frame(
    parameter = rownames(rstan_summary),
    mean = rstan_summary[, "mean"],
    sd = rstan_summary[, "sd"],
    rhat = rstan_summary[, "Rhat"],
    n_eff = rstan_summary[, "n_eff"],
    row.names = NULL,
    check.names = FALSE
  )
  ess_columns <- "n_eff"
  ess_method <- "Traditional effective sample size from rstan"
}

parameter_diagnostics$rhat_warning <-
  is.finite(parameter_diagnostics$rhat) &
  parameter_diagnostics$rhat > diagnostics_config$rhat_warning

utils::write.csv(
  parameter_diagnostics,
  file.path(
    diagnostics_config$output_dir,
    "full_fit_parameter_diagnostics.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Chain-level HMC diagnostics
# -----------------------------------------------------------------------------

sampler_parameters <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
if (length(sampler_parameters) == 0L) {
  stop("No post-warmup sampler diagnostics were found.")
}

chain_diagnostics <- do.call(rbind, lapply(
  seq_along(sampler_parameters),
  function(chain) {
    sampler <- sampler_parameters[[chain]]
    energy <- sampler[, "energy__"]
    ebfmi <- mean(diff(energy)^2) / stats::var(energy)

    data.frame(
      chain = chain,
      post_warmup_draws = nrow(sampler),
      divergences = sum(sampler[, "divergent__"]),
      treedepth_hits = sum(
        sampler[, "treedepth__"] >= diagnostics_config$max_treedepth
      ),
      max_observed_treedepth = max(sampler[, "treedepth__"]),
      mean_accept_stat = mean(sampler[, "accept_stat__"]),
      median_leapfrog_steps = stats::median(sampler[, "n_leapfrog__"]),
      max_leapfrog_steps = max(sampler[, "n_leapfrog__"]),
      mean_stepsize = mean(sampler[, "stepsize__"]),
      ebfmi = ebfmi,
      ebfmi_warning = !is.finite(ebfmi) ||
        ebfmi < diagnostics_config$ebfmi_warning
    )
  }
))

utils::write.csv(
  chain_diagnostics,
  file.path(
    diagnostics_config$output_dir,
    "full_fit_chain_diagnostics.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Concise overall summary
# -----------------------------------------------------------------------------

finite_rhat <- is.finite(parameter_diagnostics$rhat)
if (any(finite_rhat)) {
  worst_rhat_row <- which.max(ifelse(
    finite_rhat,
    parameter_diagnostics$rhat,
    -Inf
  ))
  maximum_rhat <- parameter_diagnostics$rhat[worst_rhat_row]
  maximum_rhat_parameter <- parameter_diagnostics$parameter[worst_rhat_row]
} else {
  maximum_rhat <- NA_real_
  maximum_rhat_parameter <- NA_character_
}

minimum_ess <- min(
  unlist(parameter_diagnostics[ess_columns], use.names = FALSE),
  na.rm = TRUE
)
if (!is.finite(minimum_ess)) minimum_ess <- NA_real_

overall_diagnostics <- data.frame(
  diagnostic = c(
    "chains",
    "post_warmup_draws_per_chain",
    "maximum_rhat",
    "parameters_above_rhat_threshold",
    "minimum_effective_sample_size",
    "total_divergences",
    "total_treedepth_hits",
    "minimum_ebfmi",
    "chains_below_ebfmi_threshold"
  ),
  value = c(
    nrow(chain_diagnostics),
    min(chain_diagnostics$post_warmup_draws),
    maximum_rhat,
    sum(parameter_diagnostics$rhat_warning),
    minimum_ess,
    sum(chain_diagnostics$divergences),
    sum(chain_diagnostics$treedepth_hits),
    min(chain_diagnostics$ebfmi),
    sum(chain_diagnostics$ebfmi_warning)
  ),
  detail = c(
    "",
    "",
    maximum_rhat_parameter,
    paste0("threshold = ", diagnostics_config$rhat_warning),
    ess_method,
    "",
    paste0("threshold = ", diagnostics_config$max_treedepth),
    "",
    paste0("threshold = ", diagnostics_config$ebfmi_warning)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  overall_diagnostics,
  file.path(
    diagnostics_config$output_dir,
    "full_fit_overall_diagnostics.csv"
  ),
  row.names = FALSE
)

cat("\nOverall convergence and HMC diagnostics\n")
print(overall_diagnostics, row.names = FALSE)
cat("\nChain-level HMC diagnostics\n")
print(chain_diagnostics, row.names = FALSE)

# Also run rstan's standard diagnostic messages in the Console.
rstan::check_hmc_diagnostics(fit)

message(
  "Diagnostic files written to: ",
  normalizePath(diagnostics_config$output_dir)
)
