# Leave-one-school-out estimation for ApplicationSTARV3_clean.R
#
# This script estimates J models. Fold j excludes all outcomes from school j.
# Training-dependent transformations are recomputed using the remaining J - 1
# schools. Each completed fold is saved immediately and is skipped on restart.

# -----------------------------------------------------------------------------
# LOSO configuration
# -----------------------------------------------------------------------------

loso_config <- list(
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  refresh = 100L,
  seed = 20260821L,
  cores = 4L,
  schools = NULL,             # NULL fits all; use 1L for an initial test
  overwrite = FALSE,
  save_full_stanfit = FALSE,
  output_dir = "star_output/loso"
)

if (loso_config$warmup >= loso_config$iter) {
  stop("loso_config$warmup must be smaller than loso_config$iter.")
}
dir.create(loso_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the common data preparation and Stan program without fitting the
# full-sample model or creating its posterior summaries.
old_prepare_option <- getOption("star.prepare_only")
options(star.prepare_only = TRUE)
source("ApplicationSTARV3_clean.R", local = FALSE)
options(star.prepare_only = old_prepare_option)

# -----------------------------------------------------------------------------
# Fold construction
# -----------------------------------------------------------------------------

scaling_parameters <- function(data, variables) {
  centers <- vapply(data[variables], mean, numeric(1))
  scales <- vapply(data[variables], stats::sd, numeric(1))
  if (any(!is.finite(centers)) || any(!is.finite(scales)) || any(scales <= 0)) {
    bad <- variables[!is.finite(centers) | !is.finite(scales) | scales <= 0]
    stop("Invalid training-fold scale for: ", paste(bad, collapse = ", "))
  }
  list(center = centers, scale = scales)
}

apply_scaling <- function(data, variables, parameters) {
  values <- as.matrix(data[variables])
  values <- sweep(values, 2L, parameters$center, FUN = "-")
  values <- sweep(values, 2L, parameters$scale, FUN = "/")
  colnames(values) <- variables
  values
}

make_loso_fold <- function(omitted_school, data, w_names, h_names) {
  train <- data[data$stan_school != omitted_school, , drop = FALSE]
  heldout <- data[data$stan_school == omitted_school, , drop = FALSE]
  if (nrow(train) == 0L || nrow(heldout) == 0L) {
    stop("Empty training or held-out sample for school ", omitted_school, ".")
  }

  # A category can become constant after its only school is omitted. Remove
  # such W columns within that fold because their coefficients are unidentified.
  w_names_fold <- w_names[vapply(
    train[w_names],
    function(x) is.finite(stats::sd(x)) && stats::sd(x) > 0,
    logical(1)
  )]
  if (length(w_names_fold) == 0L) {
    stop("No nonconstant W regressors in fold ", omitted_school)
  }

  age_scaling <- scaling_parameters(train, "age")
  w_scaling <- scaling_parameters(train, w_names_fold)
  h_scaling <- scaling_parameters(train, h_names)

  X_train <- cbind(
    intercept = 1,
    female = train$female,
    eth_cauc = train$eth_cauc,
    eth_afam = train$eth_afam,
    age = apply_scaling(train, "age", age_scaling)[, 1L],
    Q1 = train$Q1,
    Q2 = train$Q2,
    Q3 = train$Q3,
    free_lunch = train$free_lunch
  )
  Z_train <- X_train[, c("intercept", "free_lunch"), drop = FALSE]
  W_train <- apply_scaling(train, w_names_fold, w_scaling)
  H_train <- cbind(
    intercept = 1,
    apply_scaling(train, h_names, h_scaling)
  )
  colnames(H_train) <- c("intercept", h_names)

  control_train <- train$Tr == 0L
  if (sum(control_train) < 2L) stop("Too few training controls in fold ", omitted_school)
  outcome_mean <- mean(train$y_raw[control_train])
  outcome_sd <- stats::sd(train$y_raw[control_train])
  if (!is.finite(outcome_sd) || outcome_sd <= 0) {
    stop("Invalid training outcome SD in fold ", omitted_school)
  }
  y_train <- (train$y_raw - outcome_mean) / outcome_sd

  # Map the 78 training schools to consecutive Stan indices.
  training_school_levels <- sort(unique(train$schoolid))
  id_train <- match(train$schoolid, training_school_levels)
  J_train <- length(training_school_levels)

  stan_data_fold <- list(
    N = nrow(train),
    J = J_train,
    Kx = ncol(X_train),
    Kz = ncol(Z_train),
    Kw = ncol(W_train),
    Kh = ncol(H_train),
    id = id_train,
    y = as.numeric(y_train),
    Tr = as.integer(train$Tr),
    X = X_train,
    Z = Z_train,
    W = W_train,
    H = H_train
  )

  # Covariates of the omitted school, transformed only with training-fold
  # constants. These are stored now for the later predictive exercise.
  h_heldout <- c(
    intercept = 1,
    as.numeric(apply_scaling(heldout[1L, , drop = FALSE], h_names, h_scaling))
  )
  names(h_heldout) <- c("intercept", h_names)
  zbar_heldout <- c(
    intercept = 1,
    free_lunch = mean(heldout$free_lunch)
  )

  # Direct randomized-experiment estimate for subsequent validation. Its
  # outcome data are not included in the fitted training model.
  treated_y <- heldout$y_raw[heldout$Tr == 1L]
  control_y <- heldout$y_raw[heldout$Tr == 0L]
  if (length(treated_y) < 2L || length(control_y) < 2L) {
    direct_estimate <- c(estimate = NA_real_, se = NA_real_)
  } else {
    estimate_raw <- mean(treated_y) - mean(control_y)
    variance_raw <- stats::var(treated_y) / length(treated_y) +
      stats::var(control_y) / length(control_y)
    direct_estimate <- c(
      estimate = estimate_raw / outcome_sd,
      se = sqrt(variance_raw) / outcome_sd
    )
  }

  list(
    stan_data = stan_data_fold,
    metadata = list(
      omitted_school = omitted_school,
      original_schoolid = unique(heldout$schoolid_raw),
      training_schools = training_school_levels,
      training_N = nrow(train),
      heldout_N = nrow(heldout),
      heldout_n_treated = length(treated_y),
      heldout_n_control = length(control_y),
      outcome_center = outcome_mean,
      outcome_scale = outcome_sd,
      age_scaling = age_scaling,
      W_scaling = w_scaling,
      W_columns_omitted_as_constant = setdiff(w_names, w_names_fold),
      H_scaling = h_scaling,
      h_heldout = h_heldout,
      zbar_heldout = zbar_heldout,
      direct_ATE = direct_estimate
    )
  )
}

summarise_fit_diagnostics <- function(fit, max_treedepth) {
  fit_summary <- rstan::summary(fit)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler, function(x) sum(x[, "divergent__"]), numeric(1)
  ))
  treedepth_hits <- sum(vapply(
    sampler, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)
  ))
  ebfmi <- vapply(sampler, function(x) {
    energy <- x[, "energy__"]
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))

  list(
    max_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
    min_n_eff = min(fit_summary[, "n_eff"], na.rm = TRUE),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    ebfmi_by_chain = ebfmi,
    min_ebfmi = min(ebfmi)
  )
}

atomic_save_rds <- function(object, path) {
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  saveRDS(object, temporary_path, compress = "xz")
  if (file.exists(path) && !file.remove(path)) {
    unlink(temporary_path)
    stop("Could not replace existing checkpoint: ", path)
  }
  if (!file.rename(temporary_path, path)) {
    unlink(temporary_path)
    stop("Could not save checkpoint: ", path)
  }
}

# -----------------------------------------------------------------------------
# Compile once and estimate all folds
# -----------------------------------------------------------------------------

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, min(loso_config$cores, parallel::detectCores())))

message("Compiling the Stan model once for all LOSO folds...")
compiled_model <- rstan::stan_model(model_code = stan_code)

school_indices <- sort(unique(analysis_data$stan_school))
if (!identical(school_indices, seq_len(J))) {
  stop("Full-sample school indices are not consecutive.")
}
if (!is.null(loso_config$schools)) {
  requested_schools <- as.integer(loso_config$schools)
  if (any(!requested_schools %in% school_indices)) {
    stop("loso_config$schools contains an invalid school index.")
  }
  school_indices <- requested_schools
}

for (omitted_school in school_indices) {
  result_path <- file.path(
    loso_config$output_dir,
    sprintf("loso_school_%02d.rds", omitted_school)
  )
  fit_path <- file.path(
    loso_config$output_dir,
    sprintf("loso_school_%02d_stanfit.rds", omitted_school)
  )

  if (file.exists(result_path) && !isTRUE(loso_config$overwrite)) {
    message("Skipping completed fold ", omitted_school, " of ", J)
    next
  }

  message("Starting fold ", omitted_school, " of ", J)
  fold <- make_loso_fold(
    omitted_school = omitted_school,
    data = analysis_data,
    w_names = w_available,
    h_names = h_requested
  )

  fold_seed <- loso_config$seed + omitted_school
  elapsed <- system.time({
    fit_fold <- rstan::sampling(
      object = compiled_model,
      data = fold$stan_data,
      chains = loso_config$chains,
      iter = loso_config$iter,
      warmup = loso_config$warmup,
      seed = fold_seed,
      refresh = loso_config$refresh,
      control = list(
        adapt_delta = loso_config$adapt_delta,
        max_treedepth = loso_config$max_treedepth
      )
    )
  })

  posterior <- rstan::extract(
    fit_fold,
    pars = c("Gamma", "Sigma_delta", "sd_delta", "Omega", "sigma")
  )
  diagnostics <- summarise_fit_diagnostics(
    fit_fold,
    max_treedepth = loso_config$max_treedepth
  )

  fold_result <- list(
    metadata = fold$metadata,
    posterior = posterior,
    diagnostics = diagnostics,
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    seed = fold_seed,
    model_config = loso_config,
    column_names = list(
      X = colnames(fold$stan_data$X),
      Z = colnames(fold$stan_data$Z),
      W = colnames(fold$stan_data$W),
      H = colnames(fold$stan_data$H)
    )
  )
  atomic_save_rds(fold_result, result_path)

  if (isTRUE(loso_config$save_full_stanfit)) {
    atomic_save_rds(fit_fold, fit_path)
  }

  message(
    "Completed fold ", omitted_school, " of ", J,
    " in ", round(elapsed[["elapsed"]] / 60, 1), " minutes; ",
    "divergences = ", diagnostics$divergences,
    ", max R-hat = ", round(diagnostics$max_rhat, 3)
  )
  rm(fit_fold, fold_result, posterior)
  invisible(gc())
}

message("All LOSO fits are complete: ", normalizePath(loso_config$output_dir))
