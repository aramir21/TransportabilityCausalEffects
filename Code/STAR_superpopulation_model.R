# Shared model and utilities for the STAR super-population experiments.
# The likelihood is Gaussian and the regressors match the theoretical DGP:
# X = (1, female, age), Z = (1, free_lunch), W = avg_experience,
# H = (1, H1, H2).

star_normal_stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> Kx;
  int<lower=1> Kz;
  int<lower=1> Kw;
  int<lower=1> Kh;
  int<lower=1, upper=J> id[N];
  vector[N] y;
  int<lower=0, upper=1> Tr[N];
  matrix[N, Kx] X;
  matrix[N, Kz] Z;
  matrix[N, Kw] W;
  matrix[J, Kh] H_school;
}
parameters {
  vector[Kx] beta;
  vector[Kw] gamma;
  matrix[Kh, Kz] Gamma;
  vector<lower=0>[Kz] sd_delta;
  cholesky_factor_corr[Kz] L_Omega;
  matrix[J, Kz] delta_raw;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[Kz, Kz] L_Sigma = diag_pre_multiply(sd_delta, L_Omega);
  matrix[J, Kz] delta_mean = H_school * Gamma;
  matrix[J, Kz] delta;
  for (j in 1:J)
    delta[j] = (delta_mean[j]' + L_Sigma * delta_raw[j]')';
}
model {
  beta[1] ~ normal(0, 1.5);
  beta[2:Kx] ~ normal(0, 1);
  gamma ~ normal(0, 0.5);
  to_vector(Gamma) ~ normal(0, 0.5);
  sd_delta ~ normal(0, 0.5);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(delta_raw) ~ std_normal();
  sigma ~ normal(0, 1);

  for (n in 1:N) {
    real treatment_effect = dot_product(Z[n], delta[id[n]]);
    real mu = dot_product(X[n], beta) + dot_product(W[n], gamma)
              + Tr[n] * treatment_effect;
    y[n] ~ normal(mu, sigma);
  }
}
generated quantities {
  corr_matrix[Kz] Omega = multiply_lower_tri_self_transpose(L_Omega);
  cov_matrix[Kz] Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma);
}
"

make_star_stan_input <- function(data, school_truth, omitted_school = NULL) {
  required_data <- c("sim_school_id", "y", "treatment", "female", "age",
                     "free_lunch", "W", "H1", "H2", "pF")
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data)) {
    stop("Scenario data are missing: ", paste(missing_data, collapse = ", "))
  }
  if (!is.null(omitted_school)) {
    omitted_school <- as.integer(omitted_school)
    heldout <- data[data$sim_school_id == omitted_school, , drop = FALSE]
    train <- data[data$sim_school_id != omitted_school, , drop = FALSE]
    truth_heldout <- school_truth[
      school_truth$sim_school_id == omitted_school, , drop = FALSE]
    if (!nrow(heldout) || !nrow(train) || nrow(truth_heldout) != 1L) {
      stop("Invalid omitted school: ", omitted_school)
    }
  } else {
    train <- data
    heldout <- NULL
    truth_heldout <- NULL
  }

  school_levels <- sort(unique(train$sim_school_id))
  id <- match(train$sim_school_id, school_levels)
  first <- match(school_levels, train$sim_school_id)
  H_school <- cbind(intercept = 1, H1 = train$H1[first], H2 = train$H2[first])
  X <- cbind(intercept = 1, female = train$female, age = train$age)
  Z <- cbind(intercept = 1, free_lunch = train$free_lunch)
  W <- matrix(train$W, ncol = 1L, dimnames = list(NULL, "avg_experience"))

  stan_data <- list(
    N = nrow(train), J = length(school_levels),
    Kx = ncol(X), Kz = ncol(Z), Kw = ncol(W), Kh = ncol(H_school),
    id = as.integer(id), y = as.numeric(train$y),
    Tr = as.integer(train$treatment), X = X, Z = Z, W = W,
    H_school = H_school
  )
  stopifnot(
    identical(sort(unique(stan_data$id)), seq_len(stan_data$J)),
    all(stan_data$Tr %in% 0:1), !anyNA(unlist(stan_data))
  )

  metadata <- list(
    omitted_school = omitted_school,
    training_schools = school_levels,
    training_N = nrow(train),
    column_names = list(X = colnames(X), Z = colnames(Z), W = colnames(W),
                        H = colnames(H_school))
  )
  if (!is.null(omitted_school)) {
    metadata$heldout_N <- nrow(heldout)
    metadata$h_heldout <- c(intercept = 1, H1 = heldout$H1[1L],
                            H2 = heldout$H2[1L])
    # Primary prediction composition is population pF, not realized sample pF.
    metadata$zbar_heldout_pop <- c(intercept = 1, pF = heldout$pF[1L])
    metadata$zbar_heldout_finite <- c(intercept = 1,
                                      free_lunch = mean(heldout$free_lunch))
    metadata$tau_pop_true <- truth_heldout$tau_pop_true
    metadata$tau_finite_true <- truth_heldout$tau_finite_true
  }
  list(stan_data = stan_data, metadata = metadata)
}

summarise_star_fit_diagnostics <- function(fit, max_treedepth) {
  fit_summary <- rstan::summary(fit)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]),
                            numeric(1)))
  treedepth_hits <- sum(vapply(
    sampler, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))
  ebfmi <- vapply(sampler, function(x) {
    energy <- x[, "energy__"]
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))
  list(max_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
       min_n_eff = min(fit_summary[, "n_eff"], na.rm = TRUE),
       divergences = divergences, treedepth_hits = treedepth_hits,
       ebfmi_by_chain = ebfmi, min_ebfmi = min(ebfmi))
}

atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), "_"), dirname(path), ".tmp")
  saveRDS(object, temporary, compress = compress)
  if (file.exists(path) && !file.remove(path)) {
    unlink(temporary)
    stop("Could not replace: ", path)
  }
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save: ", path)
  }
}
