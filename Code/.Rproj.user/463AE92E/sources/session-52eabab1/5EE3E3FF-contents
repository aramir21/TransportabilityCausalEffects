rm(list = ls())
set.seed(10101)

# --- Basic structure ---
J  <- 40 # 40, 30, 20, 10 number of studies
l <- 8 # 8, 6, 4, 2
N  <- c(rep(500, l), rep(256, l), rep(145, l), rep(760, l), rep(678, l)) # sample sizes per study
NJ <- sum(N)

study_id <- unlist(lapply(1:J, function(j) rep(j, N[j])))

# Individual-level X
x1 <- rnorm(NJ)
x2 <- rnorm(NJ)
x3 <- rnorm(NJ)
X  <- cbind(x1, x2, x3)                    # Kx = 3
Kx <- ncol(X)

# Z for heterogeneous treatment effect: intercept + x1 only
Z  <- cbind(1, x1)                         # Kz = 2
Kz <- ncol(Z)

# Study-level covariates W_j   (Kw = 3 here)
Kw <- 3
W_study <- matrix(NA, nrow = J, ncol = Kw)
for (j in 1:J) {
  W_study[j, ] <- c(1, rnorm(Kw - 1, mean = j, sd = 1))
}

# Expand W to individual level (so we can store it in Data)
W_full <- W_study[study_id, ]              # NJ x Kw

# Treatment assignment
Tr <- rbinom(NJ, 1, 0.3)

# --- True hierarchical parameters ---

beta_true  <- c(0.5, 0.2, 0.3)           # baseline effect of X
gamma_true <- c(0.4, 0.1, 0.2)           # baseline effect of W on y

# Psi: Kw x Kz mapping W_j -> mean(delta_j)
Psi_true  <- matrix(c( 0.8,  0.2,
                       0.3,  0.4,
                       0.1, 0.5),
                    nrow = Kw, byrow = TRUE)

Sigma_delta_true <- matrix(c(0.15, 0.02,
                             0.02, 0.10), 2, 2)
sigma_true <- 0.5

# --- Draw study-specific deltas and generate data ---

delta_true <- matrix(NA, nrow = J, ncol = Kz)
theta      <- numeric(NJ)

library(MASS)

start <- 1
for (j in 1:J) {
  N_j <- N[j]
  ind <- start:(start + N_j - 1)
  
  # mean(delta_j | W_j)
  mean_delta_j <- W_study[j, ] %*% Psi_true   # 1 x Kz
  delta_j      <- MASS::mvrnorm(1, mean_delta_j, Sigma_delta_true)
  
  delta_true[j, ] <- delta_j
  
  theta[ind] <- Z[ind, ] %*% delta_j          # theta_ij
  start <- start + N_j
}

# Outcome: includes X*beta and W*gamma
e  <- rnorm(NJ, 0, sigma_true)
y  <- Tr * theta + X %*% beta_true + W_full %*% gamma_true + e

# Put everything in a data.frame including w's
Data <- data.frame(
  id = study_id,
  y  = as.numeric(y),
  x1 = x1, x2 = x2, x3 = x3,
  w1 = W_full[, 1],
  w2 = W_full[, 2],
  w3 = W_full[, 3],
  Tr = Tr
)

stan_hier_meta_W_gamma <- "
data {
  int<lower=1> N;                 // total observations
  int<lower=1> J;                 // number of studies
  int<lower=1> Kx;                // # individual X covariates
  int<lower=1> Kz;                // dim of z_ij (treat. effect)
  int<lower=1> Kw;                // # study-level W covariates

  int<lower=1,upper=J> id[N];     // study id
  vector[N] y;                    // outcome
  int<lower=0,upper=1> Tr[N];     // treatment indicator

  matrix[N, Kx] X;                // individual-level X
  matrix[N, Kz] Z;                // Z for theta = Z * delta_j
  matrix[J, Kw] W;                // study-level covariates
}

parameters {
  // Baseline effects
  vector[Kx] beta;                // effect of X on y
  vector[Kw] gamma;               // effect of W on y (fixed, same across j)

  // Mapping from W_j to mean(delta_j): Psi is Kw x Kz
  matrix[Kw, Kz] Psi;

  // Between-study covariance for delta_j
  vector<lower=0>[Kz] tau;              // scales
  cholesky_factor_corr[Kz] L_Omega;     // correlation

  matrix[J, Kz] delta_raw;              // non-centered
  real<lower=0> sigma;                  // residual sd
}

transformed parameters {
  matrix[J, Kz] delta;
  matrix[J, Kz] delta_mean;
  matrix[Kz, Kz] L_Sigma;

  L_Sigma = diag_pre_multiply(tau, L_Omega);
  delta_mean = W * Psi;

  for (j in 1:J) {
    // Convert row_vectors to column_vectors before adding
    delta[j] = ( (delta_mean[j])' + L_Sigma * (delta_raw[j])' )';
  }
}

model {
  // Priors
  beta       ~ normal(0, 2);
  gamma      ~ normal(0, 2);           // W effect on y
  to_vector(Psi) ~ normal(0, 1);

  tau        ~ cauchy(0, 2.5);
  L_Omega    ~ lkj_corr_cholesky(2.0);
  to_vector(delta_raw) ~ normal(0, 1);

  sigma ~ exponential(1);

  // Likelihood
  for (n in 1:N) {
    int j = id[n];
    real theta_n = dot_product(Z[n], delta[j]);   // theta_ij
    real mu_n    = Tr[n] * theta_n
                   + dot_product(X[n], beta)
                   + dot_product(W[j], gamma);    // W_j * gamma
    y[n] ~ normal(mu_n, sigma);
  }
}

generated quantities {
  cov_matrix[Kz] Sigma_delta;
  Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma);
}
"
library(rstan)

stan_data <- list(
  N  = NJ,
  J  = J,
  Kx = Kx,
  Kz = Kz,
  Kw = Kw,
  id = study_id,
  y  = as.vector(y),
  Tr = Tr,
  X  = X,
  Z  = Z,
  W  = W_study          # J x Kw
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(
  model_code = stan_hier_meta_W_gamma,
  data   = stan_data,
  chains = 1,
  iter   = 2000,
  warmup = 1000,
  seed   = 1234
)

save(fit, file = "FitV2J40.RData")
PopVal <- list(delta_true = delta_true, theta = theta, X = X, Z = Z, W = W_study)
save(PopVal, file = "PopValJ40.RData")

### Data for second example in the paper ###
HtH_inv <- solve(t(W_study)%*%W_study)
eigen(HtH_inv)
hj0 <- c(1, 0, 0)
t(hj0)%*%HtH_inv%*%hj0
hj0 <- c(0, 1, 0)
t(hj0)%*%HtH_inv%*%hj0
hj0 <- c(0, 0, 1)
t(hj0)%*%HtH_inv%*%hj0
eigen(Sigma_delta_true)
zij0 <- c(1, 0)
t(zij0)%*%Sigma_delta_true%*%zij0
zij0 <- c(0, 1)
t(zij0)%*%Sigma_delta_true%*%zij0
# post <- summary(fit, pars = c("beta", "gamma", "Psi", "Sigma_delta", "sigma"),
#                 probs = c(0.05, 0.5, 0.95))$summary
# # Add true values
# post_with_true <- cbind(true = c(beta_true, gamma_true, as.vector(t(Psi_true)),
#                                  as.vector(Sigma_delta_true), sigma_true), post)
# 
# post_with_true
# 
# 
# rstan::stan_trace(fit, pars = "beta")
# rstan::stan_trace(fit, pars = "gamma")
# rstan::stan_trace(fit, pars = "Psi")
# rstan::stan_trace(fit, pars = "Sigma_delta")
# rstan::stan_trace(fit, pars = "sigma")
# 
# post_delta <- rstan::extract(fit, pars = "delta")$delta
# dim(post_delta)
# j <- 1
# delta_j_draws <- post_delta[, j, ]  # matrix: iterations x Kz
# dim(delta_j_draws)
# apply(delta_j_draws, 2, mean)      # posterior means of δ_2 components
# apply(delta_j_draws, 2, sd)        # posterior s.d.
# plot(density(delta_j_draws[, 1]),
#      main = "Posterior of delta[2,1]",
#      xlab = "delta[2,1]")
# abline(v = mean(delta_j_draws[, 1]), col = 2, lwd = 2)
# n <- 10
# j_n <- stan_data$id[n]
# z_n <- stan_data$Z[n, ]             # vector of length Kz
# 
# theta_n_draws <- post_delta[, j_n, ] %*% z_n  # iterations x 1
# plot(density(theta_n_draws),
#      main = paste0("Posterior of theta for obs ", n),
#      xlab = "theta_n")
# theta[n]
# summary(coda::mcmc(theta_n_draws))
