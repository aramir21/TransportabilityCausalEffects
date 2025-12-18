set.seed(010101)
J <- 5                          # number of studies
N <- c(500, 256, 145, 760, 678) # sample sizes per study
NJ <- sum(N)

# study id for each individual
study_id <- unlist(sapply(1:J, function(l) rep(l, N[l])))

# individual-level covariates
x1 <- rnorm(NJ); x2 <- rnorm(NJ); x3 <- rnorm(NJ)
X  <- cbind(x1, x2, x3)

# study-level covariates (same within study)
w1 <- unlist(sapply(1:J, function(l) rep(rnorm(1, l, 1), N[l])))
w2 <- unlist(sapply(1:J, function(l) rep(rnorm(1, l, 1), N[l])))
w3 <- unlist(sapply(1:J, function(l) rep(rnorm(1, l, 1), N[l])))
w4 <- unlist(sapply(1:J, function(l) rep(rnorm(1, l, 1), N[l])))
w5 <- unlist(sapply(1:J, function(l) rep(rnorm(1, l, 1), N[l])))
W  <- cbind(w1, w2, w3, w4, w5)

# treatment indicator
Tr <- rbinom(NJ, 1, 0.3)

# z_ij used for heterogeneous treatment effect
Z <- cbind(1, x1, w1, w2)

phi   <- c(1, 0.5, 0.6, 0.3)
SIGMA <- 0.05 * diag(ncol(Z))
PHI   <- 0.05 * diag(ncol(Z))
OMEGA <- SIGMA + PHI

theta <- numeric(NJ)
delta <- matrix(NA, nrow = J, ncol = ncol(Z))

start <- 1
for (jj in 1:J) {
  deltaj <- MASS::mvrnorm(1, phi, OMEGA) # marginal N(phi, SIGMA+PHI)
  ind    <- start:(start + N[jj] - 1)
  
  theta[ind] <- Z[ind, ] %*% deltaj
  delta[jj, ] <- deltaj
  
  start <- start + N[jj]
}

mean(theta)
plot(density(theta))

# baseline effects
beta  <- rep(0.6, 3)
gamma <- c(0.5, 1, 0.3, 0.4, -0.5)

e <- rnorm(NJ)
y <- Tr * theta + X %*% beta + W %*% gamma + e

Data <- data.frame(
  id = study_id,
  y  = y,
  x1 = x1, x2 = x2, x3 = x3,
  w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5,
  Tr = Tr
)


stan_hier_micro_model <- "
data {
  int<lower=1> N;            // total observations
  int<lower=1> J;            // number of studies
  int<lower=1> Kx;           // # individual X covariates
  int<lower=1> Kw;           // # study-level W covariates
  int<lower=1> Kz;           // dimension of z_ij (heterogeneous treatment effects)

  int<lower=1,upper=J> id[N];   // study id
  vector[N] y;                   // outcome
  int<lower=0,upper=1> Tr[N];   // treatment indicator

  matrix[N, Kx] X;              // individual-level X
  matrix[N, Kw] W;              // study-level W (repeated rows)
  matrix[N, Kz] Z;              // Z matrix for theta_ij = Z * delta_j
}

parameters {
  // Fixed effects
  vector[Kx] beta;              // baseline effect of X on y
  vector[Kw] gamma;             // baseline effect of W on y

  // Hierarchical treatment-effect parameters
  vector[Kz] phi;               // global mean of delta_j
  vector<lower=0>[Kz] tau;      // scales for delta_j
  cholesky_factor_corr[Kz] L_Omega;   // correlation of delta_j

  matrix[J, Kz] delta_raw;      // non-centered delta_j (each row j)
  real<lower=0> sigma;          // residual sd
}

transformed parameters {
  matrix[J, Kz] delta;          // delta_j for each study (row j)
  matrix[Kz, Kz] L_Sigma;

  // Build Cholesky of Sigma_delta
  L_Sigma = diag_pre_multiply(tau, L_Omega);

  // Non-centered parameterization
  for (j in 1:J) {
    // delta[j] is a row_vector; RHS is vector -> transpose to row_vector
    delta[j] = (phi + L_Sigma * (delta_raw[j])')';
  }
}

model {
  // Priors
  beta  ~ normal(0, 5);
  gamma ~ normal(0, 5);
  phi   ~ normal(0, 5);

  tau        ~ cauchy(0, 2.5);
  L_Omega    ~ lkj_corr_cholesky(2.0);
  to_vector(delta_raw) ~ normal(0, 1);

  sigma ~ exponential(1);

  // Likelihood
  for (n in 1:N) {
    int j = id[n];
    real theta_n = dot_product(Z[n], delta[j]);     // heterogeneous treatment effect
    real mu_n = Tr[n] * theta_n
                + dot_product(X[n], beta)
                + dot_product(W[n], gamma);
    y[n] ~ normal(mu_n, sigma);
  }
}

generated quantities {
  cov_matrix[Kz] Sigma_delta;
  Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma);
}
"
library(rstan)

# Assuming you already ran your simulation and have:
# Data with columns: id, y, x1, x2, x3, w1, w2, w3, w4, w5, Tr

N  <- nrow(Data)
J  <- length(unique(Data$id))
id <- Data$id

# Individual-level X (all x's)
X <- as.matrix(Data[, c("x1", "x2", "x3")])
Kx <- ncol(X)

# Study-level W (all w's)
W <- as.matrix(Data[, c("w1", "w2", "w3", "w4", "w5")])
Kw <- ncol(W)

# z_ij used for heterogeneous treatment effect (you chose: 1, x1, w1, w2)
Z <- cbind(1, Data$x1, Data$w1, Data$w2)
Kz <- ncol(Z)

y  <- Data$y
Tr <- Data$Tr

stan_data <- list(
  N  = N,
  J  = J,
  Kx = Kx,
  Kw = Kw,
  Kz = Kz,
  id = id,
  y  = y,
  Tr = Tr,
  X  = X,
  W  = W,
  Z  = Z
)


fit <- stan(
  model_code = stan_hier_micro_model,
  data = stan_data,          # Your list prepared earlier
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 1234
)

print(fit, pars = c("phi", "beta", "gamma", "Sigma_delta", "sigma"))
rstan::stan_trace(fit, pars = "beta")
rstan::stan_trace(fit, pars = "gamma")
rstan::stan_trace(fit, pars = "phi")
rstan::stan_trace(fit, pars = "Sigma_delta")

post_delta <- rstan::extract(fit, pars = "delta")$delta
dim(post_delta)
j <- 1
delta_j_draws <- post_delta[, j, ]  # matrix: iterations x Kz
dim(delta_j_draws)
apply(delta_j_draws, 2, mean)      # posterior means of δ_2 components
apply(delta_j_draws, 2, sd)        # posterior s.d.
plot(density(delta_j_draws[, 1]),
     main = "Posterior of delta[2,1]",
     xlab = "delta[2,1]")
abline(v = mean(delta_j_draws[, 1]), col = 2, lwd = 2)
n <- 10
j_n <- stan_data$id[n]
z_n <- stan_data$Z[n, ]             # vector of length Kz

theta_n_draws <- post_delta[, j_n, ] %*% z_n  # iterations x 1
plot(density(theta_n_draws),
     main = paste0("Posterior of theta for obs ", n),
     xlab = "theta_n")
theta[n]
summary(coda::mcmc(theta_n_draws))
