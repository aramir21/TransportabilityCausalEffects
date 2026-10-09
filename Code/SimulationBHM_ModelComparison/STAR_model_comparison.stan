// Aggregate likelihood with known school-specific sampling variances.
// Integrate out latent school ATEs for efficient HMC.
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0, upper=1> composition_variance;
  matrix[N, K] X;
  vector[N] estimate;
  vector<lower=0>[N] sampling_se;
  vector[N] p;
  row_vector[K] x_new;
  real p_new;
}
parameters {
  vector[K] coefficients;
  vector<lower=0>[1 - composition_variance] common_sd;
  vector<lower=0>[2 * composition_variance] sd_delta;
  cholesky_factor_corr[1 + composition_variance] L_Omega;
}
transformed parameters {
  vector[N] residual_variance;
  matrix[1 + composition_variance, 1 + composition_variance] Sigma;
  if (composition_variance == 1) {
    matrix[2, 2] L_Sigma = diag_pre_multiply(sd_delta, L_Omega);
    Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    // z' Sigma z = ||z' L_Sigma||^2 avoids cancellation in the quadratic.
    for (j in 1:N)
      residual_variance[j] = square(L_Sigma[1,1] + p[j] * L_Sigma[2,1])
                            + square(p[j] * L_Sigma[2,2]);
  } else {
    Sigma[1,1] = square(common_sd[1]);
    residual_variance = rep_vector(Sigma[1,1], N);
  }
}
model {
  coefficients ~ normal(0, 0.5);
  common_sd ~ normal(0, 0.5);
  sd_delta ~ normal(0, 0.5);
  L_Omega ~ lkj_corr_cholesky(2);
  estimate ~ normal(X * coefficients, sqrt(square(sampling_se) + residual_variance));
}
generated quantities {
  real predictive_location = x_new * coefficients;
  real predictive_sd;
  if (composition_variance == 1) {
    matrix[2, 2] L_Sigma = diag_pre_multiply(sd_delta, L_Omega);
    predictive_sd = sqrt(square(L_Sigma[1,1] + p_new * L_Sigma[2,1])
                         + square(p_new * L_Sigma[2,2]));
  } else
    predictive_sd = common_sd[1];
}
