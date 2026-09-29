library(dplyr)
library(ggplot2)

# -----------------------------
# Settings
# -----------------------------
j <- 21  # Study

# -----------------------------
# Load data and fitted model
# -----------------------------
load("PopValJ40.RData")
load("FitV2J40.RData")

J40 <- fit
draws40 <- rstan::extract(J40, permuted = TRUE)

# -----------------------------
# Reconstruct study ID
# -----------------------------
J <- 40
l <- J / 5

N <- c(
  rep(35, l),
  rep(60, l),
  rep(85, l),
  rep(110, l),
  rep(135, l)
)

study_id <- rep(1:J, times = N)

# -----------------------------
# Individuals belonging to study j
# -----------------------------
idx_j <- which(study_id == j)

# -----------------------------
# Mean of Z variable in study j
# -----------------------------
mean_Z_j <- colMeans(
  PopVal[["Z"]][idx_j, , drop = FALSE]
)

# -----------------------------
# tau true study j
# -----------------------------
true_tau_j <- t(mean_Z_j) %*% PopVal[["delta_true"]][j,]
# -----------------------------
# Posterior prediction using Z40
# -----------------------------
post_delta40 <- draws40[["delta"]][, j, ]

tau40 <- as.vector(
  post_delta40 %*% mean_Z_j
)

# Posterior summaries
tau_mean <- mean(tau40)

ci_low <- as.numeric(
  quantile(tau40, 0.025)
)

ci_high <- as.numeric(
  quantile(tau40, 0.975)
)

ci_length <- ci_high - ci_low # Interval lenght

coverage <- ifelse(true_tau_j <= ci_high & true_tau_j >= ci_low, 1, 0) # Coverage == 1
  
# Data frame with posterior draws
df_tau <- tibble(
  value = tau40
)

# Histogram of posterior draws of tau
ggplot(df_tau, aes(x = value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    fill = "grey80",
    color = "white",
    linewidth = 0.3
  ) +
  
  # Posterior mean
  geom_vline(
    aes(
      xintercept = tau_mean,
      linetype = "Posterior mean"
    ),
    linewidth = 0.9,
    color = "black"
  ) +
  
  # Lower bound of 95% credible interval
  geom_vline(
    aes(
      xintercept = ci_low,
      linetype = "95% credible interval"
    ),
    linewidth = 0.8,
    color = "black"
  ) +
  
  # Upper bound of 95% credible interval
  geom_vline(
    aes(
      xintercept = ci_high,
      linetype = "95% credible interval"
    ),
    linewidth = 0.8,
    color = "black"
  ) +
  
  # True population value
  geom_vline(
    aes(
      xintercept = true_tau_j,
      linetype = "Population value"
    ),
    linewidth = 1.0,
    color = "black"
  ) +
  
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Posterior mean" = "solid",
      "95% credible interval" = "dashed",
      "Population value" = "dotted"
    )
  ) +
  
  labs(
    x = expression(tau[j]),
    y = "Density"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "plain")
  )
