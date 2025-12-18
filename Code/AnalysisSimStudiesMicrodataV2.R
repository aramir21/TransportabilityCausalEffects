set.seed(12345)

load("FitV2J40.RData")
J40 <- fit
load("FitV2J30.RData")
J30 <- fit
load("FitV2J20.RData")
J20 <- fit
load("FitV2J10.RData")
J10 <- fit
rm(fit)

draws40 <- rstan::extract(J40, permuted = TRUE)
Psi40_22 <- draws40[["Psi"]][,2,2]
draws30 <- rstan::extract(J30, permuted = TRUE)
Psi30_22 <- draws30[["Psi"]][,2,2]
draws20 <- rstan::extract(J20, permuted = TRUE)
Psi20_22 <- draws20[["Psi"]][,2,2]
draws10 <- rstan::extract(J10, permuted = TRUE)
Psi10_22 <- draws10[["Psi"]][,2,2]

# Plot Gamma_22
library(dplyr)
library(ggplot2)
library(tidyr)

truePsi_22 <- 0.4

df <- bind_rows(
  tibble(value = Psi40_22, J = "J = 40"),
  tibble(value = Psi30_22, J = "J = 30"),
  tibble(value = Psi20_22, J = "J = 20"),
  tibble(value = Psi10_22, J = "J = 10")
) %>%
  mutate(J = factor(J, levels = c("J = 40","J = 30","J = 20","J = 10")))

# 95% equal-tailed (symmetric in probability mass) credible interval per panel
ci_df <- df %>%
  group_by(J) %>%
  summarise(
    ci_low  = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    .groups = "drop"
  )

library(ggplot2)

# ---- optional: enforce consistent axes across panels ----
x_lim <- range(df$value, truePsi_22, ci_df$ci_low, ci_df$ci_high, na.rm = TRUE)

# If you want a symmetric x-window around the true value (often nicer):
# halfw <- max(abs(df$value - truePsi_22), abs(ci_df$ci_low - truePsi_22), abs(ci_df$ci_high - truePsi_22), na.rm = TRUE)
# x_lim <- truePsi_22 + c(-halfw, halfw)

# Set a common y-limit across panels (so comparisons are valid)
# Build the histogram counts per panel to get a safe upper bound:
tmp_counts <- ggplot_build(
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 80) +
    facet_wrap(~ J, ncol = 2)
)$data[[1]]$count
y_max <- max(tmp_counts, na.rm = TRUE)

ggplot(df, aes(x = value)) +
  geom_histogram(
    bins = 80,
    boundary = 0,        # stable bin alignment
    closed = "left",
    fill = "grey85",     # grayscale-safe
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    aes(xintercept = truePsi_22, linetype = "True value"),
    linewidth = 1.0,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_low, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_high, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  facet_wrap(~ J, ncol = 2) +
  scale_linetype_manual(
    name = NULL,
    values = c("True value" = "solid", "95% CI" = "dashed")
  ) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
  labs(
    # title = expression("Posterior draws of " * Gamma[22]),
    x = expression(Gamma[22]),
    y = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "plain")
  )

# Plot Sigma_22
Sigma40_22 <- draws40[["Sigma_delta"]][,2,2]
Sigma30_22 <- draws30[["Sigma_delta"]][,2,2]
Sigma20_22 <- draws20[["Sigma_delta"]][,2,2]
Sigma10_22 <- draws10[["Sigma_delta"]][,2,2]

trueSigma_22 <- 0.1

df <- bind_rows(
  tibble(value = Sigma40_22, J = "J = 40"),
  tibble(value = Sigma30_22, J = "J = 30"),
  tibble(value = Sigma20_22, J = "J = 20"),
  tibble(value = Sigma10_22, J = "J = 10")
) %>%
  mutate(J = factor(J, levels = c("J = 40","J = 30","J = 20","J = 10")))

# 95% equal-tailed (symmetric in probability mass) credible interval per panel
ci_df <- df %>%
  group_by(J) %>%
  summarise(
    ci_low  = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    .groups = "drop"
  )

# ---- optional: enforce consistent axes across panels ----
# x_lim <- range(df$value, trueSigma_22, ci_df$ci_low, ci_df$ci_high, na.rm = TRUE)
x_lim <- c(0, 0.75)
# If you want a symmetric x-window around the true value (often nicer):
# halfw <- max(abs(df$value - trueSigma_22), abs(ci_df$ci_low - trueSigma_22), abs(ci_df$ci_high - trueSigma_22), na.rm = TRUE)
# x_lim <- trueSigma_22 + c(-halfw, halfw)

# Set a common y-limit across panels (so comparisons are valid)
# Build the histogram counts per panel to get a safe upper bound:
tmp_counts <- ggplot_build(
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 80) +
    facet_wrap(~ J, ncol = 2)
)$data[[1]]$count
y_max <- max(tmp_counts, na.rm = TRUE)

ggplot(df, aes(x = value)) +
  geom_histogram(
    bins = 80,
    boundary = 0,        # stable bin alignment
    closed = "left",
    fill = "grey85",     # grayscale-safe
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    aes(xintercept = trueSigma_22, linetype = "True value"),
    linewidth = 1.0,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_low, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_high, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  facet_wrap(~ J, ncol = 2) +
  scale_linetype_manual(
    name = NULL,
    values = c("True value" = "solid", "95% CI" = "dashed")
  ) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
  labs(
    # title = expression("Posterior draws of " * Sigma[22]),
    x = expression(Sigma[22]),
    y = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "plain")
  )


# Plot theta
# --- posterior draws (delta[2,2]) ---
j <- 1 # Study
l <- 2 # Parameter
i <- 1 # Unit
post_delta40 <- draws40[["delta"]][,j,l]
load("PopValJ40.RData")
true_delta40 <- PopVal[["delta_true"]][j,l]
true_tau40 <- PopVal[["theta"]][i]

post_delta30 <- draws30[["delta"]][,j,l]
load("PopValJ30.RData")
true_delta30 <- PopVal[["delta_true"]][j,l]
true_tau30 <- PopVal[["theta"]][i]

post_delta20 <- draws20[["delta"]][,j,l]
load("PopValJ20.RData")
true_delta20 <- PopVal[["delta_true"]][j,l]
true_tau20 <- PopVal[["theta"]][i]

post_delta10 <- draws10[["delta"]][,j,l]
load("PopValJ10.RData")
true_delta10 <- PopVal[["delta_true"]][j,l]
true_tau10 <- PopVal[["theta"]][i]

# --- long df with panel-specific true values ---
df <- bind_rows(
  tibble(value = post_delta40, J = "J = 40", true = true_delta40),
  tibble(value = post_delta30, J = "J = 30", true = true_delta30),
  tibble(value = post_delta20, J = "J = 20", true = true_delta20),
  tibble(value = post_delta10, J = "J = 10", true = true_delta10)
) %>%
  mutate(J = factor(J, levels = c("J = 40","J = 30","J = 20","J = 10")))

# 95% equal-tailed CI per panel
ci_df <- df %>%
  group_by(J) %>%
  summarise(
    ci_low  = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    true    = first(true),
    .groups = "drop"
  )

# consistent axes across panels
x_lim <- range(df$value, ci_df$true, ci_df$ci_low, ci_df$ci_high, na.rm = TRUE)

tmp_counts <- ggplot_build(
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 80, boundary = 0, closed = "left") +
    facet_wrap(~ J, ncol = 2)
)$data[[1]]$count
y_max <- max(tmp_counts, na.rm = TRUE)

# ---- plot ----
ggplot(df, aes(x = value)) +
  geom_histogram(
    bins = 80,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = true, linetype = "True value"),
    linewidth = 1.0,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_low, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_high, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  facet_wrap(~ J, ncol = 2) +
  scale_linetype_manual(
    name = NULL,
    values = c("True value" = "solid", "95% CI" = "dashed")
  ) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
  labs(
    # title = expression("Posterior draws of " * theta[j2]),
    x = expression(theta[12]),
    y = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "plain")
  )

# Plot tau (treatment effect)
load("PopValJ40.RData")
true_tau40 <- PopVal[["theta"]][i]
post_delta40 <- draws40[["delta"]][,j,]
Z40 <- c(1, PopVal[["X"]][i,1])
tau40 <- post_delta40 %*% Z40

load("PopValJ30.RData")
true_tau30 <- PopVal[["theta"]][i]
post_delta30 <- draws30[["delta"]][,j,]
Z30 <- c(1, PopVal[["X"]][i,1])
tau30 <- post_delta30 %*% Z30

load("PopValJ20.RData")
true_tau20 <- PopVal[["theta"]][i]
post_delta20 <- draws20[["delta"]][,j,]
Z20 <- c(1, PopVal[["X"]][i,1])
tau20 <- post_delta20 %*% Z20

load("PopValJ10.RData")
true_tau10 <- PopVal[["theta"]][i]
post_delta10 <- draws10[["delta"]][,j,]
Z10 <- c(1, PopVal[["X"]][i,1])
tau10 <- post_delta10 %*% Z10

# --- long df (panel-specific true values) ---
df <- bind_rows(
  tibble(value = tau40, J = "J = 40", true = true_tau40),
  tibble(value = tau30, J = "J = 30", true = true_tau30),
  tibble(value = tau20, J = "J = 20", true = true_tau20),
  tibble(value = tau10, J = "J = 10", true = true_tau10)
) %>%
  mutate(J = factor(J, levels = c("J = 40","J = 30","J = 20","J = 10")))

# 95% equal-tailed CI per panel
ci_df <- df %>%
  group_by(J) %>%
  summarise(
    ci_low  = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    true    = first(true),
    .groups = "drop"
  )

# consistent axes across panels
x_lim <- range(df$value, ci_df$true, ci_df$ci_low, ci_df$ci_high, na.rm = TRUE)

tmp_counts <- ggplot_build(
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 80, boundary = 0, closed = "left") +
    facet_wrap(~ J, ncol = 2)
)$data[[1]]$count
y_max <- max(tmp_counts, na.rm = TRUE)

# ---- plot ----
ggplot(df, aes(x = value)) +
  geom_histogram(
    bins = 80,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = true, linetype = "True value"),
    linewidth = 1.0,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_low, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_high, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  facet_wrap(~ J, ncol = 2) +
  scale_linetype_manual(
    name = NULL,
    values = c("True value" = "solid", "95% CI" = "dashed")
  ) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
  labs(
    # title = bquote("Posterior draws of " * tau[.(i) * "," * .(j)]),
    x = bquote(tau[.(i) * "" * .(j)]),
    y = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "plain")
  )

library(e1071)

skewness(tau40)
skewness(tau30)
skewness(tau20)
skewness(tau10)

#### Prediction ####
x10 <- rnorm(1, 0, 1)
w10 <- rnorm(1, 0, 1); w20 <- rnorm(1, 0, 1)

Pred_tau <- function(x10, w10, w20, Gammas, Sigmas){
  # Gammas <- draws10[["Psi"]][1,,]
  # Sigmas <- draws10[["Sigma_delta"]][1,,]
  zij0 <- c(1, x10)
  hj0 <- c(1, w10, w20)
  mean_theta <- t(Gammas)%*%hj0
  thetaj0 <- MASS::mvrnorm(1, mu = mean_theta, Sigmas)
  tauj0 <-  zij0%*%thetaj0
  return(tauj0)
}

Kw <- 3
Psi_true  <- matrix(c( 0.8,  0.2,
                       0.3,  0.4,
                       0.1, 0.5),
                    nrow = Kw, byrow = TRUE)

Sigma_delta_true <- matrix(c(0.15, 0.02,
                             0.02, 0.10), 2, 2)

true_tauij0 <- Pred_tau(x10, w10, w20, Gammas = Psi_true, Sigmas = Sigma_delta_true) 

S <- 1000
taui0j040 <- sapply(1:S, function(s) {Pred_tau(x10, w10, w20, Gammas = draws40[["Psi"]][s,,], 
                                               Sigmas = draws40[["Sigma_delta"]][s,,])})

taui0j030 <- sapply(1:S, function(s) {Pred_tau(x10, w10, w20, Gammas = draws30[["Psi"]][s,,], 
                                               Sigmas = draws30[["Sigma_delta"]][s,,])})

taui0j020 <- sapply(1:S, function(s) {Pred_tau(x10, w10, w20, Gammas = draws20[["Psi"]][s,,], 
                                               Sigmas = draws20[["Sigma_delta"]][s,,])})

taui0j010 <- sapply(1:S, function(s) {Pred_tau(x10, w10, w20, Gammas = draws10[["Psi"]][s,,], 
                                            Sigmas = draws10[["Sigma_delta"]][s,,])})

# --- long df (panel-specific true values) ---
df <- bind_rows(
  tibble(value = taui0j040, J = "J = 40", true = true_tauij0),
  tibble(value = taui0j030, J = "J = 30", true = true_tauij0),
  tibble(value = taui0j020, J = "J = 20", true = true_tauij0),
  tibble(value = taui0j010, J = "J = 10", true = true_tauij0)
) %>%
  mutate(J = factor(J, levels = c("J = 40","J = 30","J = 20","J = 10")))

# 95% equal-tailed CI per panel
ci_df <- df %>%
  group_by(J) %>%
  summarise(
    ci_low  = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    true    = first(true),
    .groups = "drop"
  )

# consistent axes across panels
x_lim <- range(df$value, ci_df$true, ci_df$ci_low, ci_df$ci_high, na.rm = TRUE)

tmp_counts <- ggplot_build(
  ggplot(df, aes(x = value)) +
    geom_histogram(bins = 80, boundary = 0, closed = "left") +
    facet_wrap(~ J, ncol = 2)
)$data[[1]]$count
y_max <- max(tmp_counts, na.rm = TRUE)

# ---- plot ----
ggplot(df, aes(x = value)) +
  geom_histogram(
    bins = 80,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = true, linetype = "True value"),
    linewidth = 1.0,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_low, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_vline(
    data = ci_df,
    aes(xintercept = ci_high, linetype = "95% CI"),
    linewidth = 0.8,
    color = "black"
  ) +
  facet_wrap(~ J, ncol = 2) +
  scale_linetype_manual(
    name = NULL,
    values = c("True value" = "solid", "95% CI" = "dashed")
  ) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_max)) +
  labs(
    # title = bquote("Posterior draws of " * tau[.(i) * "," * .(j)]),
    x = bquote(tau[0]),
    y = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "plain")
  )
