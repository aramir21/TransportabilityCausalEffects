rm(list = ls())
set.seed(10101)
library(dplyr)
library(AER)
data("STAR")
STARorigin <- STAR

STAR <- STAR %>%
  mutate(
    schoolidk_raw = schoolidk,
    schoolidk = if_else(
      is.na(schoolidk_raw),
      NA_integer_,
      as.integer(factor(schoolidk_raw))
    )
  )

J <- nlevels(factor(STAR$schoolidk_raw))

N <- nrow(STAR)
STAR$id <- 1:N
STAR$smallk <- ifelse(STAR$stark == "small", 1, 0)
table(STAR$smallk)
yV1 <- data.frame(STAR$id, STAR$readk)
yV2 <- na.omit(yV1)

# Individual-level X
library(dplyr)
library(stringr)

STAR <- STAR %>%
  mutate(
    birth = as.character(birth),
    
    ## Year and age
    birth_year = as.integer(str_extract(birth, "^\\d{4}")),
    age = 1985 - birth_year,
    
    ## Quarter dummies
    Q1 = as.integer(str_detect(birth, "Q1")),
    Q2 = as.integer(str_detect(birth, "Q2")),
    Q3 = as.integer(str_detect(birth, "Q3")),
    Q4 = as.integer(str_detect(birth, "Q4")),
    
    eth_cauc = as.integer(ethnicity == "cauc"),
    eth_afam = as.integer(ethnicity == "afam"),
    male     = as.integer(gender == "male"),
    female   = as.integer(gender == "female")
  )
attach(STAR)

XV1  <- data.frame(id, 1, female, eth_cauc, eth_afam, age, Q1, Q2, Q3)
XV2 <- na.omit(XV1)

# Z for heterogeneous treatment effect: intercept + x1 only
ZV1  <- cbind(id, 1, female, eth_cauc, eth_afam)   
ZV2 <- na.omit(ZV1) 

# Study-level covariates W_j   (Kw = 3 here)
library(tidyr)
library(forcats)
library(rlang)

STAR2 <- STAR %>%
  mutate(
    tethnicityk = fct_na_value_to_level(tethnicityk, level = "NA"),
    ladderk     = fct_na_value_to_level(ladderk,     level = "NA"),
    degreek     = fct_na_value_to_level(degreek,     level = "NA"),
    schoolk     = fct_na_value_to_level(schoolk,     level = "NA")
  )

prop_wide <- function(data, var, prefix){
  var_sym <- rlang::sym(var)
  lvls <- levels(data[[var]])
  
  data %>%
    count(schoolidk, !!var_sym, name = "n") %>%
    group_by(schoolidk) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(schoolidk, !!var_sym := lvls, fill = list(n = 0, prop = 0)) %>%
    select(schoolidk, !!var_sym, prop) %>%   # <<< critical
    pivot_wider(
      id_cols = schoolidk,                   # <<< critical
      names_from  = !!var_sym,
      values_from = prop,
      names_prefix = prefix,
      values_fill  = 0
    )
}

teacher_stats_wide <- STAR2 %>%
  group_by(schoolidk) %>%
  summarise(avg_experience = mean(experiencek, na.rm = TRUE), .groups = "drop") %>%
  left_join(prop_wide(STAR2, "tethnicityk", "prop_eth_"),    by = "schoolidk") %>%
  left_join(prop_wide(STAR2, "ladderk",     "prop_ladder_"), by = "schoolidk") %>%
  left_join(prop_wide(STAR2, "degreek",     "prop_degree_"), by = "schoolidk") %>% 
  left_join(prop_wide(STAR2, "schoolk",     "prop_schoolk_"),by = "schoolidk")

summary(teacher_stats_wide)

W_study <- na.omit(teacher_stats_wide)
W_study <- W_study[,-c(3, 6, 13, 18)]


# Expand W to individual level (so we can store it in Data)
STAR <- STAR %>%
  mutate(schoolidk = as.character(schoolidk))

W_study <- W_study %>%
  mutate(schoolidk = as.character(schoolidk))

STAR_merged <- STAR %>%
  left_join(W_study, by = "schoolidk")


W_fullV1 <- STAR_merged[, c(49, 44, 61:77)]   
W_fullV2 <- na.omit(W_fullV1)
summary(W_fullV2)
W_fullV2 <- W_fullV2[, -16]

# Treatment assignment
TrV1 <- data.frame(STAR$id, STAR$smallk)
TrV2 <- na.omit(TrV1)

# common ids
ids_common <- Reduce(intersect, list(
  W_fullV2$id,
  XV2$id,
  yV2$STAR.id,
  TrV2$STAR.id
))

NJ <- length(ids_common)

# subset each dataset
W <- W_fullV2 %>%
  filter(id %in% ids_common) %>%
  arrange(id) %>%
  select(
    avg_experience,
    prop_eth_cauc, prop_eth_NA,
    prop_ladder_level1, prop_ladder_level2, prop_ladder_level3,
    prop_ladder_NA, prop_ladder_pending, prop_ladder_probation,
    prop_degree_master, `prop_degree_master+`,
    prop_degree_NA, prop_degree_specialist,
    prop_schoolk_rural, prop_schoolk_suburban, prop_schoolk_urban
  )

# Proportion of other races, mainly afroamerican, and asian, hispanic, amindian,  other,
# are the baseline
# Teachers's ladder has as baseline apprentice
# Teachers's degree has as baseline bachelor
# School's location has as baseline inner-city 

W <- cbind(1, as.matrix(W))

X <- XV2 %>%
  filter(id %in% ids_common) %>%
  arrange(id) %>%
  select(female, eth_cauc, eth_afam, age, Q1, Q2, Q3) 

# Ethnicity has as baseline other race (asian, hispanic, amindian,  other)
# Q4 is the baseline for quarters of birth
# Male is the baseline for female
X <- cbind(1, as.matrix(X))

Z <- X[,1:4]

H <- W[,c(1:3, 15:17)]

W <- W[,-1]

y  <- yV2  %>% filter(STAR.id %in% ids_common) %>% arrange(STAR.id) %>% pull(STAR.readk)
Tr <- TrV2 %>% filter(STAR.id %in% ids_common) %>% arrange(STAR.id) %>% pull(STAR.smallk)
id <- W_fullV2 %>% filter(id %in% ids_common) %>% arrange(id) %>% pull(schoolidk)

# Put everything in a data.frame including w's
Data <- data.frame(
  id = id,
  y  = y,
  Tr = Tr,
  x1 = X[,1], x2 = X[,2], x3 = X[,3], x4 = X[,4], x5 = X[,5], 
  x6 = X[,6], x7 = X[,7], x8 = X[,8],
  w1 = W[, 1], w2 = W[, 2], w3 = W[, 3], w4 = W[, 4], w5 = W[, 5],
  w6 = W[, 6], w7 = W[, 7], w8 = W[, 8], w9 = W[, 9], w10 = W[, 10],
  w11 = W[, 11], w12 = W[, 12], w13 = W[, 13], w14 = W[, 14], w15 = W[, 15],
  w16 = W[, 16],
  z1 = X[,1], z2 = X[,2], z3 = X[,3], z4 = X[,4],
  h1 = W[,1], h2 = W[,2], h3 = W[,3],
  h4 = W[,15], h5 = W[,16]
)
Kx <- ncol(X) 
Kw <- ncol(W)
Kz <- ncol(Z)
Kh <- ncol(H)

RegOLS <- lm(y ~ Tr + x2 + x3 + x4 + x5 + x6 + x7 + x8 + w1 +
               w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 +
               w11 + w12 + w13 + w14 + w15 + w16, data = Data)
summary(RegOLS)

stan_data <- list(
  N  = NJ,
  J  = J,
  Kx = Kx,
  Kz = Kz,
  Kw = Kw,
  Kh = Kh, 
  id = as.numeric(id),
  y  = as.vector(y),
  Tr = Tr,
  X  = X,
  Z  = Z,
  W  = W,
  H = H
)

stan_hier_meta_W_gamma_hetero <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> Kx;
  int<lower=1> Kz;
  int<lower=1> Kw;
  int<lower=1> Kh;

  int<lower=1,upper=J> id[N];
  vector[N] y;
  int<lower=0,upper=1> Tr[N];

  matrix[N, Kx] X;
  matrix[N, Kz] Z;

  matrix[N, Kw] W;
  matrix[N, Kh] H;
}

transformed data {
  int first_idx[J];
  matrix[J, Kh] H_study;

  for (j in 1:J) first_idx[j] = 0;

  for (n in 1:N) {
    int j = id[n];
    if (first_idx[j] == 0) first_idx[j] = n;
  }

  for (j in 1:J) {
    H_study[j] = H[first_idx[j]];
  }
}

parameters {
  vector[Kx] beta;
  vector[Kw] gamma;

  matrix[Kh, Kz] Gamma;

  vector<lower=0>[Kz] tau;
  cholesky_factor_corr[Kz] L_Omega;
  matrix[J, Kz] delta_raw;

  // sigma^2 and mixture weights
  real<lower=0> sigma2;
  vector<lower=0>[N] lambda;
}

transformed parameters {
  matrix[J, Kz] delta;
  matrix[J, Kz] delta_mean;
  matrix[Kz, Kz] L_Sigma;

  L_Sigma    = diag_pre_multiply(tau, L_Omega);
  delta_mean = H_study * Gamma;

  for (j in 1:J) {
    delta[j] = ( (delta_mean[j])' + L_Sigma * (delta_raw[j])' )';
  }
}

model {
  // constants
  real a = 0.01;
  real b = 0.01;
  real nu = 5;

  // Priors
  beta               ~ normal(0, 2);
  gamma              ~ normal(0, 2);
  to_vector(Gamma)   ~ normal(0, 1);

  tau                ~ cauchy(0, 2.5);
  L_Omega            ~ lkj_corr_cholesky(2.0);
  to_vector(delta_raw) ~ normal(0, 1);

  // sigma^2 ~ IG(a/2, b/2)
  sigma2 ~ inv_gamma(a / 2, b / 2);

  // lambda_i ~ Gamma(nu/2, nu/2)
  lambda ~ gamma(nu / 2, nu / 2);

  // Likelihood
  for (n in 1:N) {
    int j = id[n];
    real theta_n = dot_product(Z[n], delta[j]);
    real mu_n    = dot_product(X[n], beta)
                 + dot_product(W[n], gamma)
                 + Tr[n] * theta_n;

    y[n] ~ normal(mu_n, sqrt(sigma2 / lambda[n]));
  }
}

generated quantities {
  cov_matrix[Kz] Sigma_delta;
  real<lower=0> sigma;

  Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma);
  sigma = sqrt(sigma2);
}
"

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(
  model_code = stan_hier_meta_W_gamma_hetero,
  data   = stan_data,
  chains = 1,
  iter   = 2000,
  warmup = 1000,
  seed   = 1234
)

save(fit, file = "AppSTAR_Hetero.RData")

# load(file = "AppSTAR.RData")
# post <- summary(fit, pars = c("beta", "gamma", "Gamma", "Sigma_delta", "sigma"),
#                 probs = c(0.025, 0.5, 0.975))
# post

post<- rstan::extract(fit, pars = c("beta", "gamma", "Gamma", "Sigma_delta", "sigma"))
summary(coda::mcmc(post[["beta"]]))
summary(coda::mcmc(post[["gamma"]]))
summary(coda::mcmc(post[["Gamma"]][,,1]))
summary(coda::mcmc(post[["Gamma"]][,,2]))
summary(coda::mcmc(post[["Gamma"]][,,3]))
summary(coda::mcmc(post[["Gamma"]][,,4]))
summary(coda::mcmc(post[["Sigma_delta"]][,,1]))
summary(coda::mcmc(post[["Sigma_delta"]][,,2]))
summary(coda::mcmc(post[["Sigma_delta"]][,,3]))
summary(coda::mcmc(post[["Sigma_delta"]][,,4]))
summary(coda::mcmc(matrix(post[["sigma"]], fit@sim[["warmup2"]], 1)))

rstan::stan_trace(fit, pars = "beta")
rstan::stan_trace(fit, pars = "gamma")
rstan::stan_trace(fit, pars = "Gamma")
rstan::stan_trace(fit, pars = "Sigma_delta")
rstan::stan_trace(fit, pars = "sigma")

##### Trace plots #####
library(ggplot2)
library(dplyr)
library(tidyr)

mcmc_to_df <- function(x, name = "param") {
  if (is.null(dim(x))) {
    df <- data.frame(iter = seq_along(x), value = x, param = name)
  } else {
    x <- as.matrix(x)
    df <- as.data.frame(x)
    df$iter <- seq_len(nrow(df))
    df <- pivot_longer(
      df,
      cols = -iter,
      names_to = "param",
      values_to = "value"
    )
  }
  df
}

traceplot_gg <- function(df, title = NULL) {
  ggplot(df, aes(x = iter, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.3) +
    facet_wrap(~ param, scales = "free_y") +
    labs(
      title = title,
      x = "Iteration",
      y = "Value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

df_beta <- mcmc_to_df(post[["beta"]])
p_beta <- traceplot_gg(df_beta, "Trace plots: beta")
ggsave("appendix_trace_beta.pdf", p_beta, width = 8, height = 10)

df_gamma <- mcmc_to_df(post[["gamma"]])
p_gamma <- traceplot_gg(df_gamma, "Trace plots: gamma")
ggsave("appendix_trace_gamma.pdf", p_gamma, width = 8, height = 10)

for (k in 1:4) {
  df_Gamma <- mcmc_to_df(post[["Gamma"]][,,k], name = paste0("Gamma_", k))
  p <- traceplot_gg(df_Gamma, paste("Trace plots: Gamma[,,", k, "]"))
  ggsave(paste0("appendix_trace_Gamma_", k, ".pdf"),
         p, width = 8, height = 10)
}

for (k in 1:4) {
  df_Sd <- mcmc_to_df(post[["Sigma_delta"]][,,k], name = paste0("Sigma_delta_", k))
  p <- traceplot_gg(df_Sd, paste("Trace plots: Sigma[,,", k, "]"))
  ggsave(paste0("appendix_trace_Sigma_delta_", k, ".pdf"),
         p, width = 8, height = 10)
}

df_sigma <- mcmc_to_df(post[["sigma"]], "sigma")
p_sigma <- traceplot_gg(df_sigma, "Trace plot: sigma")
ggsave("appendix_trace_sigma.pdf", p_sigma, width = 6, height = 4)

##### Nice plots of tau (treatment effects) #####
library(dplyr)
library(purrr)
library(ggplot2)

post_delta <- rstan::extract(fit, pars = "delta")$delta  # iter x J x Kz

# choose 10 students (change these indices)
Ni <- 10
idx <- sample(1:NJ, Ni, replace = FALSE)

# posterior draws of tau_ij (theta) for each selected student
draws_df <- map_dfr(idx, function(i){
  j_i <- stan_data$id[i]
  z_i <- stan_data$Z[i, ]                          # length Kz
  tau_draws <- as.vector(post_delta[, j_i, ] %*% z_i)
  
  tibble(
    student = factor(paste0("i = ", i), levels = paste0("i = ", idx)),
    tau = tau_draws
  )
})

# mean + 95% equal-tail CI
sum_df <- draws_df %>%
  group_by(student) %>%
  summarise(
    mean = mean(tau),
    lo   = quantile(tau, 0.025),
    hi   = quantile(tau, 0.975),
    .groups = "drop"
  )

# histogram + vlines, 5x2 panels
ggplot(draws_df, aes(x = tau)) +
  geom_histogram(
    bins = 60,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(data = sum_df, aes(xintercept = mean, linetype = "Mean"),
             linewidth = 1.0, color = "black") +
  geom_vline(data = sum_df, aes(xintercept = lo, linetype = "95% CI"),
             linewidth = 0.8, color = "black") +
  geom_vline(data = sum_df, aes(xintercept = hi, linetype = "95% CI"),
             linewidth = 0.8, color = "black") +
  facet_wrap(~ student, ncol = 2, scales = "free_y") +
  scale_linetype_manual(values = c("95% CI" = "dashed", "Mean" = "solid")) +
  labs(x = expression(tau[ij]), y = "Count", linetype = "") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  ) + 
  scale_x_continuous(
    breaks = seq(-60, 60, by = 10),
    minor_breaks = seq(-2, 4, by = 0.25)
  )

##### Table #####
library(dplyr)
library(purrr)
library(tibble)
library(knitr)
library(kableExtra)

# safer default names if Z/H have empty colnames
z_names <- colnames(stan_data$Z)
if (is.null(z_names) || any(is.na(z_names)) || any(z_names == "")) z_names <- paste0("z", 1:Kz)

h_names <- colnames(stan_data$H)
if (is.null(h_names) || any(is.na(h_names)) || any(h_names == "")) h_names <- paste0("h", 1:Kh)

# drop z1 and h1 (and keep remaining order)
z_keep <- z_names[-1]
h_keep <- h_names[-1]

tab <- map_dfr(idx, function(i){
  
  j_i <- stan_data$id[i]
  z_i <- as.numeric(stan_data$Z[i, ])
  h_j <- round(as.numeric(stan_data$H[i, ]), 2)
  
  tau_draws <- as.vector(post_delta[, j_i, ] %*% z_i)
  
  tibble(
    i = i,
    j = j_i,
    mean = mean(tau_draws),
    lo   = quantile(tau_draws, 0.025),
    hi   = quantile(tau_draws, 0.975)
  ) %>%
    bind_cols(as_tibble(t(z_i), .name_repair = "minimal") %>% setNames(z_names)) %>%
    bind_cols(as_tibble(t(h_j), .name_repair = "minimal") %>% setNames(h_names))
})

tab <- tab %>%
  select(i, j, all_of(z_keep), all_of(h_keep), mean, lo, hi) %>%   # keep order
  mutate(
    across(
      where(is.numeric) & !c(i, j),
      ~ round(.x, 2)
    )
  )
latex_tab <- kable(
  tab,
  format = "latex",
  booktabs = TRUE,
  caption = "Posterior summaries for $\\tau_{ij}$ with regressors $z_i$ and study features $h_{j_i}$.",
  align = "r"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

save_kable(latex_tab, file = "tau_ij_table.tex")

# posterior draws of tau_ij (theta) for each selected student
idx <- 1:NJ
draws_dfAll <- sapply(1:NJ, function(i){
  j_i <- stan_data$id[i]
  z_i <- stan_data$Z[i, ]                          
  tau_draws <- as.vector(post_delta[, j_i, ] %*% z_i)
  ResPost_tau <- c(mean(tau_draws), quantile(tau_draws, c(0.025, 0.975)))}
  )

# tau_mat: 3 x N
tau_mean <- as.numeric(draws_dfAll[1, ])
tau_lo   <- as.numeric(draws_dfAll[2, ])
tau_hi   <- as.numeric(draws_dfAll[3, ])

# % of 95% CrIs that exclude 0
exclude0 <- (tau_lo > 0) | (tau_hi < 0)
pct_exclude0 <- 100 * mean(exclude0, na.rm = TRUE)
pct_exclude0

df_tau <- data.frame(tau_mean = tau_mean)

ggplot(df_tau, aes(x = tau_mean)) +
  geom_histogram(
    bins = 40,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(xintercept = 0, linewidth = 1.0) +
  labs(
    x = expression(paste("Posterior mean of ", tau[ij])),
    y = "Count",
    caption = sprintf("Share of 95%% credible intervals excluding 0: %.1f%%", pct_exclude0)
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

c(N = length(tau_mean), n_excluding0 = sum(exclude0, na.rm = TRUE), pct_excluding0 = pct_exclude0)

##### Prediction treatment effects #####
# Predictive function
Pred_tau <- function(x10, w10, Gammas, Sigmas){
  # Gammas <- post[["Gamma"]][1,,]
  # Sigmas <- post[["Sigma_delta"]][1,,]
  zij0 <- c(1, x10)
  hj0 <- c(1, w10)
  mean_theta <- t(Gammas)%*%hj0
  thetaj0 <- MASS::mvrnorm(1, mu = mean_theta, Sigmas)
  tauj0 <-  zij0%*%thetaj0
  return(tauj0)
}

S <- 1000

x10 <- c(0, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(7.95, 1, 0, 0, 0) # (Av. Exp, Prop. Cau., Rural, Urban, SubUrban)
tauij1 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                               Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(1, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(7.95, 1, 0, 0, 0) # (Av. Exp, Prop. Cau., Rural, Urban, SubUrban)
tauij2 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                            Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(0, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(4, 1, 0, 0, 0) # (Av. Exp, Prop. Cau., Rural, Urban, SubUrban)
tauij3 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                            Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(1, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(4, 1, 0, 0, 0) # (Av. Exp, Prop. Cau., Rural, Urban, SubUrban)
tauij4 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                            Sigmas = post[["Sigma_delta"]][s,,])})

### Conditional distribution of the treatment effect of a student in the sample
i <- 2623
j <- stan_data$id[i]
z_i <- stan_data$Z[i, ]
post_delta <- rstan::extract(fit, pars="delta")$delta  # iter x J x Kz
tau_id2623 <- as.vector(post_delta[, j, ] %*% z_i)

### Histograms Predictive treatment effect ###

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Put all posterior draws in one long data frame
df_tau <- tibble(
  tau = c(tauij1, tauij2, tauij3, tauij4),
  scenario = factor(rep(
    c("Male, Exp = 8",
      "Female,   Exp = 8",
      "Male, Exp = 4",
      "Female,   Exp = 4"),
    each = S
  ), levels = c("Male, Exp = 8","Female,   Exp = 8","Male, Exp = 4","Female,   Exp = 4"))
)

# Posterior mean + central (equal-tailed) 95% credible interval
sum_tau <- df_tau %>%
  group_by(scenario) %>%
  summarise(
    mean = mean(tau),
    lo   = quantile(tau, 0.025),
    hi   = quantile(tau, 0.975),
    .groups = "drop"
  )
sum_tau_long <- sum_tau %>%
  pivot_longer(cols = c(mean, lo, hi),
               names_to = "type", values_to = "x") %>%
  mutate(type = factor(type,
                       levels = c("lo","mean","hi"),
                       labels = c("95% CI","Mean","95% CI")))

# --- conditional draws for student 2623, only in the top-left facet ---
scenario_target <- levels(df_tau$scenario)[1]  # top-left facet: "Male, Exp = 8"

df_2623 <- tibble(
  tau = tau_id2623,
  scenario = factor(scenario_target, levels = levels(df_tau$scenario))
)

# conditional mean line (dotted) with legend entry
cond_mean_df <- tibble(
  scenario = factor(scenario_target, levels = levels(df_tau$scenario)),
  type = factor("Conditional mean (student 2623)",
                levels = c("95% CI", "Mean", "Conditional mean (student 2623)")),
  x = mean(tau_id2623)
)

# --- make the predictive line types also include the new legend level ---
sum_tau_long2 <- sum_tau_long %>%
  mutate(type = factor(as.character(type),
                       levels = c("95% CI", "Mean", "Conditional mean (student 2623)")))

ggplot(df_tau, aes(x = tau)) +
  # Predictive histogram (as before)
  geom_histogram(aes(y = after_stat(count)),
                 bins = 35, fill = "grey80", color = "grey80") +
  
  # Overlay: conditional histogram for student 2623 (ONLY in top-left facet)
  geom_histogram(data = df_2623,
                 aes(y = after_stat(count)),
                 bins = 35, fill = "grey50", color = "grey50",
                 alpha = 0.35) +
  
  # Predictive CI + mean (legend controlled by linetype)
  geom_vline(data = sum_tau_long2,
             aes(xintercept = x, linetype = type),
             linewidth = 1.0, color = "black") +
  
  # Conditional mean (dotted) only in top-left facet + legend entry
  geom_vline(data = cond_mean_df,
             aes(xintercept = x, linetype = type),
             linewidth = 1.1, color = "black") +
  
  scale_linetype_manual(
    values = c("95% CI" = "dashed",
               "Mean" = "solid",
               "Conditional mean (student 2623)" = "dotted")
  ) +
  
  facet_wrap(~ scenario, ncol = 2, scales = "free_y") +
  labs(x = expression(tau[0]),
       y = "Count",
       linetype = "") +
  
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal"
  )
