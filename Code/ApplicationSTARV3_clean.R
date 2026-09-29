# Hierarchical treatment-effect analysis of Project STAR
#
# This script prepares a single grade/subject sample, fits a hierarchical
# Student-t model in Stan, and reports school-level treatment effects and an
# empirical external-validity decomposition. Outcomes are expressed in
# control-group standard-deviation units in the final estimation sample.

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

config <- list(
  grade = "kindergarten",       # kindergarten, first, second, or third
  subject = "reading",         # reading or math
  seed = 10101L,
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  adapt_delta = 0.95,
  selected_school = 1L,
  fit_model = TRUE,
  save_plots = TRUE,
  # Used by ApplicationSTARV3_LOSO.R to load preprocessing and stan_code
  # without estimating the full-sample model.
  prepare_only = isTRUE(getOption("star.prepare_only", FALSE)),
  output_dir = "star_output"
)

required_packages <- c(
  "AER", "dplyr", "forcats", "ggplot2", "MASS", "rlang", "rstan",
  "stringr", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the following packages: ", paste(missing_packages, collapse = ", "))
}

if (!config$grade %in% c("kindergarten", "first", "second", "third")) {
  stop("config$grade is invalid.")
}
if (!config$subject %in% c("reading", "math")) {
  stop("config$subject is invalid.")
}
if (config$warmup >= config$iter) {
  stop("config$warmup must be smaller than config$iter.")
}

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(config$seed)

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

data("STAR", package = "AER")
STAR$id <- seq_len(nrow(STAR))

grade_suffix <- c(
  kindergarten = "k",
  first = "1",
  second = "2",
  third = "3"
)
assessment_year <- c(
  kindergarten = 1985L,
  first = 1986L,
  second = 1987L,
  third = 1988L
)

suffix <- unname(grade_suffix[config$grade])
source_names <- c(
  schoolid = paste0("schoolid", suffix),
  class_type = paste0("star", suffix),
  outcome = paste0(if (config$subject == "reading") "read" else "math", suffix),
  teacher_experience = paste0("experience", suffix),
  teacher_ethnicity = paste0("tethnicity", suffix),
  teacher_ladder = paste0("ladder", suffix),
  teacher_degree = paste0("degree", suffix),
  school_location = paste0("school", suffix),
  lunch = paste0("lunch", suffix)
)

required_columns <- c("id", "birth", "ethnicity", "gender", source_names)
missing_columns <- setdiff(required_columns, names(STAR))
if (length(missing_columns) > 0L) {
  stop("STAR is missing: ", paste(missing_columns, collapse = ", "))
}

star_selected <- STAR |>
  dplyr::transmute(
    id = .data$id,
    schoolid_raw = .data[[source_names[["schoolid"]]]],
    class_type = as.character(.data[[source_names[["class_type"]]]]),
    y_raw = as.numeric(.data[[source_names[["outcome"]]]]),
    teacher_experience = as.numeric(
      .data[[source_names[["teacher_experience"]]]]
    ),
    teacher_ethnicity = forcats::fct_na_value_to_level(
      factor(.data[[source_names[["teacher_ethnicity"]]]]), "missing"
    ),
    teacher_ladder = forcats::fct_na_value_to_level(
      factor(.data[[source_names[["teacher_ladder"]]]]), "missing"
    ),
    teacher_degree = forcats::fct_na_value_to_level(
      factor(.data[[source_names[["teacher_degree"]]]]), "missing"
    ),
    school_location = forcats::fct_na_value_to_level(
      factor(.data[[source_names[["school_location"]]]]), "missing"
    ),
    lunch_raw = as.character(.data[[source_names[["lunch"]]]]),
    birth = as.character(.data$birth),
    ethnicity = as.character(.data$ethnicity),
    gender = as.character(.data$gender)
  ) |>
  dplyr::mutate(
    # The treatment contrast is small class versus all other observed classes.
    Tr = dplyr::case_when(
      is.na(.data$class_type) ~ NA_integer_,
      .data$class_type == "small" ~ 1L,
      TRUE ~ 0L
    ),
    free_lunch = dplyr::case_when(
      is.na(.data$lunch_raw) ~ NA_integer_,
      .data$lunch_raw == "free" ~ 1L,
      .data$lunch_raw == "non-free" ~ 0L,
      TRUE ~ NA_integer_
    ),
    birth_year = as.integer(stringr::str_extract(.data$birth, "^\\d{4}")),
    age = unname(assessment_year[config$grade]) - .data$birth_year,
    Q1 = as.integer(stringr::str_detect(.data$birth, "Q1")),
    Q2 = as.integer(stringr::str_detect(.data$birth, "Q2")),
    Q3 = as.integer(stringr::str_detect(.data$birth, "Q3")),
    female = as.integer(.data$gender == "female"),
    eth_cauc = as.integer(.data$ethnicity == "cauc"),
    eth_afam = as.integer(.data$ethnicity == "afam")
  )

# Reindex nonmissing schools once so joins and Stan indices are stable.
school_levels <- sort(unique(star_selected$schoolid_raw[!is.na(star_selected$schoolid_raw)]))
star_selected$schoolid <- match(star_selected$schoolid_raw, school_levels)

# Proportions are student weighted because AER::STAR contains no teacher ID.
# One category from each factor is omitted below to avoid exact collinearity.
school_proportions <- function(data, variable, prefix) {
  variable_symbol <- rlang::sym(variable)
  factor_levels <- levels(data[[variable]])

  data |>
    dplyr::filter(!is.na(.data$schoolid)) |>
    dplyr::count(.data$schoolid, !!variable_symbol, name = "n") |>
    dplyr::group_by(.data$schoolid) |>
    dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    tidyr::complete(
      schoolid,
      !!variable_symbol := factor_levels,
      fill = list(n = 0L, proportion = 0)
    ) |>
    dplyr::select("schoolid", dplyr::all_of(variable), "proportion") |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(variable),
      values_from = "proportion",
      names_prefix = prefix,
      values_fill = 0
    )
}

school_data <- star_selected |>
  dplyr::filter(!is.na(.data$schoolid)) |>
  dplyr::group_by(.data$schoolid) |>
  dplyr::summarise(
    avg_experience = mean(.data$teacher_experience, na.rm = TRUE),
    prop_lunch_free = mean(.data$free_lunch, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    avg_experience = dplyr::na_if(.data$avg_experience, NaN),
    prop_lunch_free = dplyr::na_if(.data$prop_lunch_free, NaN)
  ) |>
  dplyr::left_join(
    school_proportions(star_selected, "teacher_ethnicity", "prop_teacher_eth_"),
    by = "schoolid"
  ) |>
  dplyr::left_join(
    school_proportions(star_selected, "teacher_ladder", "prop_ladder_"),
    by = "schoolid"
  ) |>
  dplyr::left_join(
    school_proportions(star_selected, "teacher_degree", "prop_degree_"),
    by = "schoolid"
  ) |>
  dplyr::left_join(
    school_proportions(star_selected, "school_location", "prop_school_"),
    by = "schoolid"
  )

w_requested <- c(
  "avg_experience",
  "prop_teacher_eth_cauc", "prop_teacher_eth_missing",
  "prop_ladder_level1", "prop_ladder_level2", "prop_ladder_level3",
  "prop_ladder_missing", "prop_ladder_pending", "prop_ladder_probation",
  "prop_degree_master", "prop_degree_master+", "prop_degree_missing",
  "prop_degree_specialist",
  "prop_school_rural", "prop_school_suburban", "prop_school_urban",
  "prop_lunch_free"
)
h_requested <- c("prop_teacher_eth_cauc", "prop_degree_master+")

missing_h <- setdiff(h_requested, names(school_data))
if (length(missing_h) > 0L) {
  stop("Required treatment-effect moderators are absent: ", paste(missing_h, collapse = ", "))
}

w_available <- intersect(w_requested, names(school_data))
if (length(setdiff(w_requested, w_available)) > 0L) {
  message(
    "Omitting categories absent in this sample: ",
    paste(setdiff(w_requested, w_available), collapse = ", ")
  )
}

analysis_data <- star_selected |>
  dplyr::select(
    "id", "schoolid", "schoolid_raw", "y_raw", "Tr",
    "female", "eth_cauc", "eth_afam", "age",
    "Q1", "Q2", "Q3", "free_lunch"
  ) |>
  dplyr::left_join(
    dplyr::select(
      school_data,
      "schoolid",
      dplyr::all_of(unique(c(w_available, h_requested)))
    ),
    by = "schoolid"
  ) |>
  tidyr::drop_na()

if (nrow(analysis_data) == 0L) {
  stop("No observations remain in the complete-case estimation sample.")
}

# Drop zero-variance W columns before scaling. Such coefficients are unidentified.
nonconstant <- function(x) is.numeric(x) && stats::sd(x) > 0
w_available <- w_available[vapply(analysis_data[w_available], nonconstant, logical(1))]
if (length(w_available) == 0L) {
  stop("No nonconstant school-level outcome covariates remain.")
}
if (!all(vapply(analysis_data[h_requested], nonconstant, logical(1)))) {
  stop("At least one required treatment-effect moderator is constant.")
}

# Reindex schools after all complete-case restrictions.
analysis_data <- analysis_data |>
  dplyr::arrange(.data$id) |>
  dplyr::mutate(stan_school = as.integer(factor(.data$schoolid)))

scale_vector <- function(x) {
  scaled <- as.numeric(scale(x))
  if (any(!is.finite(scaled))) stop("Cannot scale a constant or nonfinite variable.")
  scaled
}

X <- analysis_data |>
  dplyr::transmute(
    intercept = 1,
    female, eth_cauc, eth_afam,
    age = scale_vector(.data$age),
    Q1, Q2, Q3, free_lunch
  ) |>
  as.matrix()
Z <- X[, c("intercept", "free_lunch"), drop = FALSE]
W <- apply(as.matrix(analysis_data[w_available]), 2L, scale_vector)
W <- matrix(W, nrow = nrow(analysis_data), dimnames = list(NULL, w_available))
H <- cbind(
  intercept = 1,
  apply(as.matrix(analysis_data[h_requested]), 2L, scale_vector)
)
colnames(H) <- c("intercept", h_requested)

Tr <- analysis_data$Tr
y_raw <- analysis_data$y_raw
control_rows <- Tr == 0L
if (sum(control_rows) < 2L) stop("At least two controls are required.")
control_mean <- mean(y_raw[control_rows])
control_sd <- stats::sd(y_raw[control_rows])
if (!is.finite(control_sd) || control_sd <= 0) stop("Invalid control-group SD.")
y <- (y_raw - control_mean) / control_sd

id <- analysis_data$stan_school
N <- length(y)
J <- max(id)
stopifnot(
  identical(sort(unique(id)), seq_len(J)),
  nrow(X) == N, nrow(Z) == N, nrow(W) == N, nrow(H) == N,
  all(Tr %in% 0:1)
)

# H must be constant within school.
h_by_school <- split(seq_len(N), id)
stopifnot(all(vapply(
  h_by_school,
  function(rows) all(vapply(
    seq_len(ncol(H)),
    function(k) length(unique(H[rows, k])) == 1L,
    logical(1)
  )),
  logical(1)
)))

stan_data <- list(
  N = N, J = J,
  Kx = ncol(X), Kz = ncol(Z), Kw = ncol(W), Kh = ncol(H),
  id = id, y = y, Tr = Tr,
  X = X, Z = Z, W = W, H = H
)

write.csv(cbind(stan_data$id, stan_data$Tr, stan_data$H, stan_data$Z, stan_data$X[,c(1:2,5)], stan_data$W[,1]), file = "STAR_Data_Simulation.csv")

sample_description <- list(
  grade = config$grade,
  subject = config$subject,
  N = N,
  J = J,
  control_mean = control_mean,
  control_sd = control_sd,
  X_columns = colnames(X),
  Z_columns = colnames(Z),
  W_columns = colnames(W),
  H_columns = colnames(H)
)
print(sample_description)

# An unpooled OLS benchmark with unique names; treatment heterogeneity remains
# the object of the hierarchical model below.
ols_data <- data.frame(y = y, Tr = Tr, X[, -1, drop = FALSE], W, check.names = TRUE)
ols_fit <- stats::lm(y ~ ., data = ols_data)
print(summary(ols_fit))

# -----------------------------------------------------------------------------
# Hierarchical model
# -----------------------------------------------------------------------------

# For school j, delta_j = H_j Gamma + L_Sigma delta_raw_j and the individual
# treatment effect is Z_i delta_j. The Student-t likelihood is the marginalized
# form of the normal-gamma mixture, avoiding N weakly identified latent lambdas.
stan_code <- "
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
  matrix[N, Kh] H;
}
transformed data {
  int first_idx[J];
  matrix[J, Kh] H_school;
  for (j in 1:J) first_idx[j] = 0;
  for (n in 1:N) if (first_idx[id[n]] == 0) first_idx[id[n]] = n;
  for (j in 1:J) H_school[j] = H[first_idx[j]];
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
    y[n] ~ student_t(5, mu, sigma);
  }
}
generated quantities {
  corr_matrix[Kz] Omega = multiply_lower_tri_self_transpose(L_Omega);
  cov_matrix[Kz] Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma);
}
"

if (!isTRUE(config$prepare_only)) {

fit_path <- file.path(
  config$output_dir,
  paste0("STAR_", config$grade, "_", config$subject, "_fit.RData")
)

if (isTRUE(config$fit_model)) {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = max(1L, parallel::detectCores() - 1L))
  fit <- rstan::stan(
    model_code = stan_code,
    data = stan_data,
    chains = config$chains,
    iter = config$iter,
    warmup = config$warmup,
    seed = config$seed,
    control = list(adapt_delta = config$adapt_delta)
  )
  save(fit, sample_description, file = fit_path)
} else {
  if (!file.exists(fit_path)) stop("Saved fit not found: ", fit_path)
  load(fit_path) # loads fit and sample_description
}

print(
  fit,
  pars = c(
    "beta", "gamma", "Gamma", "sd_delta", "Omega", "Sigma_delta", "sigma"
  ),
  probs = c(0.025, 0.5, 0.975)
)
rstan::check_hmc_diagnostics(fit)

post <- rstan::extract(
  fit,
  pars = c("beta", "gamma", "Gamma", "delta", "Sigma_delta", "sigma")
)

summarise_draws <- function(draws) {
  c(
    mean = mean(draws),
    median = stats::median(draws),
    lower_95 = unname(stats::quantile(draws, 0.025)),
    upper_95 = unname(stats::quantile(draws, 0.975)),
    Pr_positive = mean(draws > 0),
    Pr_negative = mean(draws < 0)
  )
}

summarise_matrix_draws <- function(draws, names) {
  result <- as.data.frame(t(apply(draws, 2L, summarise_draws)))
  cbind(term = names, result, row.names = NULL)
}

beta_summary <- summarise_matrix_draws(post$beta, colnames(X))
gamma_summary <- summarise_matrix_draws(post$gamma, colnames(W))
Gamma_summaries <- lapply(seq_len(ncol(Z)), function(k) {
  output <- summarise_matrix_draws(post$Gamma[, , k], colnames(H))
  output$Z_moderator <- colnames(Z)[k]
  output[, c("Z_moderator", setdiff(names(output), "Z_moderator"))]
})
Gamma_summary <- dplyr::bind_rows(Gamma_summaries)

# Retain the unique elements of the symmetric treatment-effect covariance
# matrix. Diagonal entries are variances; off-diagonal entries are covariances.
Sigma_indices <- which(
  lower.tri(matrix(NA_real_, nrow = ncol(Z), ncol = ncol(Z)), diag = TRUE),
  arr.ind = TRUE
)
Sigma_delta_summary <- do.call(rbind, lapply(seq_len(nrow(Sigma_indices)), function(i) {
  row_index <- Sigma_indices[i, "row"]
  column_index <- Sigma_indices[i, "col"]
  draw_summary <- summarise_draws(
    post$Sigma_delta[, row_index, column_index]
  )

  data.frame(
    parameter = paste0(
      "Sigma_delta[", row_index, ",", column_index, "]"
    ),
    row_coefficient = colnames(Z)[row_index],
    column_coefficient = colnames(Z)[column_index],
    element = if (row_index == column_index) "variance" else "covariance",
    as.list(draw_summary),
    row.names = NULL,
    check.names = FALSE
  )
}))

utils::write.csv(beta_summary, file.path(config$output_dir, "beta_summary.csv"), row.names = FALSE)
utils::write.csv(gamma_summary, file.path(config$output_dir, "gamma_summary.csv"), row.names = FALSE)
utils::write.csv(Gamma_summary, file.path(config$output_dir, "Gamma_summary.csv"), row.names = FALSE)
utils::write.csv(
  Sigma_delta_summary,
  file.path(config$output_dir, "Sigma_delta_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# School treatment effects and new-school prediction
# -----------------------------------------------------------------------------

school_rows <- split(seq_len(N), factor(id, levels = seq_len(J)))
Zbar_school <- do.call(rbind, lapply(school_rows, function(rows) {
  colMeans(Z[rows, , drop = FALSE])
}))
H_school <- do.call(rbind, lapply(school_rows, function(rows) H[rows[1L], , drop = FALSE]))

# Draws x schools: average effect for the observed composition in each school.
tau_schools <- vapply(seq_len(J), function(j) {
  rowSums(post$delta[, j, , drop = FALSE][, 1L, ] *
            matrix(Zbar_school[j, ], nrow = nrow(post$delta), ncol = ncol(Z), byrow = TRUE))
}, numeric(dim(post$delta)[1L]))
colnames(tau_schools) <- paste0("school_", seq_len(J))

school_tau_summary <- as.data.frame(t(apply(tau_schools, 2L, summarise_draws)))
school_tau_summary <- cbind(school = seq_len(J), school_tau_summary, row.names = NULL)
utils::write.csv(
  school_tau_summary,
  file.path(config$output_dir, "school_treatment_effects.csv"),
  row.names = FALSE
)

j_selected <- config$selected_school
if (!j_selected %in% seq_len(J)) stop("config$selected_school is not observed.")
tau_observed <- tau_schools[, j_selected]

# Posterior prediction for a new school with the selected observed school's H
# and mean Z, including a new draw of its unexplained random coefficient.
set.seed(config$seed)
n_draws <- dim(post$Gamma)[1L]
tau_new <- vapply(seq_len(n_draws), function(s) {
  Gamma_s <- matrix(post$Gamma[s, , ], nrow = ncol(H), ncol = ncol(Z))
  Sigma_s <- matrix(post$Sigma_delta[s, , ], nrow = ncol(Z), ncol = ncol(Z))
  mean_delta <- drop(t(Gamma_s) %*% H_school[j_selected, ])
  delta_new <- MASS::mvrnorm(1L, mu = mean_delta, Sigma = Sigma_s)
  drop(crossprod(Zbar_school[j_selected, ], delta_new))
}, numeric(1))

selected_school_summary <- rbind(
  observed_school = summarise_draws(tau_observed),
  new_comparable_school = summarise_draws(tau_new)
)
utils::write.csv(
  selected_school_summary,
  file.path(config$output_dir, paste0("school_", j_selected, "_comparison.csv"))
)

# -----------------------------------------------------------------------------
# Empirical external-validity decomposition
# -----------------------------------------------------------------------------

# For each posterior draw, V_sys is the empirical variance across the observed
# schools of E[tau_j | H_j, Zbar_j]. V_res averages Zbar_j' Sigma Zbar_j.
external_validity <- t(vapply(seq_len(n_draws), function(s) {
  Gamma_s <- matrix(post$Gamma[s, , ], nrow = ncol(H), ncol = ncol(Z))
  Sigma_s <- matrix(post$Sigma_delta[s, , ], nrow = ncol(Z), ncol = ncol(Z))
  conditional_mean <- rowSums((H_school %*% Gamma_s) * Zbar_school)
  V_sys <- mean((conditional_mean - mean(conditional_mean))^2)
  V_res <- mean(rowSums((Zbar_school %*% Sigma_s) * Zbar_school))
  c(V_sys = V_sys, V_res = V_res, V_total = V_sys + V_res,
    R_EV2 = V_sys / (V_sys + V_res))
}, numeric(4)))
external_validity <- as.data.frame(external_validity)
external_validity_summary <- as.data.frame(t(vapply(
  external_validity,
  function(x) summarise_draws(x)[c("mean", "median", "lower_95", "upper_95")],
  numeric(4)
)))
utils::write.csv(
  external_validity_summary,
  file.path(config$output_dir, "external_validity_summary.csv")
)

# -----------------------------------------------------------------------------
# Essential figures
# -----------------------------------------------------------------------------

if (isTRUE(config$save_plots)) {
  comparison_data <- rbind(
    data.frame(tau = tau_observed, distribution = "Observed school"),
    data.frame(tau = tau_new, distribution = "New comparable school")
  )
  p_comparison <- ggplot2::ggplot(
    comparison_data,
    ggplot2::aes(x = .data$tau, fill = .data$distribution)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      position = "identity",
      bins = 40,
      alpha = 0.55,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(
      values = c(
        "Observed school" = "grey35",
        "New comparable school" = "grey75"
      )
    ) +
    ggplot2::labs(
      title = paste("Posterior ATE distributions for school", j_selected),
      subtitle = paste(
        "Observed school and out-of-sample school with the same H and mean Z"
      ),
      x = "Treatment effect (control-group SD units)", y = "Posterior density",
      fill = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(
    file.path(config$output_dir, paste0("school_", j_selected, "_comparison.pdf")),
    p_comparison, width = 8, height = 6
  )
  
  p_school_means <- ggplot2::ggplot(
    school_tau_summary,
    ggplot2::aes(x = .data$mean)
  ) +
    ggplot2::geom_histogram(
      bins = 25,
      fill = "grey65",
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "Posterior mean ATEs across observed schools",
      subtitle = paste(J, "schools in the estimation sample"),
      x = "Posterior mean school ATE (control-group SD units)",
      y = "Number of schools"
    ) +
    ggplot2::theme_classic()
  ggplot2::ggsave(
    file.path(config$output_dir, "observed_school_mean_ATEs.pdf"),
    p_school_means, width = 8, height = 6
  )
  
  p_external_validity <- ggplot2::ggplot(
    external_validity,
    ggplot2::aes(x = .data$R_EV2)
  ) +
    ggplot2::geom_histogram(bins = 40, fill = "grey50", color = "white") +
    ggplot2::labs(
      title = expression("Posterior distribution of " * R[EV]^2),
      subtitle = paste("Empirical population of", J, "observed schools"),
      x = expression(R[EV]^2), y = "Posterior draws"
    ) +
    ggplot2::theme_classic()
  ggplot2::ggsave(
    file.path(config$output_dir, "external_validity_R2.pdf"),
    p_external_validity, width = 8, height = 6
  )
}

message("Analysis complete. Outputs written to: ", normalizePath(config$output_dir))

}
