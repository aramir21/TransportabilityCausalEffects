rm(list = ls())
set.seed(10101)
library(dplyr)
library(tidyr)
library(stringr)
library(rlang)
library(AER)
data("STAR")
STARorigin <- STAR
N <- nrow(STAR)
STAR$id <- 1:N

select_star_sample <- function(
    data,
    grade = c("kindergarten", "first", "second", "third"),
    subject = c("reading", "math"),
    complete_cases = FALSE
) {
  grade <- match.arg(grade)
  subject <- match.arg(subject)

  # Create a stable observation identifier before filtering.
  if (!"id" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(id = dplyr::row_number())
  }

  # Map grade names to STAR variable suffixes.
  grade_suffix <- c(
    kindergarten = "k",
    first = "1",
    second = "2",
    third = "3"
  )

  suffix <- unname(grade_suffix[grade])

  # Construct grade-specific variable names.
  school_var <- paste0("schoolid", suffix)
  treatment_var <- paste0("star", suffix)

  outcome_var <- paste0(
    if (subject == "reading") "read" else "math",
    suffix
  )

  experience_var <- paste0("experience", suffix)
  teacher_ethnicity_var <- paste0("tethnicity", suffix)
  ladder_var <- paste0("ladder", suffix)
  degree_var <- paste0("degree", suffix)
  school_location_var <- paste0("school", suffix)
  lunch_var <- paste0("lunch", suffix)

  # Verify that all required variables exist.
  required_vars <- c(
    school_var,
    treatment_var,
    outcome_var,
    experience_var,
    teacher_ethnicity_var,
    ladder_var,
    degree_var,
    school_location_var,
    lunch_var
  )

  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0L) {
    stop(
      "The following variables are missing from data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Construct generic variables for the selected grade and subject.
  result <- data %>%
    dplyr::mutate(
      selected_grade = grade,
      selected_subject = subject,

      # School identifier.
      schoolid_raw = .data[[school_var]],

      schoolid = dplyr::if_else(
        is.na(schoolid_raw),
        NA_integer_,
        as.integer(factor(schoolid_raw))
      ),

      # Treatment assignment.
      class_type = as.character(.data[[treatment_var]]),

      Tr = dplyr::case_when(
        is.na(class_type) ~ NA_integer_,
        class_type == "small" ~ 1L,
        TRUE ~ 0L
      ),

      # Outcome.
      y_raw = as.numeric(.data[[outcome_var]]),

      # Teacher experience.
      teacher_experience =
        as.numeric(.data[[experience_var]]),

      # Teacher ethnicity.
      teacher_ethnicity =
        forcats::fct_na_value_to_level(
          factor(.data[[teacher_ethnicity_var]]),
          level = "NA"
        ),

      # Teacher career-ladder status.
      teacher_ladder =
        forcats::fct_na_value_to_level(
          factor(.data[[ladder_var]]),
          level = "NA"
        ),

      # Teacher degree.
      teacher_degree =
        forcats::fct_na_value_to_level(
          factor(.data[[degree_var]]),
          level = "NA"
        ),

      # School location.
      school_location =
        forcats::fct_na_value_to_level(
          factor(.data[[school_location_var]]),
          level = "NA"
        ),

      # Categorical lunch status, including an explicit missing category.
      lunch_status =
        forcats::fct_na_value_to_level(
          factor(.data[[lunch_var]]),
          level = "NA"
        ),

      # Binary free-lunch indicator.
      free_lunch = dplyr::case_when(
        is.na(.data[[lunch_var]]) ~ NA_integer_,
        as.character(.data[[lunch_var]]) == "free" ~ 1L,
        as.character(.data[[lunch_var]]) == "non-free" ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  # Optionally restrict to the basic complete-case sample.
  if (complete_cases) {
    result <- result %>%
      tidyr::drop_na(
        id,
        schoolid,
        Tr,
        y_raw
      )
  }

  # Store original STAR variable names.
  attr(result, "grade") <- grade
  attr(result, "subject") <- subject
  attr(result, "school_variable") <- school_var
  attr(result, "treatment_variable") <- treatment_var
  attr(result, "outcome_variable") <- outcome_var
  attr(result, "experience_variable") <- experience_var
  attr(result, "teacher_ethnicity_variable") <-
    teacher_ethnicity_var
  attr(result, "ladder_variable") <- ladder_var
  attr(result, "degree_variable") <- degree_var
  attr(result, "school_location_variable") <-
    school_location_var
  attr(result, "lunch_variable") <- lunch_var

  result
}

prepare_star_outcome <- function(
    data,
    standardize = c("control", "full_sample", "none")
) {
  standardize <- match.arg(standardize)

  # Check that select_star_sample() created the required variables
  required_vars <- c("id", "Tr", "y_raw")
  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from data: ",
      paste(missing_vars, collapse = ", "),
      ". Run select_star_sample() first."
    )
  }

  # Restrict the data to observations with an ID, treatment status,
  # outcome, and selected-grade school identifier.
  result <- data %>%
    tidyr::drop_na(id, schoolid, Tr, y_raw)

  # Calculate the reference mean and standard deviation
  if (standardize == "control") {
    reference_sample <- result$Tr == 0

    if (sum(reference_sample) < 2) {
      stop("At least two control observations are required.")
    }

    reference_mean <- mean(result$y_raw[reference_sample])
    reference_sd <- sd(result$y_raw[reference_sample])
    reference_group <- "Control group: Tr == 0"

  } else if (standardize == "full_sample") {
    reference_mean <- mean(result$y_raw)
    reference_sd <- sd(result$y_raw)
    reference_group <- "Full estimation sample"

  } else {
    reference_mean <- 0
    reference_sd <- 1
    reference_group <- "No standardization"
  }

  if (!is.finite(reference_sd) || reference_sd <= 0) {
    stop(
      "The reference-group standard deviation must be positive and finite."
    )
  }

  # Under standardize = "none", this simply gives y = y_raw.
  result <- result %>%
    mutate(
      y = (y_raw - reference_mean) / reference_sd
    )

  # Store the standardization constants for interpretation and replication
  attr(result, "standardization") <- standardize
  attr(result, "reference_group") <- reference_group
  attr(result, "reference_mean") <- reference_mean
  attr(result, "reference_sd") <- reference_sd

  result
}

# Step 1: Select grade and subject.
STAR_selected <- select_star_sample(
  data = STAR,
  grade = "kindergarten",
  subject = "reading", # reading or math
  complete_cases = FALSE
)

# Step 2: Construct the complete-case sample and standardize the outcome
STAR_outcome <- prepare_star_outcome(
  data = STAR_selected,
  standardize = "control"
)

# ---------------------------------------------------------------------------
# Individual-level covariates
# ---------------------------------------------------------------------------

# The STAR cohorts were observed in successive school years.  Consequently,
# age at assessment must change with the selected grade rather than always
# being calculated for the kindergarten year (1985).
assessment_year <- c(
  kindergarten = 1985L,
  first = 1986L,
  second = 1987L,
  third = 1988L
)

selected_grade <- unique(STAR_selected$selected_grade)

if (length(selected_grade) != 1L ||
    !selected_grade %in% names(assessment_year)) {
  stop("STAR_selected must contain exactly one valid selected_grade.")
}

individual_data <- STAR_selected %>%
  mutate(
    birth_character = as.character(birth),
    birth_year = as.integer(
      str_extract(birth_character, "^\\d{4}")
    ),
    age = unname(assessment_year[selected_grade]) - birth_year,
    Q1 = as.integer(str_detect(birth_character, "Q1")),
    Q2 = as.integer(str_detect(birth_character, "Q2")),
    Q3 = as.integer(str_detect(birth_character, "Q3")),
    Q4 = as.integer(str_detect(birth_character, "Q4")),
    eth_cauc = as.integer(ethnicity == "cauc"),
    eth_afam = as.integer(ethnicity == "afam"),
    female = as.integer(gender == "female")
  ) %>%
  select(
    id,
    female,
    eth_cauc,
    eth_afam,
    age,
    Q1,
    Q2,
    Q3,
    free_lunch
  ) %>%
  drop_na()

# Q4, male, ethnicities other than cauc/afam and no free lunch are the reference categories.

# ---------------------------------------------------------------------------
# School-level covariates for the selected grade
# ---------------------------------------------------------------------------

# Calculate the within-school proportion in every level of a categorical
# variable.  Because STAR has one row per student, teacher-characteristic
# proportions are student weighted unless the data are first reduced to one
# row per distinct teacher.
prop_wide <- function(data, var, prefix) {
  if (!var %in% names(data)) {
    stop("Variable not found: ", var)
  }

  if (!is.factor(data[[var]])) {
    stop("Variable must be a factor: ", var)
  }

  var_sym <- rlang::sym(var)
  lvls <- levels(data[[var]])

  data %>%
    filter(!is.na(schoolid)) %>%
    count(schoolid, !!var_sym, name = "n") %>%
    group_by(schoolid) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(
      schoolid,
      !!var_sym := lvls,
      fill = list(n = 0L, prop = 0)
    ) %>%
    select(schoolid, !!var_sym, prop) %>%
    pivot_wider(
      id_cols = schoolid,
      names_from = !!var_sym,
      values_from = prop,
      names_prefix = prefix,
      values_fill = 0
    )
}

teacher_stats_wide <- STAR_selected %>%
  filter(!is.na(schoolid)) %>%
  group_by(schoolid) %>%
  summarise(
    avg_experience = {
      value <- mean(teacher_experience, na.rm = TRUE)
      if (is.nan(value)) NA_real_ else value
    },
    .groups = "drop"
  ) %>%
  left_join(
    prop_wide(
      STAR_selected,
      "teacher_ethnicity",
      "prop_teacher_eth_"
    ),
    by = "schoolid"
  ) %>%
  left_join(
    prop_wide(
      STAR_selected,
      "teacher_ladder",
      "prop_ladder_"
    ),
    by = "schoolid"
  ) %>%
  left_join(
    prop_wide(
      STAR_selected,
      "teacher_degree",
      "prop_degree_"
    ),
    by = "schoolid"
  ) %>%
  left_join(
    prop_wide(
      STAR_selected,
      "school_location",
      "prop_school_"
    ),
    by = "schoolid"
  )

# Student free-lunch composition at the school level
student_stats_wide <- STAR_selected %>%
  filter(!is.na(schoolid)) %>%
  group_by(schoolid) %>%
  summarise(
    prop_lunch_free = {
      value <- mean(free_lunch, na.rm = TRUE)
      if (is.nan(value)) NA_real_ else value
    },
    .groups = "drop"
  )

# Combine teacher, school-location, and student-composition variables
school_stats_wide <- teacher_stats_wide %>%
  left_join(
    student_stats_wide,
    by = "schoolid"
  )

# Select covariates by name instead of relying on column positions. Categories
# absent from the selected grade are omitted because an all-zero regressor has
# no identified likelihood contribution. The reference categories are:
#   teacher ethnicity: non-caucasian;
#   ladder: apprentice;
#   degree: bachelor;
#   school location: inner city;
#   lunch status: non-free lunch.
w_covariates_requested <- c(
  "avg_experience",
  "prop_teacher_eth_cauc",
  "prop_teacher_eth_NA",
  "prop_ladder_level1",
  "prop_ladder_level2",
  "prop_ladder_level3",
  "prop_ladder_NA",
  "prop_ladder_pending",
  "prop_ladder_probation",
  "prop_degree_master",
  "prop_degree_master+",
  "prop_degree_NA",
  "prop_degree_specialist",
  "prop_school_rural",
  "prop_school_suburban",
  "prop_school_urban",
  "prop_lunch_free"
)

missing_w_covariates <- setdiff(
  w_covariates_requested,
  names(school_stats_wide)
)

if (length(missing_w_covariates) > 0L) {
  message(
    "Omitting categories absent from the selected grade: ",
    paste(missing_w_covariates, collapse = ", ")
  )
}

W_study <- school_stats_wide %>%
  select(
    schoolid,
    any_of(w_covariates_requested)
  ) %>%
  drop_na()

# Attach the selected-grade school characteristics to each student.
W_full <- STAR_selected %>%
  select(id, schoolid) %>%
  left_join(W_study, by = "schoolid") %>%
  drop_na()

# Selected outcome and treatment are already grade/subject generic.
outcome_data <- STAR_selected %>%
  select(id, y_raw, Tr) %>%
  drop_na()

# ---------------------------------------------------------------------------
# Common estimation sample and aligned model matrices
# ---------------------------------------------------------------------------

ids_common <- Reduce(
  intersect,
  list(
    W_full$id,
    individual_data$id,
    outcome_data$id
  )
)

ids_common <- sort(ids_common)
NJ <- length(ids_common)

if (NJ == 0L) {
  stop("No observations remain in the common estimation sample.")
}

aligned_individual <- individual_data %>%
  filter(id %in% ids_common) %>%
  arrange(id)

aligned_school <- W_full %>%
  filter(id %in% ids_common) %>%
  arrange(id)

aligned_outcome <- outcome_data %>%
  filter(id %in% ids_common) %>%
  arrange(id)

stopifnot(
  identical(aligned_individual$id, ids_common),
  identical(aligned_school$id, ids_common),
  identical(aligned_outcome$id, ids_common)
)

# X contains individual-level predictors; its first column is the intercept.
X <- aligned_individual %>%
  select(
    female,
    eth_cauc,
    eth_afam,
    age,
    Q1,
    Q2,
    Q3,
    free_lunch
  ) %>%
  mutate(
    age = as.numeric(scale(age))
  ) %>%
  as.matrix()

X <- cbind(intercept = 1, X)
X[, "age"] <- as.numeric(scale(X[, "age"]))

# Z determines treatment-effect heterogeneity: intercept and free lunch indicator.
Z <- X[, c("intercept", "free_lunch"), drop = FALSE]

# W contains school-level predictors of the outcome.  Its intercept is omitted
# because the common intercept is already included in X.
W <- aligned_school %>%
  select(-id, -schoolid) %>%
  as.matrix()

W <- scale(W)

# H contains school-level predictors of the heterogeneous treatment effects.
# The order is intercept, teacher caucasian proportion, and teacher master+ proportion.

h_covariates <- c(
  "prop_teacher_eth_cauc",
  "prop_degree_master+"
)

missing_h_covariates <- setdiff(h_covariates, names(aligned_school))

if (length(missing_h_covariates) > 0L) {
  stop(
    "The following required H covariates are absent from the selected grade: ",
    paste(missing_h_covariates, collapse = ", ")
  )
}

H <- aligned_school %>%
  select(all_of(h_covariates)) %>%
  as.matrix()
H <- cbind(intercept = 1, H)
H[, -1] <- scale(H[, -1, drop = FALSE])

id <- aligned_school$schoolid

stopifnot(
  length(id) == nrow(H)
)

H_check <- data.frame(
  id = aligned_school$schoolid,
  H,
  check.names = FALSE
) %>%
  group_by(id) %>%
  summarise(
    across(everything(), n_distinct),
    .groups = "drop"
  )

stopifnot(
  all(as.matrix(H_check[, -1, drop = FALSE]) == 1)
)

Tr <- aligned_outcome$Tr
y_raw <- aligned_outcome$y_raw

# Standardize using the control group in the final estimation sample.  The
# treatment effects are therefore measured in control-group SD units.
control_mean <- mean(y_raw[Tr == 0])
control_sd <- sd(y_raw[Tr == 0])

stopifnot(
  sum(Tr == 0) >= 2L,
  is.finite(control_mean),
  is.finite(control_sd),
  control_sd > 0
)

y <- (y_raw - control_mean) / control_sd

# Reindex schools only after imposing all estimation-sample restrictions.
# This ensures that the Stan indices are consecutive and every j in 1:J occurs.
id_factor <- factor(aligned_school$schoolid)
id <- as.integer(id_factor)
J <- nlevels(id_factor)

stopifnot(
  length(id) == NJ,
  identical(sort(unique(id)), seq_len(J)),
  length(y) == NJ,
  length(Tr) == NJ,
  nrow(X) == NJ,
  nrow(Z) == NJ,
  nrow(W) == NJ,
  nrow(H) == NJ
)

# A convenient combined data frame for descriptive checks and OLS benchmarks.
Data <- data.frame(
  id = id,
  original_id = ids_common,
  y = y,
  y_raw = y_raw,
  Tr = Tr,
  X,
  W,
  Z,
  H,
  check.names = FALSE
)

# Constants needed to report standardized effects or recover the raw-score
# scale later.
score_standardization <- list(
  grade = selected_grade,
  subject = unique(STAR_selected$selected_subject),
  reference_group = "Tr == 0 in the final estimation sample",
  mean = control_mean,
  sd = control_sd
)

Kx <- ncol(X) 
Kw <- ncol(W)
Kz <- ncol(Z)
Kh <- ncol(H)

# Construct the OLS benchmark from the columns actually present.
# x1 is excluded because lm() includes an intercept by default.
ols_predictors <- c(
  "Tr",
  colnames(X)[-1],
  colnames(W)
)

ols_formula <- reformulate(
  termlabels = ols_predictors,
  response = "y"
)

RegOLS <- lm(
  formula = ols_formula,
  data = Data
)

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

# delta is theta in our notation, and theta is tau in our notation.
# The model is hierarchical in delta (theta) with a multivariate normal
# distribution whose mean depends on study-level covariates H and whose
# covariance matrix is Sigma_delta. The treatment effect for student i in
# study j is tau_ij = z_i' * delta_j, where z_i are the individual-level
# covariates for student i.
# tau is psi in our notation

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
  beta               ~ normal(0, 10);
  gamma              ~ normal(0, 10);
  to_vector(Gamma)   ~ normal(0, 10);

  tau                ~ cauchy(0, 2.5);
  L_Omega            ~ lkj_corr_cholesky(2.0);
  to_vector(delta_raw) ~ std_normal();

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
options(mc.cores = parallel::detectCores()-4)

fit <- stan(
  model_code = stan_hier_meta_W_gamma_hetero,
  data   = stan_data,
  # Multiple independent chains are required for split-R-hat and for
  # diagnosing whether the chains explore the same posterior distribution.
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 1234,
  control = list(adapt_delta = 0.95)
)

grade_file_label <- c(
  kindergarten = "Kinder",
  first = "First",
  second = "Second",
  third = "Third"
)

subject_file_label <- c(
  reading = "Reading",
  math = "Math"
)

grade_selected <- unique(STAR_selected$selected_grade)
subject_selected <- unique(STAR_selected$selected_subject)

stopifnot(
  length(grade_selected) == 1L,
  length(subject_selected) == 1L,
  grade_selected %in% names(grade_file_label),
  subject_selected %in% names(subject_file_label)
)

fit_filename <- paste0(
  "AppSTAR_HeteroV2",
  unname(grade_file_label[grade_selected]),
  unname(subject_file_label[subject_selected]),
  ".RData"
)

save(fit, file = fit_filename)

message("Stan fit saved to: ", fit_filename)

# load(file = "AppSTAR_HeteroV2KinderReading.RData")
# load(file = "AppSTAR_HeteroV2KinderMath.RData")
# post <- summary(fit, pars = c("beta", "gamma", "Gamma", "Sigma_delta", "sigma"),
#                 probs = c(0.025, 0.5, 0.975))

post <- rstan::extract(
  fit,
  pars = c("beta", "gamma", "Gamma", "Sigma_delta", "sigma")
)

colnames(post[["beta"]]) <- colnames(X)
colnames(post[["gamma"]]) <- colnames(W)

Gamma_1 <- post[["Gamma"]][, , 1]
Gamma_2 <- post[["Gamma"]][, , 2]

colnames(Gamma_1) <- colnames(H)
colnames(Gamma_2) <- colnames(H)

summary(coda::mcmc(post[["beta"]]))
summary(coda::mcmc(post[["gamma"]]))

cat("\nGamma for:", colnames(Z)[1], "\n")
summary(coda::mcmc(Gamma_1))

cat("\nGamma for:", colnames(Z)[2], "\n")
summary(coda::mcmc(Gamma_2))

summary(coda::mcmc(post[["Sigma_delta"]][,,1]))
summary(coda::mcmc(post[["Sigma_delta"]][,,2]))
summary(coda::mcmc(matrix(post[["sigma"]], ncol = 1)))

Gamma_1 <- post[["Gamma"]][, , 1]

Gamma_1_summary <- t(
  apply(
    Gamma_1,
    MARGIN = 2,
    FUN = function(x) {
      c(
        posterior_mean = mean(x),
        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975)),
        Pr_positive = mean(x > 0),
        Pr_negative = mean(x < 0)
      )
    }
  )
)

Gamma_1_summary <- as.data.frame(Gamma_1_summary)

if (!is.null(colnames(H))) {
  Gamma_1_summary$H_covariate <- colnames(H)
} else {
  Gamma_1_summary$H_covariate <- paste0(
    "h",
    seq_len(nrow(Gamma_1_summary))
  )
}

Gamma_1_summary <- Gamma_1_summary %>%
  select(
    H_covariate,
    posterior_mean,
    lower_95,
    upper_95,
    Pr_positive,
    Pr_negative
  )

Gamma_1_summary

Gamma_1_summary %>%
  mutate(
    Pr_direction = pmax(Pr_positive, Pr_negative)
  ) %>%
  arrange(desc(Pr_direction))

Gamma_2 <- post[["Gamma"]][, , 2]

Gamma_2_summary <- t(
  apply(
    Gamma_2,
    MARGIN = 2,
    FUN = function(x) {
      c(
        posterior_mean = mean(x),
        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975)),
        Pr_positive = mean(x > 0),
        Pr_negative = mean(x < 0)
      )
    }
  )
)

Gamma_2_summary <- as.data.frame(Gamma_2_summary)

if (!is.null(colnames(H))) {
  Gamma_2_summary$H_covariate <- colnames(H)
} else {
  Gamma_2_summary$H_covariate <- paste0(
    "h",
    seq_len(nrow(Gamma_2_summary))
  )
}

Gamma_2_summary <- Gamma_2_summary %>%
  select(
    H_covariate,
    posterior_mean,
    lower_95,
    upper_95,
    Pr_positive,
    Pr_negative
  )

Gamma_2_summary

Gamma_2_summary %>%
  mutate(
    Pr_direction = pmax(Pr_positive, Pr_negative)
  ) %>%
  arrange(desc(Pr_direction))


# Report convergence and HMC diagnostics.  In particular, inspect split-R-hat,
# bulk/tail effective sample sizes, divergences, treedepth, and E-BFMI before
# interpreting posterior summaries.
print(
  fit,
  pars = c("beta", "gamma", "Gamma", "Sigma_delta", "sigma"),
  probs = c(0.025, 0.5, 0.975)
)
rstan::check_hmc_diagnostics(fit)

rstan::stan_trace(fit, pars = "beta")
rstan::stan_trace(fit, pars = "gamma")
rstan::stan_trace(fit, pars = "Gamma")
rstan::stan_trace(fit, pars = "Sigma_delta")
rstan::stan_trace(fit, pars = "sigma")

##### Trace plots #####
library(ggplot2)
library(dplyr)
library(tidyr)

# Unlike rstan::extract(..., permuted = TRUE), as.array() retains the actual
# iteration and chain dimensions.  This ordering is essential for trace plots.
mcmc_to_df <- function(fit, pars) {
  draws <- as.array(fit, pars = pars)
  draw_names <- dimnames(draws)

  stopifnot(length(dim(draws)) == 3L)

  expand.grid(
    iter = seq_len(dim(draws)[1]),
    chain = seq_len(dim(draws)[2]),
    param_index = seq_len(dim(draws)[3])
  ) %>%
    mutate(
      value = as.vector(draws),
      param = draw_names[[3]][param_index],
      chain = factor(chain)
    ) %>%
    select(iter, chain, param, value)
}

traceplot_gg <- function(df, title = NULL) {
  ggplot(df, aes(x = iter, y = value, color = chain, group = chain)) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    facet_wrap(~ param, scales = "free_y") +
    labs(
      title = title,
      x = "Iteration",
      y = "Value",
      color = "Chain"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

df_beta <- mcmc_to_df(fit, "beta")
p_beta <- traceplot_gg(df_beta, "Trace plots: beta")
ggsave("appendix_trace_beta_v1.pdf", p_beta, width = 8, height = 10)

df_gamma <- mcmc_to_df(fit, "gamma")
p_gamma <- traceplot_gg(df_gamma, "Trace plots: gamma")
ggsave("appendix_trace_gamma.pdf", p_gamma, width = 8, height = 10)

for (k in seq_len(Kz)) {
  gamma_pars <- paste0("Gamma[", seq_len(Kh), ",", k, "]")
  df_Gamma <- mcmc_to_df(fit, gamma_pars)
  p <- traceplot_gg(df_Gamma, paste0("Trace plots: Gamma[,", k, "]"))
  ggsave(paste0("appendix_trace_Gamma_v1_", k, ".pdf"),
         p, width = 8, height = 10)
}

for (k in seq_len(Kz)) {
  sigma_delta_pars <- paste0("Sigma_delta[", seq_len(Kz), ",", k, "]")
  df_Sd <- mcmc_to_df(fit, sigma_delta_pars)
  p <- traceplot_gg(
    df_Sd,
    paste0("Trace plots: Sigma_delta[,", k, "]")
  )
  ggsave(paste0("appendix_trace_Sigma_delta_", k, ".pdf"),
         p, width = 8, height = 10)
}

df_sigma <- mcmc_to_df(fit, "sigma")
p_sigma <- traceplot_gg(df_sigma, "Trace plot: sigma")
ggsave("appendix_trace_sigma_v1.pdf", p_sigma, width = 6, height = 4)

# ##### Nice plots of tau (treatment effects) #####
# library(dplyr)
# library(purrr)
# library(ggplot2)

# post_delta <- rstan::extract(fit, pars = "delta")$delta  # iter x J x Kz

# # choose 10 students (change these indices)
# Ni <- 10
# idx <- sample(1:NJ, Ni, replace = FALSE)

# # posterior draws of tau_ij (theta) for each selected student
# posterior_tau_i <- function(i, post_delta, stan_data) {
#   j_i <- as.integer(stan_data$id[i])
#   z_i <- as.numeric(stan_data$Z[i, ])

#   S <- dim(post_delta)[1]
#   Kz <- stan_data$Kz

#   stopifnot(
#     j_i >= 1L,
#     j_i <= dim(post_delta)[2],
#     length(z_i) == Kz,
#     dim(post_delta)[3] == Kz
#   )

#   # Retain an S x Kz matrix even when Kz = 1
#   delta_j_draws <- matrix(
#     post_delta[, j_i, , drop = FALSE],
#     nrow = S,
#     ncol = Kz
#   )

#   drop(delta_j_draws %*% z_i)
# }

# draws_df <- purrr::map_dfr(idx, function(i) {
#   tau_draws <- posterior_tau_i(
#     i = i,
#     post_delta = post_delta,
#     stan_data = stan_data
#   )

#   tibble::tibble(
#     student = factor(
#       paste0("i = ", i),
#       levels = paste0("i = ", idx)
#     ),
#     tau = tau_draws
#   )
# })

# # mean + 95% equal-tail CI
# sum_df <- draws_df %>%
#   group_by(student) %>%
#   summarise(
#     mean = mean(tau),
#     lo   = quantile(tau, 0.025),
#     hi   = quantile(tau, 0.975),
#     .groups = "drop"
#   )

# # histogram + vlines, 5x2 panels
# ggplot(draws_df, aes(x = tau)) +
#   geom_histogram(
#     bins = 60,
#     boundary = 0,
#     closed = "left",
#     fill = "grey85",
#     color = "white",
#     linewidth = 0.25
#   ) +
#   geom_vline(data = sum_df, aes(xintercept = mean, linetype = "Mean"),
#              linewidth = 1.0, color = "black") +
#   geom_vline(data = sum_df, aes(xintercept = lo, linetype = "95% CI"),
#              linewidth = 0.8, color = "black") +
#   geom_vline(data = sum_df, aes(xintercept = hi, linetype = "95% CI"),
#              linewidth = 0.8, color = "black") +
#   facet_wrap(~ student, ncol = 2, scales = "free_y") +
#   scale_linetype_manual(values = c("95% CI" = "dashed", "Mean" = "solid")) +
#   labs(x = expression(tau[ij]), y = "Count", linetype = "") +
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "bottom",
#     panel.grid = element_blank()
#   ) + 
#   scale_x_continuous(
#     breaks = seq(-60, 60, by = 10),
#     minor_breaks = seq(-2, 4, by = 0.25)
#   )

# ggsave(
#   filename = "school_treatment_effects.pdf",
#   plot = last_plot(),
#   width = 8,
#   height = 10
# )


#### Students in a particular school ####
# Select a school using its Stan index j
# j_selected <- 128L

# Row indices of all students belonging to that school
posterior_school_average_tau <- function(
    j_selected,
    stan_data,
    delta_draws
) {
  school_rows <- which(stan_data$id == j_selected)

  if (length(school_rows) == 0L) {
    stop("The selected school does not occur in stan_data$id.")
  }

  S <- dim(delta_draws)[1]
  J <- dim(delta_draws)[2]
  Kz <- stan_data$Kz

  stopifnot(
    j_selected >= 1L,
    j_selected <= J,
    dim(delta_draws)[3] == Kz,
    ncol(stan_data$Z) == Kz
  )

  # Average student-level treatment-effect modifiers in the school
  z_bar <- colMeans(
    stan_data$Z[school_rows, , drop = FALSE]
  )

  # Posterior draws of the selected school's random coefficients
  delta_selected <- matrix(
    delta_draws[, j_selected, , drop = FALSE],
    nrow = S,
    ncol = Kz
  )

  # Posterior draws of the observed school's average treatment effect
  tau_draws <- drop(
    delta_selected %*% z_bar
  )

  list(
    school = j_selected,
    student_rows = school_rows,
    n_students = length(school_rows),
    z_bar = z_bar,
    delta_draws = delta_selected,
    tau_draws = tau_draws
  )
}

j_selected <- 1L # School 25 is different compared to schools 1, 5, and 128

post_delta <- rstan::extract(
  fit,
  pars = "delta"
)$delta

school_posterior <- posterior_school_average_tau(
  j_selected = j_selected,
  stan_data = stan_data,
  delta_draws = post_delta
)

tau_school <- school_posterior$tau_draws

quantile(
  tau_school,
  probs = c(0.025, 0.5, 0.975)
)


# ##### No use ##############

# school_prediction <- predict_school_average_tau(
#   j_selected = j_selected,
#   stan_data = stan_data,
#   Gamma_draws = post[["Gamma"]],
#   Sigma_draws = post[["Sigma_delta"]],
#   seed = 1234
# )

# tau_school <- school_prediction$tau_draws

# school_tau_summary <- c(
#   posterior_mean = mean(tau_school),
#   posterior_median = median(tau_school),
#   lower_95 = unname(quantile(tau_school, 0.025)),
#   upper_95 = unname(quantile(tau_school, 0.975)),
#   Pr_positive = mean(tau_school > 0),
#   Pr_negative = mean(tau_school < 0)
# )

# school_tau_summary
# ##################################

library(ggplot2)

p_school_tau <- ggplot(
  data.frame(tau = as.numeric(tau_school)),
  aes(x = tau)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    fill = "grey70",
    color = "white",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    title = paste(
      "Posterior treatment effect: observed school",
      j_selected
    ),
    subtitle = paste(
      "Number of students =",
      school_posterior$n_students
    ),
    x = expression(tau[j]),
    y = "Posterior density"
  ) +
  theme_classic()

# Display on screen
grDevices::windows(
  width = 10,
  height = 7
)

print(p_school_tau)

# Save as a vector PDF
ggsave(
  filename = paste0(
    "p_school_tau_school_",
    j_selected,
    ".pdf"
  ),
  plot = p_school_tau,
  width = 8,
  height = 6,
  device = cairo_pdf
)

# All schools in the estimation sample
school_ids <- sort(unique(stan_data$id))

post_delta <- rstan::extract(
  fit,
  pars = "delta"
)$delta

school_posteriors <- lapply(
  school_ids,
  function(j) {
    posterior_school_average_tau(
      j_selected = j,
      stan_data = stan_data,
      delta_draws = post_delta
    )
  }
)

names(school_posteriors) <- paste0(
  "school_",
  school_ids
)

# Rows = posterior draws; columns = observed schools
tau_schools <- vapply(
  school_posteriors,
  function(x) x$tau_draws,
  numeric(dim(post_delta)[1])
)

colnames(tau_schools) <- paste0(
  "school_",
  school_ids
)

dim(tau_schools)

school_tau_summary <- t(
  vapply(
    seq_along(school_ids),
    function(j) {
      x <- tau_schools[, j]

      c(
        posterior_mean = mean(x),
        posterior_median = median(x),
        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975)),
        Pr_positive = mean(x > 0),
        Pr_negative = mean(x < 0)
      )
    },
    numeric(6)
  )
)

school_tau_summary <- data.frame(
  school = school_ids,
  school_tau_summary,
  row.names = NULL,
  check.names = FALSE
)

school_tau_summary

# ###### No use
school_ids <- sort(unique(stan_data$id))

set.seed(1234)  # Set once, outside the loop

school_posteriors <- lapply(
  school_ids,
  function(j) {
    posterior_school_average_tau(
      j_selected = j,
      stan_data = stan_data,
      delta_draws = post_delta
    )
  }
)

names(school_posteriors) <- paste0(
  "school_",
  school_ids
)

# names(school_predictions) <- paste0("school_", school_ids)

# # Posterior draws: rows = draws, columns = schools
# tau_school <- sapply(
#   school_predictions, 
#   function(x) x$tau_draws
# )

# colnames(tau_school) <- paste0("school_", school_ids)

# # Summary of posterior draws for each school
# school_tau_summary <- t(
#   apply(tau_school, 2, function(x) {
#     c(
#       posterior_mean = mean(x),
#       posterior_median = median(x),
#       lower_95 = unname(quantile(x, 0.025)),
#       upper_95 = unname(quantile(x, 0.975)),
#       Pr_positive = mean(x > 0),
#       Pr_negative = mean(x < 0)
#     )
#   })
# )

# school_tau_summary <- data.frame(
#   school = school_ids,
#   school_tau_summary,
#   row.names = NULL
# )

# school_tau_summary
##################################################

# Histogram posterior means all schools
library(ggplot2)

p_ate_schools <- ggplot(
  school_tau_summary,
  aes(x = posterior_mean)
) +
  geom_histogram(
    bins = 25,
    fill = "grey75",
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
    labs(
    title = "Posterior distribution average treatment effect (ATE) across schools",
    x = "Posterior mean of school ATE",
    y = "Number of schools"
  ) +
  theme_classic()

dev.off()

grDevices::windows(
  width = 10,
  height = 7
)
print(p_ate_schools)

ggsave(
  filename = "p_ate_schools.pdf",
  plot = p_ate_schools,
  width = 8,
  height = 6
)

##### Predictive posterior draws of tau_j* same average regressors as selected shool j but outside sample

selected_rows <- which(stan_data$id == j_selected)

if (length(selected_rows) == 0L) {
  stop("j_selected is not present in stan_data$id.")
}

h_new <- as.numeric(
  stan_data$H[selected_rows[1], ]
)

zbar_new <- colMeans(
  stan_data$Z[selected_rows, , drop = FALSE]
)

n_draws <- dim(post[["Gamma"]])[1]
Kh <- stan_data$Kh
Kz <- stan_data$Kz

set.seed(10101)

tau_new_school <- vapply(
  seq_len(n_draws),
  function(s) {
    Gamma_s <- matrix(
      post[["Gamma"]][s, , , drop = FALSE],
      nrow = Kh,
      ncol = Kz
    )

    Sigma_s <- matrix(
      post[["Sigma_delta"]][s, , , drop = FALSE],
      nrow = Kz,
      ncol = Kz
    )

    mean_delta_s <- drop(
      t(Gamma_s) %*% h_new
    )

    delta_new_s <- MASS::mvrnorm(
      n = 1,
      mu = mean_delta_s,
      Sigma = Sigma_s
    )

    drop(crossprod(zbar_new, delta_new_s))
  },
  numeric(1)
)

new_school_summary <- c(
  posterior_mean = mean(tau_new_school),
  posterior_median = median(tau_new_school),
  lower_95 = unname(quantile(tau_new_school, 0.025)),
  upper_95 = unname(quantile(tau_new_school, 0.975)),
  Pr_positive = mean(tau_new_school > 0),
  Pr_above_015 = mean(tau_new_school > 0.15)
)

new_school_summary

# # Rows belonging to the selected observed school
# selected_rows <- which(id == j_selected)

# if (length(selected_rows) == 0L) {
#   stop("j_selected is not present in id.")
# }

# # New school has the same H covariates and average Z composition
# # as the selected school
# h_new <- H[selected_rows[1], ]

# zbar_new <- colMeans(
#   Z[selected_rows, , drop = FALSE]
# )

# n_draws <- dim(post[["Gamma"]])[1]

# set.seed(10101)

# tau_new_school <- vapply(
#   seq_len(n_draws),
#   function(s) {
#     Gamma_s <- post[["Gamma"]][s, , ]
#     Sigma_s <- post[["Sigma_delta"]][s, , ]

#     mean_delta_s <- drop(
#       t(Gamma_s) %*% h_new
#     )

#     delta_new_s <- MASS::mvrnorm(
#       n = 1,
#       mu = mean_delta_s,
#       Sigma = Sigma_s
#     )

#     drop(
#       crossprod(zbar_new, delta_new_s)
#     )
#   },
#   numeric(1)
# )

# # Posterior predictive summary
# quantile(
#   tau_new_school,
#   probs = c(0.025, 0.5, 0.975)
# )

# mean(tau_new_school)
# mean(tau_new_school > 0)
# mean(tau_new_school > 0.15)

# Posterior predictive distribution
tau_comparison <- bind_rows(
  data.frame(
    tau = as.numeric(tau_school),
    distribution = "Observed school"
  ),
  data.frame(
    tau = as.numeric(tau_new_school),
    distribution = "New school with same regressors"
  )
)

p_tau_comparison <- ggplot(
  tau_comparison,
  aes(
    x = tau,
    fill = distribution
  )
) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    bins = 40,
    alpha = 0.55,
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Observed school" = "grey35",
      "New school with same regressors" = "grey75"
    )
  ) +
  labs(
    title = "Observed-school posterior versus new-school prediction",
    subtitle = paste(
      "Both schools have the same regressors as school",
      j_selected
    ),
    x = expression(tau[j]),
    y = "Posterior density",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

print(p_tau_comparison)

ggsave(
  filename = paste0(
    "tau_comparison_school_",
    j_selected,
    ".pdf"
  ),
  plot = p_tau_comparison,
  width = 8,
  height = 6,
  device = cairo_pdf
)

comparison_summary <- data.frame(
  school = c("Observed school", "New school"),
  sd = c(
    sd(tau_school),
    sd(tau_new_school)
  ),
  interval_width_95 = c(
    diff(quantile(tau_school, c(0.025, 0.975))),
    diff(quantile(tau_new_school, c(0.025, 0.975)))
  )
)

comparison_summary


############ Predictive external validity R2 ############
# -------------------------------------------------------------------------
# Empirical external-validity decomposition
# -------------------------------------------------------------------------

stopifnot(
  length(id) == nrow(H),
  length(id) == nrow(Z)
)

# Rows corresponding to each observed school
school_rows <- split(
  seq_along(id),
  factor(id, levels = unique(id))
)

school_ids <- unique(id)
J_observed <- length(school_rows)

# One H vector per school
H_school <- do.call(
  rbind,
  lapply(
    school_rows,
    function(rows) H[rows[1], , drop = FALSE]
  )
)

# Average Z vector within each school
Zbar_school <- do.call(
  rbind,
  lapply(
    school_rows,
    function(rows) {
      colMeans(Z[rows, , drop = FALSE])
    }
  )
)

rownames(H_school) <- school_ids
rownames(Zbar_school) <- school_ids

# Confirm dimensions
stopifnot(
  nrow(H_school) == J_observed,
  nrow(Zbar_school) == J_observed,
  ncol(H_school) == dim(post[["Gamma"]])[2],
  ncol(Zbar_school) == dim(post[["Gamma"]])[3],
  ncol(Zbar_school) == dim(post[["Sigma_delta"]])[2]
)

n_draws <- dim(post[["Gamma"]])[1]

EV_draws <- t(
  vapply(
    seq_len(n_draws),
    function(s) {
      Gamma_s <- post[["Gamma"]][s, , ]
      Sigma_s <- post[["Sigma_delta"]][s, , ]

      # E[tau_j | H_j, Zbar_j, Gamma_s]
      conditional_mean <- rowSums(
        (H_school %*% Gamma_s) * Zbar_school
      )

      # Empirical-population variance: denominator J, not J - 1
      V_sys <- mean(
        (conditional_mean - mean(conditional_mean))^2
      )

      # Zbar_j' Sigma_s Zbar_j for every school
      residual_variance <- rowSums(
        (Zbar_school %*% Sigma_s) * Zbar_school
      )

      V_res <- mean(residual_variance)

      R_EV2 <- V_sys / (V_sys + V_res)

      c(
        V_sys = V_sys,
        V_res = V_res,
        V_total = V_sys + V_res,
        R_EV2 = R_EV2
      )
    },
    numeric(4)
  )
)

EV_draws <- as.data.frame(EV_draws)

EV_summary <- t(
  vapply(
    EV_draws,
    function(x) {
      c(
        posterior_mean = mean(x),
        posterior_median = median(x),
        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975))
      )
    },
    numeric(4)
  )
)

EV_summary

p_R_EV2 <- ggplot(
  EV_draws,
  aes(x = R_EV2)
) +
  geom_histogram(
    bins = 40,
    fill = "grey50",
    color = "white"
  ) +
  labs(
    title = expression(
      "Posterior distribution of " * R[EV]^2
    ),
    subtitle = paste(
      "Empirical population of",
      J_observed,
      "observed schools"
    ),
    x = expression(R[EV]^2),
    y = "Posterior draws"
  ) +
  theme_classic()

print(p_R_EV2)
ggsave(
  filename = "posterior_R_EV2.pdf",
  plot = p_R_EV2,
  width = 8,
  height = 6
)

############ Up to this point for the paper ############
###########################################################

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
idx <- seq_len(NJ)
draws_dfAll <- sapply(seq_len(NJ), function(i){
  j_i <- stan_data$id[i]
  z_i <- stan_data$Z[i, ]                          
  tau_draws <- as.vector(post_delta[, j_i, ] %*% z_i)
  c(mean(tau_draws), quantile(tau_draws, c(0.025, 0.975)))
})

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

# Use all retained posterior draws rather than assuming a particular number of
# chains or post-warmup iterations.
S <- dim(post[["Gamma"]])[1]

# The school profile is rural and otherwise matches student 2623's observed
# school covariates, making the top-left predictive/conditional comparison
# interpretable.  The lower panels change average experience only.
x10 <- c(0, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(7.95, 1, 1, 0, 0) # (Av. Exp, Prop. Cau., Rural, Suburban, Urban)
tauij1 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                               Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(1, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(7.95, 1, 1, 0, 0) # (Av. Exp, Prop. Cau., Rural, Suburban, Urban)
tauij2 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                            Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(0, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(4, 1, 1, 0, 0) # (Av. Exp, Prop. Cau., Rural, Suburban, Urban)
tauij3 <- sapply(1:S, function(s) {Pred_tau(x10, w10, Gammas = post[["Gamma"]][s,,], 
                                            Sigmas = post[["Sigma_delta"]][s,,])})

x10 <- c(1, 1, 0) # (Female, caucasian, Afro. Amer.)
w10 <- c(4, 1, 1, 0, 0) # (Av. Exp, Prop. Cau., Rural, Suburban, Urban)
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
      "Female, Exp = 8",
      "Male, Exp = 4",
      "Female, Exp = 4"),
    each = S
  ), levels = c("Male, Exp = 8", "Female, Exp = 8",
                "Male, Exp = 4", "Female, Exp = 4"))
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
