# This block assumes that select_star_sample() has already been defined and
# that STAR_selected has been created, for example:
#
# STAR_selected <- select_star_sample(
#   data = STAR,
#   grade = "second",
#   subject = "math",
#   complete_cases = FALSE
# )

library(dplyr)
library(tidyr)
library(stringr)
library(rlang)

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
    birth_year = as.integer(str_extract(birth_character, "^\\d{4}")),
    age = unname(assessment_year[selected_grade]) - birth_year,
    Q1 = as.integer(str_detect(birth_character, "Q1")),
    Q2 = as.integer(str_detect(birth_character, "Q2")),
    Q3 = as.integer(str_detect(birth_character, "Q3")),
    Q4 = as.integer(str_detect(birth_character, "Q4")),
    eth_cauc = as.integer(ethnicity == "cauc"),
    eth_afam = as.integer(ethnicity == "afam"),
    female = as.integer(gender == "female")
  ) %>%
  select(id, female, eth_cauc, eth_afam, age, Q1, Q2, Q3) %>%
  drop_na()

# Q4, male, and ethnicities other than cauc/afam are the reference categories.

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

# Select covariates by name instead of relying on column positions.  Some
# categories (for example, ladder = "pending") may not occur in every grade.
# Such columns are omitted from W rather than filled with zeros, because an
# all-zero regressor has no identified likelihood contribution.  The excluded
# categories are the regression reference levels:
#   teacher ethnicity: non-caucasian categories;
#   ladder: apprentice;
#   degree: bachelor;
#   school location: inner city.
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
  "prop_school_urban"
)

missing_w_covariates <- setdiff(
  w_covariates_requested,
  names(teacher_stats_wide)
)

if (length(missing_w_covariates) > 0L) {
  message(
    "Omitting categories absent from the selected grade: ",
    paste(missing_w_covariates, collapse = ", ")
  )
}

W_study <- teacher_stats_wide %>%
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
  select(female, eth_cauc, eth_afam, age, Q1, Q2, Q3) %>%
  as.matrix()
X <- cbind(1, X)
colnames(X) <- paste0("x", seq_len(ncol(X)))

# Z determines treatment-effect heterogeneity: intercept, female, caucasian,
# and African American indicators.
Z <- X[, seq_len(4), drop = FALSE]
colnames(Z) <- paste0("z", seq_len(ncol(Z)))

# W contains school-level predictors of the outcome.  Its intercept is omitted
# because the common intercept is already included in X.
W <- aligned_school %>%
  select(-id, -schoolid) %>%
  as.matrix()
colnames(W) <- paste0("w", seq_len(ncol(W)))

# H contains school-level predictors of the heterogeneous treatment effects.
# The order is intercept, average experience, teacher caucasian proportion,
# and the rural/suburban/urban indicators.
h_covariates <- c(
  "avg_experience",
  "prop_teacher_eth_cauc",
  "prop_school_rural",
  "prop_school_suburban",
  "prop_school_urban"
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
H <- cbind(1, H)
colnames(H) <- paste0("h", seq_len(ncol(H)))

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

# Construct the OLS benchmark from the columns that are actually present.
# x1 is excluded because it is the intercept column; lm() includes an
# intercept by default.  This construction remains valid when a grade lacks an
# optional school-level category and W therefore has fewer columns.
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

# Constants needed to report standardized effects or recover the raw-score
# scale later.
score_standardization <- list(
  grade = selected_grade,
  subject = unique(STAR_selected$selected_subject),
  reference_group = "Tr == 0 in the final estimation sample",
  mean = control_mean,
  sd = control_sd
)
