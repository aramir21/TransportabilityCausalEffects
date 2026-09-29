# Two-stage OLS analysis of school-specific treatment effects
#
# Prerequisite:
# Run the grade/subject preparation code first so that Data, X, and H exist.
# In particular, Data must contain id, y, Tr, the columns of X, and the
# columns of H.
#
# Stage 1:
#   Within each school j, estimate
#
#     y_ij = alpha_j + Tr_ij * tau_j + x_ij' beta_j + error_ij.
#
# Stage 2:
#   Regress the estimated school effects on school characteristics:
#
#     tau_hat_j = h_j' Gamma + second-stage error_j.
#
# This is a descriptive two-stage OLS exercise. It does not propagate the
# first-stage estimation uncertainty as fully as the hierarchical Stan model.

# ---------------------------------------------------------------------------
# Validate the required objects
# ---------------------------------------------------------------------------

required_objects <- c("Data", "X", "H")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0L) {
  stop(
    "Run the grade/subject preparation code first. Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

required_data_columns <- c("id", "y", "Tr")
missing_data_columns <- setdiff(required_data_columns, names(Data))

if (length(missing_data_columns) > 0L) {
  stop(
    "Data is missing required columns: ",
    paste(missing_data_columns, collapse = ", ")
  )
}

if (nrow(Data) != nrow(X) || nrow(Data) != nrow(H)) {
  stop("Data, X, and H must have the same number and ordering of rows.")
}

# ---------------------------------------------------------------------------
# Stage 1: estimate a treatment effect separately for each school
# ---------------------------------------------------------------------------

# X includes an intercept in its first column. The school-specific lm()
# regressions include their own intercept, so x1 must not be included again.
x_covariates <- colnames(X)[-1]
x_covariates <- intersect(x_covariates, names(Data))

stage1_formula <- reformulate(
  termlabels = c("Tr", x_covariates),
  response = "y"
)

estimate_school_effect <- function(school_data, formula) {
  school_id <- unique(school_data$id)

  if (length(school_id) != 1L) {
    stop("Each first-stage data set must contain exactly one school.")
  }

  # A within-school treatment effect cannot be estimated unless the school
  # contains observations in both treatment states.
  observed_treatment <- unique(school_data$Tr[!is.na(school_data$Tr)])

  if (length(observed_treatment) < 2L) {
    return(data.frame(
      id = school_id,
      tau_hat = NA_real_,
      tau_se = NA_real_,
      n_school = nrow(school_data),
      residual_df = NA_real_,
      first_stage_r_squared = NA_real_,
      status = "No within-school treatment variation"
    ))
  }

  fit_j <- tryCatch(
    lm(formula = formula, data = school_data),
    error = function(e) e
  )

  if (inherits(fit_j, "error")) {
    return(data.frame(
      id = school_id,
      tau_hat = NA_real_,
      tau_se = NA_real_,
      n_school = nrow(school_data),
      residual_df = NA_real_,
      first_stage_r_squared = NA_real_,
      status = conditionMessage(fit_j)
    ))
  }

  coefficient_table <- coef(summary(fit_j))

  if (!"Tr" %in% rownames(coefficient_table) ||
      !is.finite(coefficient_table["Tr", "Estimate"])) {
    return(data.frame(
      id = school_id,
      tau_hat = NA_real_,
      tau_se = NA_real_,
      n_school = nrow(school_data),
      residual_df = df.residual(fit_j),
      first_stage_r_squared = summary(fit_j)$r.squared,
      status = "Treatment coefficient not estimable"
    ))
  }

  data.frame(
    id = school_id,
    tau_hat = unname(coefficient_table["Tr", "Estimate"]),
    tau_se = unname(coefficient_table["Tr", "Std. Error"]),
    n_school = nrow(school_data),
    residual_df = df.residual(fit_j),
    first_stage_r_squared = summary(fit_j)$r.squared,
    status = "Estimated"
  )
}

data_by_school <- split(Data, Data$id, drop = TRUE)

school_effects <- do.call(
  rbind,
  lapply(
    data_by_school,
    estimate_school_effect,
    formula = stage1_formula
  )
)

rownames(school_effects) <- NULL

if (!any(school_effects$status == "Estimated")) {
  stop("No school-specific treatment effects could be estimated.")
}

# Report schools omitted from the second stage.
stage1_status <- table(school_effects$status)
print(stage1_status)

# ---------------------------------------------------------------------------
# Construct one H row per school
# ---------------------------------------------------------------------------

h_names <- colnames(H)

if (is.null(h_names) || any(h_names == "")) {
  h_names <- paste0("h", seq_len(ncol(H)))
}

H_individual <- data.frame(
  id = Data$id,
  H,
  check.names = FALSE
)
names(H_individual)[-1] <- h_names

# Verify the key model assumption that H is constant within each school.
H_unique_counts <- aggregate(
  H_individual[, h_names, drop = FALSE],
  by = list(id = H_individual$id),
  FUN = function(x) length(unique(x))
)

if (!all(as.matrix(H_unique_counts[, h_names, drop = FALSE]) == 1L)) {
  stop("At least one H covariate varies within a school.")
}

H_study <- H_individual[
  !duplicated(H_individual$id),
  c("id", h_names),
  drop = FALSE
]

# ---------------------------------------------------------------------------
# Stage 2: regress estimated treatment effects on H
# ---------------------------------------------------------------------------

stage2_data <- merge(
  school_effects,
  H_study,
  by = "id",
  all.x = TRUE,
  sort = TRUE
)

stage2_data <- stage2_data[
  stage2_data$status == "Estimated" &
    is.finite(stage2_data$tau_hat),
  ,
  drop = FALSE
]

# h1 is the intercept column of H. lm() supplies an intercept automatically,
# so only the remaining columns enter as explicit regressors.
h_covariates <- h_names[-1]

if (nrow(stage2_data) <= length(h_covariates) + 1L) {
  stop(
    "Too few schools for the requested second-stage regression: ",
    nrow(stage2_data), " usable schools and ",
    length(h_covariates) + 1L, " coefficients."
  )
}

stage2_formula <- reformulate(
  termlabels = h_covariates,
  response = "tau_hat"
)

RegOLS_stage2 <- lm(
  formula = stage2_formula,
  data = stage2_data
)

# Requested second-stage OLS result
print(summary(RegOLS_stage2))

# A sensitivity analysis using inverse first-stage sampling variances. This is
# not identical to a random-effects meta-regression, but it accounts partially
# for differences in the precision of the estimated school effects.
stage2_wls_data <- stage2_data[
  is.finite(stage2_data$tau_se) & stage2_data$tau_se > 0,
  ,
  drop = FALSE
]

RegWLS_stage2 <- NULL

if (nrow(stage2_wls_data) > length(h_covariates) + 1L) {
  RegWLS_stage2 <- lm(
    formula = stage2_formula,
    data = stage2_wls_data,
    weights = 1 / tau_se^2
  )

  message("\nInverse-variance weighted sensitivity analysis:")
  print(summary(RegWLS_stage2))
}

# ---------------------------------------------------------------------------
# Save results with grade- and subject-specific filenames when available
# ---------------------------------------------------------------------------

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

grade_selected <- if (exists("selected_grade", inherits = TRUE)) {
  unique(get("selected_grade", inherits = TRUE))
} else {
  "SelectedGrade"
}

subject_selected <- if (exists("STAR_selected", inherits = TRUE) &&
                        "selected_subject" %in% names(STAR_selected)) {
  unique(STAR_selected$selected_subject)
} else {
  "SelectedSubject"
}

grade_label <- if (length(grade_selected) == 1L &&
                   grade_selected %in% names(grade_file_label)) {
  unname(grade_file_label[grade_selected])
} else {
  "SelectedGrade"
}

subject_label <- if (length(subject_selected) == 1L &&
                     subject_selected %in% names(subject_file_label)) {
  unname(subject_file_label[subject_selected])
} else {
  "SelectedSubject"
}

result_file <- paste0(
  "STAR_TwoStageOLS_",
  grade_label,
  subject_label,
  ".rds"
)

summary_file <- paste0(
  "STAR_TwoStageOLS_",
  grade_label,
  subject_label,
  "_summary.txt"
)

two_stage_results <- list(
  first_stage_formula = stage1_formula,
  school_effects = school_effects,
  first_stage_status = stage1_status,
  second_stage_data = stage2_data,
  second_stage_formula = stage2_formula,
  second_stage_ols = RegOLS_stage2,
  second_stage_wls = RegWLS_stage2
)

saveRDS(two_stage_results, file = result_file)

capture.output(
  {
    cat("First-stage status\n")
    print(stage1_status)
    cat("\nSecond-stage OLS\n")
    print(summary(RegOLS_stage2))

    if (!is.null(RegWLS_stage2)) {
      cat("\nInverse-variance weighted sensitivity analysis\n")
      print(summary(RegWLS_stage2))
    }
  },
  file = summary_file
)

message("Two-stage results saved to: ", result_file)
message("Text summary saved to: ", summary_file)

# Identifing potential common factors that may explain the heterogeneity in treatment effects across schools
school_effects <- school_effects %>%
  mutate(
    t_stat = tau_hat / tau_se,
    p_value = 2 * pt(
      abs(t_stat),
      df = residual_df,
      lower.tail = FALSE
    ),

    # Adjustment for testing effects in many schools
    p_value_fdr = p.adjust(
      p_value,
      method = "BH"
    ),

    significant_5pct = p_value < 0.05,
    significant_fdr = p_value_fdr < 0.05,

    effect_group = case_when(
      significant_fdr & tau_hat > 0 ~ "Positive",
      significant_fdr & tau_hat < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  )

table(school_effects$effect_group)

x_names <- colnames(X)[-1]
w_names <- colnames(W)

school_profiles <- Data %>%
  group_by(id) %>%
  summarise(
    # Means of individual characteristics describe student composition
    across(
      all_of(x_names),
      ~ mean(.x, na.rm = TRUE),
      .names = "mean_{.col}"
    ),

    # W is constant within school, so take its first value
    across(
      all_of(w_names),
      first
    ),

    n_students = n(),
    treatment_share = mean(Tr),
    .groups = "drop"
  )

school_analysis <- school_effects %>%
  left_join(
    school_profiles,
    by = "id"
  ) %>%
  filter(
    status == "Estimated",
    is.finite(tau_hat),
    is.finite(tau_se),
    tau_se > 0
  )

profile_comparison <- school_analysis %>%
  group_by(effect_group) %>%
  summarise(
    across(
      c(
        starts_with("mean_x"),
        all_of(w_names)
      ),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    ),
    n_schools = n(),
    .groups = "drop"
  )

profile_comparison

school_predictors <- c(
  paste0("mean_", x_names),
  w_names
)

school_formula <- reformulate(
  termlabels = school_predictors,
  response = "tau_hat"
)

school_characteristics_ols <- lm(
  formula = school_formula,
  data = school_analysis
)

summary(school_characteristics_ols)

school_characteristics_wls <- lm(
  formula = school_formula,
  data = school_analysis,
  weights = 1 / tau_se^2
)

summary(school_characteristics_wls)

# install.packages("metafor")
library(metafor)

RegMeta_stage2 <- rma(
  yi = tau_hat,
  sei = tau_se,
  mods = ~ h2 + h3 + h4,
  method = "REML",
  data = stage2_data
)

summary(RegMeta_stage2)

second_stage_results <- as.data.frame(
  coef(summary(RegOLS_stage2)),
  check.names = FALSE
)

second_stage_results$term <- rownames(second_stage_results)
rownames(second_stage_results) <- NULL

second_stage_results <- second_stage_results %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  select(
    term,
    estimate,
    std_error,
    t_value,
    p_value
  )

second_stage_results$p_fdr <- NA_real_

slope_rows <- second_stage_results$term != "(Intercept)"

second_stage_results$p_fdr[slope_rows] <- p.adjust(
  second_stage_results$p_value[slope_rows],
  method = "BH"
)

second_stage_results <- second_stage_results %>%
  mutate(
    significant_5pct = p_value < 0.05,
    significant_fdr_5pct = p_fdr < 0.05
  )

second_stage_results

second_stage_results %>%
  filter(
    term != "(Intercept)",
    significant_fdr_5pct
  ) %>%
  arrange(p_fdr)
