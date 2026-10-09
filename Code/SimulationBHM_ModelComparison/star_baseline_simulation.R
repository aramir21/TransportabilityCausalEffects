# STAR-calibrated hierarchical treatment-effect simulation
#
# The CSV contains repeated column names (three columns named "intercept").
# They are intentionally ignored: only substantively named columns are used.

read_star_population <- function(csv_file = "STAR_Data_Simulation.csv") {
  dat <- read.csv(csv_file, check.names = FALSE)

  required <- c(
    "school_id", "prop_teacher_eth_cauc", "prop_degree_master+",
    "avg_experience", "free_lunch"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("CSV is missing required columns: ", paste(missing, collapse = ", "))
  }

  school_vars <- c(
    "prop_teacher_eth_cauc", "prop_degree_master+", "avg_experience"
  )
  varying <- vapply(
    school_vars,
    function(v) any(tapply(dat[[v]], dat$school_id, function(x) {
      length(unique(x[!is.na(x)])) != 1L
    })),
    logical(1)
  )
  if (any(varying)) {
    stop(
      "These nominally school-level variables vary within school: ",
      paste(names(varying)[varying], collapse = ", ")
    )
  }

  first <- !duplicated(dat$school_id)
  schools <- dat[first, c("school_id", school_vars)]
  pF <- aggregate(free_lunch ~ school_id, data = dat, FUN = mean)
  names(pF)[2L] <- "pF"
  schools <- merge(schools, pF, by = "school_id", sort = TRUE)

  if (anyNA(schools) || any(schools$pF < 0 | schools$pF > 1)) {
    stop("Invalid or missing values in the reconstructed school population.")
  }
  schools
}

simulate_star <- function(
    schools,
    J = 10L,
    N_j = 50L,
    treatment_prob = 0.3,
    Gamma_0 = matrix(
      c(0.152091, 0.058677,
        0.110434, -0.065810,
        0.103012, -0.067449),
      nrow = 2L, ncol = 3L
    ),
    Sigma_0 = diag(0.124636, 2L),
    beta_0 = c(0, 0, 0),
    delta_0 = 0,
    sigma_u = 1
) {
  stopifnot(
    J >= 1L, N_j >= 1L,
    treatment_prob >= 0, treatment_prob <= 1,
    all(dim(Gamma_0) == c(2L, 3L)),
    all(dim(Sigma_0) == c(2L, 2L)),
    length(beta_0) == 3L,
    sigma_u > 0
  )
  if (max(abs(Sigma_0 - t(Sigma_0))) > 1e-12) {
    stop("Sigma_0 must be symmetric.")
  }
  chol_Sigma <- tryCatch(
    chol(Sigma_0),
    error = function(e) stop("Sigma_0 must be positive definite.")
  )

  # Joint resampling preserves dependence among H1, H2, W, and pF.
  sampled <- schools[sample.int(nrow(schools), J, replace = TRUE), ]
  sampled$sim_school_id <- seq_len(J)

  age_values <- c(-2.787070, -0.551523, 1.684024, 3.919572)
  age_probs <- c(0.001040, 0.757282, 0.235610, 0.006068)

  student_list <- vector("list", J)
  truth_list <- vector("list", J)

  for (j in seq_len(J)) {
    h_j <- c(
      1,
      sampled$prop_teacher_eth_cauc[j],
      sampled[["prop_degree_master+"]][j]
    )
    v_j <- as.vector(t(chol_Sigma) %*% rnorm(2L))
    theta_j <- as.vector(Gamma_0 %*% h_j + v_j)

    free_lunch <- rbinom(N_j, 1L, sampled$pF[j])
    female <- rbinom(N_j, 1L, 0.486997)
    age <- sample(age_values, N_j, replace = TRUE, prob = age_probs)
    treatment <- rbinom(N_j, 1L, treatment_prob)
    tau <- theta_j[1L] + free_lunch * theta_j[2L]
    u <- rnorm(N_j, 0, sigma_u)
    y <- treatment * tau +
      beta_0[1L] + beta_0[2L] * female + beta_0[3L] * age +
      sampled$avg_experience[j] * delta_0 + u

    student_list[[j]] <- data.frame(
      sim_school_id = j,
      source_school_id = sampled$school_id[j],
      student_id = seq_len(N_j),
      H1 = sampled$prop_teacher_eth_cauc[j],
      H2 = sampled[["prop_degree_master+"]][j],
      W = sampled$avg_experience[j],
      pF = sampled$pF[j],
      free_lunch = free_lunch,
      female = female,
      age = age,
      treatment = treatment,
      tau = tau,
      y = y
    )

    mu_j <- as.numeric(c(1, sampled$pF[j]) %*% Gamma_0 %*% h_j)
    truth_list[[j]] <- data.frame(
      sim_school_id = j,
      source_school_id = sampled$school_id[j],
      mu_systematic = mu_j,
      tau_pop = sum(c(1, sampled$pF[j]) * theta_j),
      tau_true = sum(c(1, mean(free_lunch)) * theta_j),
      theta_intercept = theta_j[1L],
      theta_free_lunch = theta_j[2L]
    )
  }

  list(
    students = do.call(rbind, student_list),
    school_truth = do.call(rbind, truth_list),
    sampled_schools = sampled,
    settings = list(J = J, N_j = N_j, treatment_prob = treatment_prob)
  )
}

check_star_calibration <- function(
    schools,
    n_schools = 100000L,
    seed = 20260827L,
    Gamma_0 = matrix(
      c(0.152091, 0.058677,
        0.110434, -0.065810,
        0.103012, -0.067449),
      nrow = 2L, ncol = 3L
    ),
    Sigma_0 = diag(0.124636, 2L)
) {
  set.seed(seed)
  s <- schools[sample.int(nrow(schools), n_schools, replace = TRUE), ]
  h <- rbind(1, s$prop_teacher_eth_cauc, s[["prop_degree_master+"]])
  zbar <- cbind(1, s$pF)
  mu <- rowSums(zbar * t(Gamma_0 %*% h))
  residual_variance <- rowSums((zbar %*% Sigma_0) * zbar)
  v <- matrix(rnorm(2L * n_schools), ncol = 2L) %*% chol(Sigma_0)
  residual_ate <- rowSums(zbar * v)
  tau_pop <- mu + residual_ate

  data.frame(
    quantity = c(
      "mean_tau", "var_systematic", "mean_residual_variance",
      "R2_EV", "sd_tau", "cor_H1_pF"
    ),
    simulated = c(
      mean(tau_pop), var(mu), mean(residual_variance),
      var(mu) / (var(mu) + mean(residual_variance)),
      sd(tau_pop), cor(s$prop_teacher_eth_cauc, s$pF)
    ),
    target = c(0.1895, 0.008501, 0.161522, 0.05, 0.412339, -0.4526)
  )
}

# ---- Small smoke test ----------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(12345)
  star_schools <- read_star_population("STAR_Data_Simulation.csv")
  small <- simulate_star(
    schools = star_schools,
    J = 10L,
    N_j = 50L,
    treatment_prob = 0.3
  )

  stopifnot(
    nrow(star_schools) == 79L,
    nrow(small$students) == 10L * 50L,
    nrow(small$school_truth) == 10L,
    all(small$students$free_lunch %in% 0:1),
    all(small$students$treatment %in% 0:1),
    !anyNA(small$students),
    !anyNA(small$school_truth)
  )

  cat("Small simulation passed all structural checks.\n")
  print(small$settings)
  print(head(small$students))
  print(small$school_truth)

  cat("\nLarge-population calibration check:\n")
  print(check_star_calibration(star_schools))
}
