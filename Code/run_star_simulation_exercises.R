# STAR super-population Monte Carlo design
# One new 80-school x 450-student master population is generated per
# replication. Nine nested counterfactual scenarios share that realization.

source("star_baseline_simulation.R")

star_scenario_table <- function() {
  data.frame(
    scenario = c("J_10", "J_20", "baseline", "J_80", "N_50", "N_450",
                 "c_0", "c_0.5", "c_2"),
    experiment = c(rep("number_of_schools", 4L),
                   rep("within_school_size", 2L),
                   rep("residual_heterogeneity", 3L)),
    J = c(10L, 20L, 40L, 80L, 40L, 40L, 40L, 40L, 40L),
    N_j = c(150L, 150L, 150L, 150L, 50L, 450L, 150L, 150L, 150L),
    c_scale = c(1, 1, 1, 1, 1, 1, 0, 0.5, 2),
    stringsAsFactors = FALSE
  )
}

generate_star_master <- function(
    schools, J_max = 80L, N_max = 450L,
    treatment_prob = 0.3, female_prob = 0.487,
    Gamma_0 = matrix(c(0.152, 0.059, 0.110, -0.066, 0.103, -0.067),
                     nrow = 2L, ncol = 3L),
    Sigma_0 = diag(0.124636, 2L), beta_0 = c(0, 0, 0),
    delta_0 = 0, sigma_u = 1) {
  stopifnot(J_max >= 80L, N_max >= 450L,
            treatment_prob >= 0, treatment_prob <= 1,
            female_prob >= 0, female_prob <= 1,
            all(dim(Gamma_0) == c(2L, 3L)),
            all(dim(Sigma_0) == c(2L, 2L)),
            length(beta_0) == 3L, sigma_u > 0)

  # Joint empirical draw of (H1, H2, W, pF), once per replication.
  school <- schools[sample.int(nrow(schools), J_max, replace = TRUE), ]
  rownames(school) <- NULL
  school$sim_school_id <- seq_len(J_max)

  # Standardized heterogeneity shared across all values of c.
  epsilon <- matrix(rnorm(2L * J_max), nrow = J_max, ncol = 2L)
  colnames(epsilon) <- c("epsilon_intercept", "epsilon_free_lunch")
  school <- cbind(school, epsilon)

  age_values <- c(-2.787070, -0.551523, 1.684024, 3.919572)
  age_probs <- c(0.001040, 0.757282, 0.235610, 0.006068)
  n_master <- J_max * N_max
  school_index <- rep(seq_len(J_max), each = N_max)
  student_index <- rep(seq_len(N_max), times = J_max)
  free_lunch <- rbinom(n_master, 1L, school$pF[school_index])
  female <- rbinom(n_master, 1L, female_prob)
  age <- sample(age_values, n_master, replace = TRUE, prob = age_probs)
  treatment <- rbinom(n_master, 1L, treatment_prob)
  u <- rnorm(n_master, 0, sigma_u)

  students <- data.frame(
    sim_school_id = school_index,
    source_school_id = school$school_id[school_index],
    student_id = student_index,
    H1 = school$prop_teacher_eth_cauc[school_index],
    H2 = school[["prop_degree_master+"]][school_index],
    W = school$avg_experience[school_index],
    pF = school$pF[school_index],
    free_lunch = free_lunch, female = female, age = age,
    treatment = treatment, u = u
  )

  list(
    schools = school, students = students,
    parameters = list(Gamma_0 = Gamma_0, Sigma_0 = Sigma_0,
                      beta_0 = beta_0, delta_0 = delta_0,
                      sigma_u = sigma_u,
                      treatment_prob = treatment_prob,
                      female_prob = female_prob)
  )
}

construct_star_scenario <- function(master, J, N_j, c_scale) {
  stopifnot(J <= 80L, N_j <= 450L, c_scale >= 0)
  school <- master$schools[seq_len(J), , drop = FALSE]
  keep <- master$students$sim_school_id <= J &
    master$students$student_id <= N_j
  dat <- master$students[keep, , drop = FALSE]
  rownames(dat) <- NULL
  par <- master$parameters

  # v_j(c) = sqrt(c) Sigma_0^(1/2) epsilon_j; epsilon is never redrawn.
  sigma_root <- t(chol(par$Sigma_0))
  epsilon <- as.matrix(school[, c("epsilon_intercept",
                                   "epsilon_free_lunch")])
  v <- sqrt(c_scale) * t(sigma_root %*% t(epsilon))
  h <- rbind(1, school$prop_teacher_eth_cauc,
             school[["prop_degree_master+"]])
  theta <- t(par$Gamma_0 %*% h) + v
  colnames(theta) <- c("theta_intercept", "theta_free_lunch")

  j <- dat$sim_school_id
  dat$tau <- theta[j, 1L] + dat$free_lunch * theta[j, 2L]
  dat$y <- dat$treatment * dat$tau + par$beta_0[1L] +
    par$beta_0[2L] * dat$female + par$beta_0[3L] * dat$age +
    par$delta_0 * dat$W + dat$u

  mu <- rowSums(cbind(1, school$pF) * t(par$Gamma_0 %*% h))
  sample_pF <- vapply(seq_len(J), function(k) {
    mean(dat$free_lunch[dat$sim_school_id == k])
  }, numeric(1))
  truth <- data.frame(
    sim_school_id = seq_len(J), source_school_id = school$school_id,
    c_scale = c_scale, mu_systematic = mu,
    tau_pop_true = theta[, 1L] + school$pF * theta[, 2L],
    tau_finite_true = theta[, 1L] + sample_pF * theta[, 2L],
    sample_pF = sample_pF,
    v_intercept = v[, 1L], v_free_lunch = v[, 2L],
    theta_intercept = theta[, 1L], theta_free_lunch = theta[, 2L]
  )
  list(data = dat, school_truth = truth)
}

construct_all_star_scenarios <- function(master, K = 10L) {
  settings <- star_scenario_table()
  if (K < 1L || K > min(settings$J)) {
    stop("K must be between 1 and ", min(settings$J), ".")
  }
  scenarios <- vector("list", nrow(settings))
  names(scenarios) <- settings$scenario
  for (s in seq_len(nrow(settings))) {
    q <- settings[s, ]
    scenario <- construct_star_scenario(master, q$J, q$N_j, q$c_scale)
    # These K schools are each omitted separately in LOSO, not jointly.
    scenario$loso_schools <- sample.int(q$J, K, replace = FALSE)
    scenario$settings <- q
    scenarios[[s]] <- scenario
  }
  scenarios
}

validate_common_random_numbers <- function(master, scenarios) {
  b <- scenarios$baseline$data
  same_values <- function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE))
  stopifnot(
    nrow(master$students) == 80L * 450L,
    nrow(b) == 40L * 150L,
    same_values(scenarios$J_10$data, b[b$sim_school_id <= 10L, ]),
    same_values(scenarios$J_20$data, b[b$sim_school_id <= 20L, ]),
    same_values(scenarios$N_50$data, b[b$student_id <= 50L, ]),
    same_values(scenarios$J_80$data[
      scenarios$J_80$data$sim_school_id <= 40L, ], b),
    same_values(scenarios$N_450$data[
      scenarios$N_450$data$student_id <= 150L, ], b),
    all(scenarios$c_0$school_truth$v_intercept == 0),
    all(scenarios$c_0$school_truth$v_free_lunch == 0)
  )
  fixed <- setdiff(names(b), c("tau", "y"))
  stopifnot(
    same_values(scenarios$c_0$data[fixed], b[fixed]),
    same_values(scenarios$c_0.5$data[fixed], b[fixed]),
    same_values(scenarios$c_2$data[fixed], b[fixed])
  )
  TRUE
}

run_star_exercises <- function(
    csv_file = "STAR_Data_Simulation.csv", R_rep = 100L, K = 10L,
    seed = 20260828L, analysis_fun = NULL, save_scenarios = FALSE,
    output_dir = "STAR_superpopulation_results",
    output_file = "STAR_superpopulation_summary_R100.rds") {
  stopifnot(R_rep >= 1L, K >= 1L)
  schools <- read_star_population(csv_file)
  settings <- star_scenario_table()
  if (K > min(settings$J)) stop("K cannot exceed 10.")
  if (save_scenarios && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  results <- vector("list", R_rep)

  for (r in seq_len(R_rep)) {
    # Regenerate the complete master population across replications.
    set.seed(seed + r - 1L)
    master <- generate_star_master(schools)
    scenarios <- construct_all_star_scenarios(master, K)
    validate_common_random_numbers(master, scenarios)

    if (is.function(analysis_fun)) {
      results[[r]] <- lapply(scenarios, function(q) analysis_fun(
        data = q$data, school_truth = q$school_truth,
        loso_schools = q$loso_schools, settings = q$settings,
        replication = r))
    } else {
      # Compact output when no estimation method has yet been supplied.
      results[[r]] <- lapply(scenarios, function(q) list(
        settings = q$settings, loso_schools = q$loso_schools,
        mean_tau_pop_true = mean(q$school_truth$tau_pop_true),
        sd_tau_pop_true = sd(q$school_truth$tau_pop_true),
        treatment_rate = mean(q$data$treatment),
        n_observations = nrow(q$data)))
    }

    # Optional streaming storage prevents accumulation of large datasets.
    if (save_scenarios) {
      saveRDS(list(replication = r, master = master, scenarios = scenarios),
              file.path(output_dir, sprintf("replication_%04d.rds", r)))
    }
    message("Completed replication ", r, " of ", R_rep)
  }

  output <- list(scenario_table = settings, R_rep = R_rep, K = K,
                 seed = seed, results = results,
                 saved_full_scenarios = save_scenarios,
                 scenario_directory = if (save_scenarios) output_dir else NULL)
  saveRDS(output, output_file)
  message("Saved Monte Carlo summary to: ",
          normalizePath(output_file, mustWork = FALSE))
  invisible(output)
}
