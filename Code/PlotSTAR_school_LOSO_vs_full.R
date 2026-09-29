# Histogram comparison of a school's latent ATE before and after observing it
#
# Requires ApplicationSTARV3_LOSO_metrics.R to have been run. This script does
# not estimate Stan models; it uses the saved latent-ATE comparison draws.

plot_config <- list(
  school_index = 1L,
  draws_path = file.path(
    "star_output",
    "loso_latent_metrics",
    "LOSO_vs_full_latent_ATE_draws.rds"
  ),
  output_dir = file.path("star_output", "loso_latent_metrics"),
  bins = 40L
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install the ggplot2 package before running this script.")
}
if (!file.exists(plot_config$draws_path)) {
  stop(
    "Latent comparison draws not found. Run ",
    "ApplicationSTARV3_LOSO_metrics.R first."
  )
}

comparison_draws <- readRDS(plot_config$draws_path)
available_schools <- vapply(comparison_draws, function(x) x$fold, integer(1))

if (!plot_config$school_index %in% available_schools) {
  stop("The selected school index is not present in the comparison draws.")
}

school_result <- comparison_draws[[
  match(plot_config$school_index, available_schools)
]]

plot_data <- rbind(
  data.frame(
    tau = as.numeric(school_result$loso_latent_ATE),
    distribution = "LOSO predictive: before observing school"
  ),
  data.frame(
    tau = as.numeric(school_result$full_posterior_latent_ATE),
    distribution = "Full posterior: after observing school"
  )
)

distribution_order <- c(
  "LOSO predictive: before observing school",
  "Full posterior: after observing school"
)
plot_data$distribution <- factor(
  plot_data$distribution,
  levels = distribution_order
)

mean_data <- stats::aggregate(
  tau ~ distribution,
  data = plot_data,
  FUN = mean
)

school_label <- paste0(
  "STAR school ID ", school_result$original_schoolid
)

p_school_loso_full <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = .data$tau, fill = .data$distribution)
) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    position = "identity",
    bins = plot_config$bins,
    alpha = 0.55,
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_vline(
    data = mean_data,
    ggplot2::aes(
      xintercept = .data$tau,
      linetype = .data$distribution
    ),
    color = "black",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "LOSO predictive: before observing school" = "grey75",
      "Full posterior: after observing school" = "grey30"
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "LOSO predictive: before observing school" = "dashed",
      "Full posterior: after observing school" = "solid"
    )
  ) +
  ggplot2::labs(
    title = paste(
      "Latent ATE prediction before and after observing school",
      plot_config$school_index
    ),
    subtitle = school_label,
    x = "School average treatment effect (control-group SD units)",
    y = "Posterior density",
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(face = "bold")
  )

dir.create(plot_config$output_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- sprintf(
  "school_%02d_LOSO_vs_full_latent_ATE",
  plot_config$school_index
)

ggplot2::ggsave(
  filename = file.path(plot_config$output_dir, paste0(output_stem, ".pdf")),
  plot = p_school_loso_full,
  width = 8,
  height = 6
)

utils::write.csv(
  data.frame(
    distribution = as.character(mean_data$distribution),
    posterior_mean = mean_data$tau,
    posterior_sd = vapply(
      split(plot_data$tau, plot_data$distribution),
      stats::sd,
      numeric(1)
    )[as.character(mean_data$distribution)],
    lower_95 = vapply(
      split(plot_data$tau, plot_data$distribution),
      function(x) unname(stats::quantile(x, 0.025)),
      numeric(1)
    )[as.character(mean_data$distribution)],
    upper_95 = vapply(
      split(plot_data$tau, plot_data$distribution),
      function(x) unname(stats::quantile(x, 0.975)),
      numeric(1)
    )[as.character(mean_data$distribution)]
  ),
  file = file.path(plot_config$output_dir, paste0(output_stem, "_summary.csv")),
  row.names = FALSE
)

print(p_school_loso_full)
message(
  "Figure written to: ",
  normalizePath(
    file.path(plot_config$output_dir, paste0(output_stem, ".pdf"))
  )
)
