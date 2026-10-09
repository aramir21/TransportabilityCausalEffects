# Presentation of the oracle LOSO simulation results.
# Run with Rscript from any directory, or source from Code or SimulationBHM.
# Optional command-line arguments: metrics directory, output directory.
# Dependencies: ggplot2, dplyr, tidyr, readr, gtable (installed with ggplot2).

show_mc_bars <- TRUE
mc_multiplier <- stats::qnorm(0.975)

packages <- c("ggplot2", "dplyr", "tidyr", "readr", "gtable")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1),
                                   quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (!is.null(source_file)) {
  dirname(normalizePath(source_file))
} else if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (dir.exists("SimulationBHM/STAR_superpopulation_metrics")) {
  normalizePath("SimulationBHM")
} else {
  getwd()
}
input_dir <- if (length(args) >= 1L) args[1] else
  file.path(script_dir, "STAR_superpopulation_metrics")
output_dir <- if (length(args) >= 2L) args[2] else
  file.path(input_dir, "presentation")

# Inspect every scenario directory; generated output is excluded on reruns.
scenario_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
scenario_dirs <- scenario_dirs[normalizePath(scenario_dirs, mustWork = FALSE) !=
                               normalizePath(output_dir, mustWork = FALSE)]
files <- file.path(scenario_dirs, "LOSO_scenario_average_metrics.csv")
if (!length(files)) stop("No scenario directories found in ", input_dir)
if (any(!file.exists(files))) {
  stop("Missing scenario metric files: ",
       paste(files[!file.exists(files)], collapse = ", "))
}
required_columns <- c("scenario", "J", "N_j", "c_scale", "metric",
                      "mean_over_replications", "sd_across_replications",
                      "monte_carlo_se", "replications")
metrics <- dplyr::bind_rows(lapply(files, function(file) {
  x <- readr::read_csv(file, show_col_types = FALSE)
  if (!all(required_columns %in% names(x))) stop("Invalid columns in ", file)
  if (dplyr::n_distinct(x$scenario) != 1L ||
      any(x$scenario != basename(dirname(file)))) {
    stop("Scenario label disagrees with directory: ", file)
  }
  x
}))
levels <- c(50L, 80L, 90L, 95L)
required_metrics <- c("log_predictive_score", "RMSE",
                      paste0("coverage_", levels), paste0("mean_length_", levels))
if (anyDuplicated(metrics[c("scenario", "metric")])) {
  stop("Duplicate scenario-metric pairs.")
}
settings <- dplyr::distinct(metrics, scenario, J, N_j, c_scale, replications)
if (anyDuplicated(settings$scenario) || nrow(settings) != 9L) {
  stop("Expected nine scenarios with constant settings within each scenario.")
}
for (scenario_name in settings$scenario) {
  x <- dplyr::filter(metrics, scenario == scenario_name, metric %in% required_metrics)
  if (!setequal(x$metric, required_metrics) ||
      any(!is.finite(x$mean_over_replications)) ||
      any(!is.finite(x$monte_carlo_se) | x$monte_carlo_se < 0)) {
    stop("Missing or invalid presentation metrics for ", scenario_name)
  }
}
if (any(!is.finite(settings$J) | settings$J <= 0 |
        !is.finite(settings$N_j) | settings$N_j <= 0 |
        !is.finite(settings$c_scale) | settings$c_scale < 0 |
        !is.finite(settings$replications) | settings$replications < 2)) {
  stop("Invalid scenario settings or replication counts.")
}
if (any(grepl("^coverage_", metrics$metric) &
        (metrics$mean_over_replications < 0 | metrics$mean_over_replications > 1)) ||
    any(metrics$metric %in% c("RMSE", paste0("mean_length_", levels)) &
        metrics$mean_over_replications < 0)) stop("Invalid coverage or length/RMSE.")
baseline <- dplyr::filter(settings, scenario == "baseline")
if (nrow(baseline) != 1L) stop("Exactly one baseline scenario is required.")

# Select comparative statics by parameter values, not experiment labels:
# the single baseline observation belongs to each of the three panels.
panels <- dplyr::bind_rows(
  metrics |> dplyr::filter(N_j == baseline$N_j, c_scale == baseline$c_scale) |>
    dplyr::mutate(panel = "(a) Number of studies: J", x = J),
  metrics |> dplyr::filter(J == baseline$J, N_j == baseline$N_j) |>
    dplyr::mutate(panel = "(b) Heterogeneity: c", x = c_scale),
  metrics |> dplyr::filter(J == baseline$J, c_scale == baseline$c_scale) |>
    dplyr::mutate(panel = "(c) Within-study sample size: N[j]", x = N_j)
) |>
  dplyr::mutate(panel = factor(panel, levels = c("(a) Number of studies: J",
    "(b) Heterogeneity: c", "(c) Within-study sample size: N[j]")))
if (!setequal(settings$scenario, panels$scenario)) {
  stop("Some scenarios change more than one baseline parameter.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_figure <- function(metric_name, y_label, nominal = NULL) {
  x <- dplyr::filter(panels, metric == metric_name) |>
    dplyr::mutate(mc_lower = mean_over_replications - mc_multiplier * monte_carlo_se,
                  mc_upper = mean_over_replications + mc_multiplier * monte_carlo_se)
  if (!is.null(nominal)) {
    x <- x |> dplyr::mutate(mc_lower = pmax(0, mc_lower),
                            mc_upper = pmin(1, mc_upper))
  }
  y_limits <- if (!is.null(nominal)) c(0, 1) else if (show_mc_bars) {
    range(x$mc_lower, x$mc_upper)
  } else range(x$mean_over_replications)
  # Separate plots give each panel its own exact breaks while preserving
  # numeric spacing and a common y scale. Facet-wide pooled breaks would
  # incorrectly put J values (e.g. 40 and 80) on the N_j axis.
  plots <- lapply(levels(x$panel), function(panel_name) {
    panel_data <- dplyr::filter(x, panel == panel_name)
    p <- ggplot2::ggplot(panel_data,
                         ggplot2::aes(x = x, y = mean_over_replications)) +
    ggplot2::geom_line(linewidth = 0.5, colour = "#244A64") +
    (if (show_mc_bars) ggplot2::geom_linerange(
      ggplot2::aes(ymin = mc_lower, ymax = mc_upper),
      linewidth = 0.45, colour = "#244A64")) +
    ggplot2::geom_point(size = 1.15, colour = "#244A64") +
    ggplot2::geom_point(data = dplyr::filter(panel_data, scenario == "baseline"),
               shape = 21, size = 1.5, fill = "white", colour = "#244A64") +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_x_continuous(breaks = sort(unique(panel_data$x))) +
    ggplot2::scale_y_continuous(limits = y_limits,
      breaks = if (!is.null(nominal)) seq(0, 1, 0.2) else ggplot2::waiver()) +
    ggplot2::labs(x = NULL, y = if (panel_name == levels(x$panel)[1]) y_label else NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = "grey95"),
          plot.caption = ggplot2::element_text(hjust = 0))
    if (!is.null(nominal)) {
      p <- p + ggplot2::geom_hline(yintercept = nominal, linetype = "dashed",
                                  colour = "grey40")
    }
    ggplot2::ggplotGrob(p)
  })
  # Equalize grob columns so the three plotting areas have equal widths.
  widths <- do.call(grid::unit.pmax, lapply(plots, function(p) p$widths))
  for (i in seq_along(plots)) {
    plots[[i]]$widths <- widths
  }
  caption <- if (show_mc_bars) {
    "Means +/- 1.96 Monte Carlo SE (approximate 95% MC intervals). Open circles: baseline."
  } else "Means over simulation replications. Open circles: baseline."
  # ggplotGrob objects share the name "layout". Putting them directly in
  # grobTree can replace siblings with the same name during grid drawing.
  # Assign explicit cells and unique names instead of nested viewports.
  combined <- gtable::gtable(
    widths = grid::unit(rep(1, 3), "null"),
    heights = grid::unit(c(0.92, 0.08), "null"),
    name = "simulation_figure")
  combined <- gtable::gtable_add_grob(combined, grobs = plots,
    t = rep(1L, 3), l = seq_len(3), clip = "off",
    name = paste0("simulation_panel_", seq_len(3)))
  gtable::gtable_add_grob(combined, grid::textGrob(
    caption, x = 0.02, y = 0.5, just = "left",
    gp = grid::gpar(fontsize = 9)),
    t = 2L, l = 1L, r = 3L, clip = "off", name = "simulation_caption")
}
save_figure <- function(p, name) {
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".pdf")), p,
         width = 11, height = 3.8, units = "in")
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), p,
         width = 11, height = 3.8, units = "in", dpi = 300)
}
save_figure(make_figure("log_predictive_score", "Log predictive score"),
            "figure_1_log_predictive_score")
save_figure(make_figure("RMSE", "RMSE"), "figure_2_RMSE")
save_figure(make_figure("coverage_95", "95% predictive coverage", 0.95),
            "figure_3_coverage_95")
save_figure(make_figure("mean_length_95", "95% predictive interval length"),
            "figure_4_interval_length_95")
for (level in levels) {
  save_figure(make_figure(paste0("coverage_", level),
                         paste0(level, "% predictive coverage"), level / 100),
              paste0("appendix_coverage_", level))
  save_figure(make_figure(paste0("mean_length_", level),
                         paste0(level, "% predictive interval length")),
              paste0("appendix_interval_length_", level))
}

# Baseline first, then J, c, N_j, sorted numerically within each dimension.
scenario_order <- settings |>
  dplyr::mutate(dimension = dplyr::case_when(scenario == "baseline" ~ 0L,
                               J != baseline$J ~ 1L,
                               c_scale != baseline$c_scale ~ 2L,
                               TRUE ~ 3L)) |>
  dplyr::arrange(dimension, J, c_scale, N_j) |>
  dplyr::pull(scenario)
main_table <- metrics |>
  dplyr::filter(metric %in% c("log_predictive_score", "RMSE", "coverage_95",
                       "mean_length_95")) |>
  dplyr::select(scenario, metric, mean_over_replications) |>
  tidyr::pivot_wider(names_from = metric, values_from = mean_over_replications) |>
  dplyr::arrange(match(scenario, scenario_order)) |>
  dplyr::select(Scenario = scenario, `Log predictive score` = log_predictive_score,
         RMSE, `95% coverage` = coverage_95,
         `95% interval length` = mean_length_95)
readr::write_csv(main_table, file.path(output_dir, "table_main.csv"))
display_table <- main_table
display_table[-1] <- lapply(display_table[-1], function(x) sprintf("%.4f", x))
markdown <- c(paste(names(display_table), collapse = " | "),
              paste(rep("---", ncol(display_table)), collapse = " | "),
              apply(display_table, 1, paste, collapse = " | "))
writeLines(markdown, file.path(output_dir, "table_main.md"))
appendix_table <- metrics |>
  dplyr::filter(metric %in% c(paste0("coverage_", levels),
                       paste0("mean_length_", levels))) |>
  dplyr::select(scenario, metric, mean_over_replications) |>
  tidyr::pivot_wider(names_from = metric, values_from = mean_over_replications) |>
  dplyr::arrange(match(scenario, scenario_order)) |>
  dplyr::select(scenario, dplyr::all_of(as.vector(rbind(paste0("coverage_", levels),
                                        paste0("mean_length_", levels)))))
readr::write_csv(appendix_table, file.path(output_dir, "table_appendix_intervals.csv"))
readr::write_csv(settings |> dplyr::arrange(match(scenario, scenario_order)),
                 file.path(output_dir, "scenario_settings.csv"))
readr::write_csv(metrics |> dplyr::filter(metric %in% required_metrics) |>
                   dplyr::arrange(match(scenario, scenario_order), metric),
                 file.path(output_dir, "presentation_metrics_with_mcse.csv"))
writeLines(c(
  "Predictive target: the true superpopulation treatment effect in a held-out study.",
  "Log predictive scores use natural logarithms; higher is better on the same outcome scale.",
  "RMSE is the supplied mean of replication-level RMSEs, not a pooled RMSE.",
  "Coverage measures calibration; interval length measures sharpness.",
  "Intervals are the equal-tailed predictive intervals computed by the metrics script.",
  "Dashed coverage lines indicate nominal coverage; open circles indicate baseline.",
  "All plots share a y-axis scale across panels within a figure.",
  "Tables retain full numerical precision in CSV; Markdown displays four decimals.",
  paste0("Monte Carlo uncertainty bars enabled: ", show_mc_bars, "."),
  "Bars show mean +/- 1.96 MCSE: approximate pointwise 95% intervals assuming independent replications and a normal approximation for the replication mean.",
  "Coverage MC intervals are truncated to [0,1]; near-boundary normal approximations may be inaccurate.",
  "MC intervals concern uncertainty in simulation means, not predictive uncertainty or variation across replications.",
  "MCSE and replication counts are also exported in presentation_metrics_with_mcse.csv.",
  "Intervals for scenario differences require the covariance of paired replication results; overlap of individual bars is not a formal test of differences.",
  paste0("Baseline: J=", baseline$J, ", c=", baseline$c_scale,
         ", N_j=", baseline$N_j, ".")
), file.path(output_dir, "presentation_notes.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
print(main_table)
message("Presentation saved to: ", normalizePath(output_dir))
