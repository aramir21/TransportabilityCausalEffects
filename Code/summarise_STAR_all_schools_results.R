# Audit the completed 79-fold comparison and summarise saved HMC diagnostics.
output_dir <- "star_output/model_comparison"
s <- read.csv(file.path(output_dir, "comparison_summary.csv"))
d <- read.csv(file.path(output_dir, "comparison_school_metrics.csv"))
hmc <- read.csv(file.path(output_dir, "comparison_diagnostics.csv"))
stopifnot(nrow(d) == 316L, nrow(hmc) == 316L, nrow(s) == 4L,
          !anyDuplicated(d[c("model", "fold")]),
          all(is.finite(d$log_predictive_score)), all(s$schools == 79L))
for (model in c("M0", "M1", "M2", "M3")) {
  school <- d[d$model == model, ]
  row <- s[s$model == model, ]
  stopifnot(identical(sort(school$fold), seq_len(79L)),
            isTRUE(all.equal(row$mean_log_score, mean(school$log_predictive_score))),
            isTRUE(all.equal(row$RMSE_reference_means, sqrt(mean(school$squared_mean_error)))),
            isTRUE(all.equal(row$posterior_expected_RMSE,
                             sqrt(mean(school$posterior_expected_squared_error)))))
  for (level in c(50, 80, 90, 95)) {
    stopifnot(isTRUE(all.equal(row[[paste0("posterior_coverage_", level)]],
                              mean(school[[paste0("posterior_coverage_", level)]]))),
              isTRUE(all.equal(row[[paste0("mean_interval_length_", level)]],
                                mean(school[[paste0("length_", level)]]))))
  }
}
diagnostics_summary <- do.call(rbind, lapply(c("M0", "M1", "M2", "M3"), function(model) {
  x <- hmc[hmc$model == model, ]
  stopifnot(nrow(x) == 79L, all(is.finite(as.matrix(x[c("max_rhat", "min_n_eff",
    "divergences", "treedepth_hits", "min_ebfmi")]))))
  data.frame(model = model, folds = nrow(x), max_rhat = max(x$max_rhat),
    min_n_eff = min(x$min_n_eff), divergences = sum(x$divergences),
    treedepth_hits = sum(x$treedepth_hits), min_ebfmi = min(x$min_ebfmi),
    folds_rhat_above_1_01 = sum(x$max_rhat > 1.01),
    folds_ess_below_400 = sum(x$min_n_eff < 400),
    folds_ebfmi_below_0_3 = sum(x$min_ebfmi < .3))
}))
write.csv(diagnostics_summary, file.path(output_dir, "comparison_diagnostics_summary.csv"),
          row.names = FALSE)
print(s[c("model", "schools", "mean_log_score", "RMSE_reference_means", "posterior_expected_RMSE")],
      row.names = FALSE)
print(read.csv(file.path(output_dir, "comparison_coverage.csv")), row.names = FALSE)
print(diagnostics_summary, row.names = FALSE)
cat("All 316 model-school rows and reported averages verified.\n")
