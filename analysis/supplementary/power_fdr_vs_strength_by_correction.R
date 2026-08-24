# Supplementary. Data generation process 1 (simple), scenario 2: empirical
# power and false discovery proportion prior to and after multiple
# testing correction (Benjamini-Hochberg, Bonferroni, Benjamini-Yekutieli,
# unadjusted), as a function of average surrogate strength, for a fixed
# sample size n = 50.

library(ggplot2)
library(dplyr)
library(tidyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "supplementary", "power_fdr_vs_strength_by_correction.rds")
figure_path <- fs::path(output_path, "figures", "supplementary", "power_fdr_vs_strength_by_correction")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
n <- 50
p <- 500
prop_valid <- 0.1
n_sim <- 500
corr <- 0

sigma_grid <- c(0.01, 0.65, 1.8, 3, 5.5, 9, 15, 30, 68, 244)

metrics_by_correction <- cache_rds(results_path, {
  do.call(rbind, lapply(sigma_grid, function(sigma) {
    avg_us <- round(mean(calc_truth(
      p = p, prop_valid = prop_valid, valid_sigma = sigma, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple"
    )$us_true), 2)

    res <- simulate_screening_metrics(
      n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = n_sim,
      valid_sigma = sigma, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", alpha = alpha, n_cores = SIMULATION_N_CORES
    )

    data.frame(
      avg_us = avg_us,
      avg_fdr_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_fdr"]],
      avg_tpr_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_tpr"]],
      avg_fdr_bonf = res[["metrics"]][["metrics_bonf"]][["avg_fdr"]],
      avg_tpr_bonf = res[["metrics"]][["metrics_bonf"]][["avg_tpr"]],
      avg_fdr_bh = res[["metrics"]][["metrics_bh"]][["avg_fdr"]],
      avg_tpr_bh = res[["metrics"]][["metrics_bh"]][["avg_tpr"]],
      avg_fdr_by = res[["metrics"]][["metrics_by"]][["avg_fdr"]],
      avg_tpr_by = res[["metrics"]][["metrics_by"]][["avg_tpr"]]
    )
  }))
})

plot_power_fdr_by_correction(metrics_by_correction, out_path = figure_path)
