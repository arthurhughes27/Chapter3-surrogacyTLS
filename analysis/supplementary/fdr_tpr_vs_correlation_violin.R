# Supplementary. Data generation process 1 (simple), scenario 2: violin
# plots of empirical power and false discovery proportion prior to
# multiple testing correction, for a fixed sample size n = 50 and average
# surrogate strength U_S bar ~ 0.9, across different levels of
# inter-predictor correlation.

library(ggplot2)
library(dplyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "supplementary", "fdr_tpr_vs_correlation_violin.rds")
figure_path <- fs::path(output_path, "figures", "supplementary", "fdr_tpr_vs_correlation_violin")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
valid_sigma <- 1.8 # corresponds to avg surrogate strength ~ 0.9
n <- 50
p <- 500
prop_valid <- 0.1
n_sim <- 5

corr_grid <- seq(0, 0.5, 0.1)

metrics_fdr_tpr_vs_corr <- do.call(rbind, checkpoint_grid(
  path = results_path,
  grid = corr_grid,
  key_fn = function(corr) paste0("corr=", corr),
  label_fn = function(corr) sprintf("correlation = %.1f", corr),
  compute_one = function(corr) {
    res <- simulate_screening_metrics(
      n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = n_sim,
      valid_sigma = valid_sigma, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", alpha = alpha, n_cores = SIMULATION_N_CORES
    )
    df <- res[["per_replicate"]] %>% filter(correction == "unadjusted")
    data.frame(correlation = corr, fdr = df$fdr, tpr = df$tpr)
  }
))

plot_power_fdr_violin_by_corr(metrics_fdr_tpr_vs_corr, out_path = figure_path)

rm(list = ls())
