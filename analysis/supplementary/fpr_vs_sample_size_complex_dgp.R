# Supplementary. Data generation process 2 (complex), scenario 1: boxplots
# of observed false positive rates against different sample sizes in the
# uncorrelated, null-only setting (prop_valid = 0). The nominal
# significance level alpha = 0.05 is plotted as a dashed reference line.

library(ggplot2)
library(dplyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "supplementary", "fpr_vs_sample_size_complex_dgp.rds")
figure_path <- fs::path(output_path, "figures", "supplementary", "fpr_vs_sample_size_complex_dgp")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
p <- 500
prop_valid <- 0
n_sim <- 5
corr <- 0

n_grid <- seq(10, 100, 10)

metrics_fpr_vs_n_complex <- do.call(rbind, checkpoint_grid(
  path = results_path,
  grid = n_grid,
  key_fn = function(n) paste0("n=", n),
  label_fn = function(n) sprintf("sample size n = %d", n),
  compute_one = function(n) {
    res <- simulate_screening_metrics(
      n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = n_sim,
      valid_sigma = NA, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "complex", alpha = alpha, n_cores = SIMULATION_N_CORES
    )
    df <- res[["per_replicate"]] %>% filter(correction == "unadjusted")
    data.frame(n = n, fpr = df$fpr)
  }
))

plot_fpr_boxplot(metrics_fpr_vs_n_complex, x_var = "n", x_lab = "Sample size", out_path = figure_path)

rm(list = ls())
