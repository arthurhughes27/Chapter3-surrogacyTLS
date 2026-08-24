# Data generation process 1 (simple), scenario 1: boxplots of observed
# false positive rates against different sample sizes in the
# uncorrelated, null-only setting (prop_valid = 0). The nominal
# significance level alpha = 0.05 is plotted as a dashed reference line.

library(ggplot2)
library(dplyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "simulation", "fpr_vs_sample_size.rds")
figure_path <- fs::path(output_path, "figures", "simulation", "fpr_vs_sample_size")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
p <- 500
prop_valid <- 0
n_sim <- 500
corr <- 0

n_grid <- seq(10, 100, 10)

metrics_fpr_vs_n <- cache_rds(results_path, {
  do.call(rbind, lapply(n_grid, function(n) {
    res <- simulate_screening_metrics(
      n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = n_sim,
      valid_sigma = NA, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", alpha = alpha, n_cores = SIMULATION_N_CORES
    )
    df <- res[["per_replicate"]] %>% filter(correction == "unadjusted")
    data.frame(n = n, fpr = df$fpr)
  }))
})

plot_fpr_boxplot(metrics_fpr_vs_n, x_var = "n", x_lab = "Sample size", out_path = figure_path)
