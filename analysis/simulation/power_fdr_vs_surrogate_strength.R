# Data generation process 1 (simple), scenario 2: empirical power and
# false discovery proportion prior to multiple testing correction, as a
# function of average surrogate strength (U_S bar), for three different
# sample sizes.

library(ggplot2)
library(dplyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "simulation", "power_fdr_vs_surrogate_strength.rds")
figure_path <- fs::path(output_path, "figures", "simulation", "power_fdr_vs_surrogate_strength")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
p <- 500
prop_valid <- 0.1
n_sim <- 500
corr <- 0

# valid_sigma values chosen (by calc_truth(), see analysis/simulation/README
# or the original thesis text) to give average surrogate strengths spread
# roughly evenly across (0.5, 1).
sigma_grid <- c(0.01, 0.65, 1.8, 3, 5.5, 9, 15, 30, 68, 244)
n_grid <- c(30, 50, 100)

# Flatten the (sigma, n) grid to one list of pairs so each combination is
# checkpointed independently.
combos <- expand.grid(sigma = sigma_grid, n = n_grid, KEEP.OUT.ATTRS = FALSE)
combo_grid <- split(combos, seq_len(nrow(combos)))

metrics_strength <- do.call(rbind, checkpoint_grid(
  path = results_path,
  grid = combo_grid,
  key_fn = function(combo) sprintf("sigma=%s_n=%s", combo$sigma, combo$n),
  label_fn = function(combo) sprintf("sigma = %s, n = %s", combo$sigma, combo$n),
  compute_one = function(combo) {
    sigma <- combo$sigma
    n <- combo$n

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
    data.frame(n = n, sigma_S = sigma, avg_us = avg_us,
               avg_fdr = res[["metrics"]][["metrics_unadjusted"]][["avg_fdr"]],
               avg_tpr = res[["metrics"]][["metrics_unadjusted"]][["avg_tpr"]])
  }
))

plot_power_fdr_vs_strength(metrics_strength, out_path = figure_path)
