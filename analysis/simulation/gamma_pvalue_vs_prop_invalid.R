# Data generation process 1 (simple): distribution of the evaluation-stage
# p-value for gamma_S, as a function of the proportion of invalid
# surrogates making up gamma_S (a combination of 20 candidates). Sample
# size n = 50, valid surrogate strength U_Sj ~ 0.9 (calibrated via
# valid_sigma below). The nominal significance level alpha = 0.05 is
# plotted as a dashed reference line. Desired power for the composite
# surrogate is fixed at 80%.

library(ggplot2)
library(dplyr)
library(SurrogateRank)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "simulation", "gamma_pvalue_vs_prop_invalid.rds")
figure_path <- fs::path(output_path, "figures", "simulation", "gamma_pvalue_vs_prop_invalid")

y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
p <- 20
corr <- 0
n1 <- 25; n0 <- 25
valid_sigma <- 1.8 # corresponds to avg U_Y ~ 0.9 for this design
prop_invalid_grid <- seq(0, 1, 0.1)
n_sim <- 5

p_evaluate <- do.call(rbind, checkpoint_grid(
  path = results_path,
  grid = prop_invalid_grid,
  key_fn = function(prop_invalid) paste0("prop_invalid=", prop_invalid),
  label_fn = function(prop_invalid) sprintf("prop. invalid = %.1f", prop_invalid),
  compute_one = function(prop_invalid) {
    p_values <- simulate_gamma_pvalues(
      n1 = n1, n0 = n0, p = p, prop_invalid = prop_invalid,
      valid_sigma = valid_sigma, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", n_sim = n_sim, n_cores = SIMULATION_N_CORES
    )
    data.frame(prop_invalid = prop_invalid, p_value = p_values)
  }
))

plot_gamma_pvalue_boxplot(p_evaluate, out_path = figure_path)
