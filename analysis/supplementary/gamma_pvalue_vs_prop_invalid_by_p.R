# Supplementary. Data generation process 1 (simple): distribution of the
# evaluation-stage p-value for gamma_S as a function of the proportion of
# invalid surrogates making up gamma_S, for two different numbers of
# candidates: a) 100 candidates, b) 10 candidates. Sample size n = 50,
# valid surrogate strength U_Sj ~ 0.9. The nominal significance level
# alpha = 0.05 is plotted as a dashed reference line. Desired power for
# the composite surrogate is fixed at 80%.
#
# NOTE: in the original simulation_combined.Rmd, panel b) (p = 10) only
# simulated grid points 10-11 of 11 (`for (i in 10:length(prop_invalid_grid))`
# instead of `for (i in 1:length(prop_invalid_grid))`), leaving most of
# that panel's x-axis unsimulated. That has been fixed here so the full
# grid is simulated for both panels.

library(ggplot2)
library(dplyr)
library(SurrogateRank)
library(Rsurrogate)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path_a <- fs::path(output_path, "results", "supplementary", "gamma_pvalue_vs_prop_invalid_p100.rds")
results_path_b <- fs::path(output_path, "results", "supplementary", "gamma_pvalue_vs_prop_invalid_p10.rds")
figure_path <- fs::path(output_path, "figures", "supplementary", "gamma_pvalue_vs_prop_invalid_by_p")

y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
corr <- 0
n1 <- 25; n0 <- 25
valid_sigma <- 1.8 # corresponds to avg U_Y ~ 0.9 for this design
prop_invalid_grid <- seq(0, 1, 0.1)
n_sim <- 500

run_gamma_grid <- function(p) {
  do.call(rbind, lapply(prop_invalid_grid, function(prop_invalid) {
    p_values <- simulate_gamma_pvalues(
      n1 = n1, n0 = n0, p = p, prop_invalid = prop_invalid,
      valid_sigma = valid_sigma, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", n_sim = n_sim, n_cores = SIMULATION_N_CORES
    )
    data.frame(prop_invalid = prop_invalid, p_value = p_values)
  }))
}

p_evaluate_p100 <- cache_rds(results_path_a, run_gamma_grid(p = 100))
p_evaluate_p10 <- cache_rds(results_path_b, run_gamma_grid(p = 10))

plot_gamma_pvalue_boxplot_pair(p_evaluate_p100, p_evaluate_p10, out_path = figure_path)
