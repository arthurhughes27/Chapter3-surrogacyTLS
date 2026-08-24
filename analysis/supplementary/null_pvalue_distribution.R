# Supplementary. Data generation process 1 (simple): distribution of raw
# p-values under the null hypothesis. Sample size n = 50, predictors
# generated without correlation, across 1000 simulations.

library(ggplot2)
library(dplyr)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "simulation_helpers.R"))
source(fs::path("R", "simulation_plots.R"))

set.seed(SIMULATION_SEED)

results_path <- fs::path(output_path, "results", "supplementary", "null_pvalue_distribution.rds")
figure_path <- fs::path(output_path, "figures", "supplementary", "null_pvalue_distribution")

alpha <- 0.05
y1_mean <- 3; y1_sd <- 1
y0_mean <- 0; y0_sd <- 1
n <- 50
p <- 500
prop_valid <- 0
n_sim <- 1000
corr <- 0

null_pvalues <- cache_rds(results_path, {
  res <- simulate_screening_metrics(
    n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = n_sim,
    valid_sigma = NA, corr = corr,
    y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
    mode = "simple", alpha = alpha, n_cores = SIMULATION_N_CORES
  )
  # Pool p-values across every candidate of every replicate: since
  # prop_valid = 0, all p candidates in all n_sim replicates are null, so
  # this reproduces the original figure's pooled null p-value distribution
  # (n_sim x p values in total) computed in a single parallelised call
  # rather than n_sim separate single-replicate calls.
  data.frame(p = as.vector(res[["p_values"]][["p_unadjusted"]]))
})

plot_null_pvalue_histogram(null_pvalues, out_path = figure_path)
