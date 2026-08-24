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

# This figure has no natural outer grid (a single n_sim = 1000 call), so
# it's split into chunks of replicates purely for checkpointing: each
# chunk is saved as soon as it completes, rather than losing all 1000
# replicates to a crash near the end.
n_chunks <- 10
chunk_size <- n_sim / n_chunks

null_pvalues <- do.call(rbind, checkpoint_grid(
  path = results_path,
  grid = seq_len(n_chunks),
  key_fn = function(chunk) paste0("chunk=", chunk),
  label_fn = function(chunk) sprintf("replicates %d-%d of %d",
                                      (chunk - 1) * chunk_size + 1, chunk * chunk_size, n_sim),
  compute_one = function(chunk) {
    res <- simulate_screening_metrics(
      n1 = n / 2, n0 = n / 2, p = p, prop_valid = prop_valid, n_sim = chunk_size,
      valid_sigma = NA, corr = corr,
      y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
      mode = "simple", alpha = alpha, n_cores = SIMULATION_N_CORES
    )
    # Pool p-values across every candidate of every replicate in this
    # chunk: since prop_valid = 0, all p candidates in all replicates are
    # null, so this reproduces the original figure's pooled null p-value
    # distribution (n_sim x p values in total).
    data.frame(p = as.vector(res[["p_values"]][["p_unadjusted"]]))
  }
))

plot_null_pvalue_histogram(null_pvalues, out_path = figure_path)
