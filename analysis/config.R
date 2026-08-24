# Shared configuration sourced by analysis scripts.
# Centralises the random seed and standard input/output paths so they are
# defined once rather than repeated in every script.

RISE_SEED <- 18042025
SIMULATION_SEED <- 23102024

processed_data_path <- fs::path("data")
output_path <- fs::path("output")

# Number of cores available for the parallelised simulation replicates in
# analysis/simulation/ and analysis/supplementary/.
SIMULATION_N_CORES <- max(1, parallel::detectCores() - 1)
