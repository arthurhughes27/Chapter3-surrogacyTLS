# Shared configuration sourced by analysis scripts.
# Centralises the random seed and standard input/output paths so they are
# defined once rather than repeated in every script.

RISE_SEED <- 18042025
SIMULATION_SEED <- 23102024

processed_data_path <- fs::path("data")
output_path <- fs::path("output")

SIMULATION_N_CORES <- getOption("simulation.n_cores", max(1, min(parallel::detectCores() - 1, 9)))
