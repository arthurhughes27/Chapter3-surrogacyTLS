# Shared configuration sourced by analysis scripts.
# Centralises the random seed and standard input/output paths so they are
# defined once rather than repeated in every script.

RISE_SEED <- 18042025
SIMULATION_SEED <- 23102024

processed_data_path <- fs::path("data")
output_path <- fs::path("output")

# Number of cores for the parallelised simulation replicates in
# analysis/simulation/ and analysis/supplementary/ (via pbmcapply::pbmclapply
# in R/simulation_helpers.R). Each forked replicate worker does real,
# memory-heavy work (generating up to 500 candidates and running hundreds
# of per-candidate tests), so this deliberately does NOT default to
# detectCores() - 1: on a many-core machine that forks far more workers
# than there is RAM to comfortably hold, which is what crashed a system
# previously. Capped at 4 by default; raise it explicitly (e.g. in a
# script, before sourcing this file is too late - set
# options(simulation.n_cores = N) then reference that) if your machine has
# the memory to spare.
#
# rise.screen()/rise.evaluate() themselves are always called with
# n.cores = 1 inside simulate_gamma_pvalues() (see R/simulation_helpers.R)
# specifically to avoid double parallelisation: the outer pbmclapply over
# replicates already uses all the parallelism we want, so the per-marker
# screening inside a single replicate must stay single-core.
SIMULATION_N_CORES <- getOption("simulation.n_cores", max(1, min(parallel::detectCores() - 1, 4)))
