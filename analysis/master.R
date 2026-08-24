# Master script to run the full pipeline end-to-end: preprocessing of the
# raw study data, the RISE application analyses that depend on it, and the
# (independent) RISE simulation studies.

source(fs::path("analysis", "preprocessing", "preprocessing_master.R"))
source(fs::path("analysis", "application", "application_master.R"))
source(fs::path("analysis", "simulation", "simulation_master.R"))
source(fs::path("analysis", "supplementary", "supplementary_master.R"))
