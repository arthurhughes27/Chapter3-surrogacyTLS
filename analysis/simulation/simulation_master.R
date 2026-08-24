# Master script to run all main-text simulation studies in order.

source(fs::path("analysis", "simulation", "fpr_vs_sample_size.R"))

gc()

source(fs::path("analysis", "simulation", "fpr_vs_correlation.R"))

gc()
source(fs::path("analysis", "simulation", "power_fdr_vs_surrogate_strength.R"))

gc()

source(fs::path("analysis", "simulation", "gamma_pvalue_vs_prop_invalid.R"))

gc()
