# Master script to run all supplementary material in order: simulation
# studies, then application significant-markers tables (which require
# analysis/application/application_master.R to have already run, since
# they read its saved screening_*.rds files).

source(fs::path("analysis", "supplementary", "fdr_tpr_vs_correlation_violin.R"))

gc()

source(fs::path("analysis", "supplementary", "power_fdr_vs_strength_by_correction.R"))

gc()

source(fs::path("analysis", "supplementary", "fpr_vs_sample_size_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "fpr_vs_correlation_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "power_fdr_vs_strength_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "gamma_pvalue_vs_prop_invalid_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "null_pvalue_distribution.R"))

gc()

source(fs::path("analysis", "supplementary", "gamma_pvalue_vs_prop_invalid_by_p.R"))

gc()

source(fs::path("analysis", "supplementary", "RISE-SDY1276-Male.R"))

gc()

source(fs::path("analysis", "supplementary", "significant_markers_table_ad26mva.R"))
source(fs::path("analysis", "supplementary", "significant_markers_table_rvsv.R"))
source(fs::path("analysis", "supplementary", "significant_markers_table_sdy1276.R"))

