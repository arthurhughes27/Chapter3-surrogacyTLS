# Master script to run all supplementary material in order:
#   simulation/            - supplementary simulation studies
#   supplementary_applications/ - application variants not in the main text
#   sensitivity_analyses/  - fixed-epsilon screening sensitivity tables
#   full_tables/           - full significant-markers tables per application
#     (require application_master.R, and supplementary_applications/, to
#     have already run, since they read those scripts' saved screening_*.rds)
#
# Figures/results/tables from every subfolder still all save into the
# same flat output/{figures,results,tables}/supplementary/ locations.

source(fs::path("analysis", "supplementary", "simulation", "fdr_tpr_vs_correlation_violin.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "power_fdr_vs_strength_by_correction.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "fpr_vs_sample_size_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "fpr_vs_correlation_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "power_fdr_vs_strength_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "gamma_pvalue_vs_prop_invalid_complex_dgp.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "null_pvalue_distribution.R"))

gc()

source(fs::path("analysis", "supplementary", "simulation", "gamma_pvalue_vs_prop_invalid_by_p.R"))

gc()

source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-SDY1276-Male.R"))

gc()

source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-Ad26MVA-Reversed.R"))
source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-rVSV-Reversed.R"))

gc()

source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-SDY1276-Florida.R"))
source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-SDY1276-Brisbane10.R"))
source(fs::path("analysis", "supplementary", "supplementary_applications", "RISE-SDY1276-Brisbane59.R"))

gc()

source(fs::path("analysis", "supplementary", "sensitivity_analyses", "epsilon_sensitivity_tiv_female.R"))
source(fs::path("analysis", "supplementary", "sensitivity_analyses", "epsilon_sensitivity_tiv_male.R"))
source(fs::path("analysis", "supplementary", "sensitivity_analyses", "epsilon_sensitivity_ad26mva.R"))
source(fs::path("analysis", "supplementary", "sensitivity_analyses", "epsilon_sensitivity_rvsv.R"))

gc()

source(fs::path("analysis", "supplementary", "full_tables", "significant_markers_table_ad26mva.R"))
source(fs::path("analysis", "supplementary", "full_tables", "significant_markers_table_rvsv.R"))
source(fs::path("analysis", "supplementary", "full_tables", "significant_markers_table_sdy1276.R"))
source(fs::path("analysis", "supplementary", "full_tables", "significant_markers_table_sdy1276_male.R"))
