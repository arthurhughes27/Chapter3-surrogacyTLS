# Master script to run all supplementary simulation studies in order.

source(fs::path("analysis", "supplementary", "fdr_tpr_vs_correlation_violin.R"))
source(fs::path("analysis", "supplementary", "power_fdr_vs_strength_by_correction.R"))
source(fs::path("analysis", "supplementary", "fpr_vs_sample_size_complex_dgp.R"))
source(fs::path("analysis", "supplementary", "fpr_vs_correlation_complex_dgp.R"))
source(fs::path("analysis", "supplementary", "power_fdr_vs_strength_complex_dgp.R"))
source(fs::path("analysis", "supplementary", "gamma_pvalue_vs_prop_invalid_complex_dgp.R"))
source(fs::path("analysis", "supplementary", "null_pvalue_distribution.R"))
source(fs::path("analysis", "supplementary", "gamma_pvalue_vs_prop_invalid_by_p.R"))
