# Supplementary. Full table of significant screening-stage markers for
# the Ad26MVA application (analysis/application/RISE-Ad26MVA.R): the
# PREVAC-screened, EBOVAC2-evaluated Ad26MVA/placebo surrogate signature.
# Marker (capitalised), delta with a 90% CI, adjusted p-value, and
# normalised screening weight, sorted by adjusted p-value.
#
# Requires analysis/application/RISE-Ad26MVA.R to have been run first,
# since it depends on that script's saved screening_ad26mva.rds.

library(dplyr)
library(xtable)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "application_tables.R"))

results_path <- fs::path(output_path, "results", "application")
table_path <- fs::path(output_path, "tables", "supplementary", "significant_markers_ad26mva.tex")

screening <- readRDS(fs::path(results_path, "screening_ad26mva.rds"))

build_significant_markers_table(
  screening_metrics = screening$screening.metrics,
  screening_weights = screening$screening.weights,
  significant_markers = screening$significant.markers,
  out_path = table_path,
  caption = "Significant screening-stage markers for the Ad26MVA/placebo surrogate signature (screened in PREVAC, evaluated in EBOVAC2).",
  label = "tab:significant-markers-ad26mva"
)
