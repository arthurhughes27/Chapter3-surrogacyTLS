# Supplementary. Full table of significant screening-stage markers for
# the SDY1276 (TIV) application (analysis/application/RISE-SDY1276.R).
# Marker (capitalised), delta with a 90% CI, adjusted p-value, and
# normalised screening weight, sorted by adjusted p-value.
#
# Requires analysis/application/RISE-SDY1276.R to have been run first,
# since it depends on that script's saved screening_sdy1276.rds.

library(dplyr)
library(xtable)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "application_tables.R"))

results_path <- fs::path(output_path, "results", "application")
table_path <- fs::path(output_path, "tables", "supplementary", "significant_markers_sdy1276.tex")

screening <- readRDS(fs::path(results_path, "screening_sdy1276.rds"))

build_significant_markers_table(
  screening_metrics = screening$screening.metrics,
  screening_weights = screening$screening.weights,
  significant_markers = screening$significant.markers,
  out_path = table_path,
  caption = "Significant screening-stage markers for the SDY1276 (TIV) surrogate signature.",
  label = "tab:significant-markers-sdy1276"
)
