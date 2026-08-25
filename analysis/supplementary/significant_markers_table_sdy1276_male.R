# Supplementary. Full table of significant screening-stage markers for
# the SDY1276 (TIV) application restricted to sex == "Male"
# (analysis/supplementary/RISE-SDY1276-Male.R). Marker (capitalised),
# delta with a 90% CI, adjusted p-value, and normalised screening
# weight, sorted by adjusted p-value.
#
# Requires analysis/supplementary/RISE-SDY1276-Male.R to have been run
# first, since it depends on that script's saved screening_sdy1276_male.rds.

library(dplyr)
library(xtable)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "application_tables.R"))

results_path <- fs::path(output_path, "results", "supplementary")
table_path <- fs::path(output_path, "tables", "supplementary", "significant_markers_sdy1276_male.tex")

screening <- readRDS(fs::path(results_path, "screening_sdy1276_male.rds"))

build_significant_markers_table(
  screening_metrics = screening$screening.metrics,
  screening_weights = screening$screening.weights,
  significant_markers = screening$significant.markers,
  out_path = table_path,
  caption = "Significant screening-stage markers for the SDY1276 (TIV) surrogate signature, restricted to male participants.",
  label = "tab:significant-markers-sdy1276-male"
)
