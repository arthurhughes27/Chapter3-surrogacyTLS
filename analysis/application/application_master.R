# Master script to run all application analyses in order

source(fs::path("analysis", "application", "RISE-SDY1276.R"))
source(fs::path("analysis", "application", "RISE-Ad26MVA.R"))
source(fs::path("analysis", "application", "RISE-rVSV.R"))

# Cross-study comparison of the resulting surrogate marker signatures
source(fs::path("analysis", "application", "RISE_signatures.R"))
