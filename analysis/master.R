# Master script to run the full pipeline end-to-end: preprocessing of the
# raw study data, followed by the RISE application analyses that depend
# on it.

source(fs::path("analysis", "preprocessing", "preprocessing_master.R"))
source(fs::path("analysis", "application", "application_master.R"))
