# =============================================================================
# Effective sample sizes by analysis plan
#
# For each planned RISE analysis, an "eligible" participant is one with
# every piece of data the analysis plan requires (no missing timepoints):
#
#   SDY1276 TIV:                 GE at day 0 & day 1;
#                                 antibody at day 0 & day 28 (Female only)
#   PREVAC Ad26/MVA + placebo:   GE at day 7;
#                                 antibody at day 365
#   PREVAC rVSV + placebo:       GE at day 7;
#                                 antibody at day 180
#   EBOVAC2 Ad26/MVA:            GE at day 0 & day 7;
#                                 antibody at day 0 & day 365
#   Hamburg rVSV:                GE at day 0 & day 7;
#                                 antibody at day 0 & day 180
#
# PREVAC Ad26/MVA + placebo and PREVAC rVSV + placebo pool two
# study-vaccine arms into a single treated-vs-untreated comparison, so
# the reported N is the pooled total across both arms.
# =============================================================================

# ---- Libraries ----
library(tidyverse)
library(fs)
library(xtable)

# ---- Paths ----
processed_data_path       <- fs::path("data")
descriptive_tables_folder <- fs::path("output", "tables", "descriptive")
fs::dir_create(descriptive_tables_folder)

# ---- Load harmonised clinical data ----
df_clinical_all <- readRDS(fs::path(processed_data_path, "df_clinical_all.rds"))

# ---- Helpers ----

# Participants with a transcriptomic (GE) sample at a given timepoint
ge_participants <- function(t) {
  df_clinical_all %>%
    filter(time == t) %>%
    distinct(participant_id) %>%
    pull(participant_id)
}

# One row per participant, carrying study_vaccine, sex, and the antibody
# columns needed below (values are repeated across a participant's rows
# in df_clinical_all, so distinct() collapses them safely)
ab_df <- df_clinical_all %>%
  distinct(participant_id, study_vaccine, sex, ab_p_0, ab_p_28, ab_p_180, ab_p_365)

# Restrict to the given study-vaccine group(s) (and, optionally, sex),
# then require every listed GE timepoint and antibody column to be
# non-missing
eligible_participants <- function(groups, ge_times = character(0), ab_cols = character(0),
                                  sex_filter = NULL) {
  d <- ab_df %>% filter(study_vaccine %in% groups)
  if (!is.null(sex_filter)) {
    d <- d %>% filter(sex == sex_filter)
  }
  for (t in ge_times) {
    d <- d %>% filter(participant_id %in% ge_participants(t))
  }
  for (col in ab_cols) {
    d <- d %>% filter(!is.na(.data[[col]]))
  }
  d
}

# ---- Analysis plans ----

analysis_plans <- list(
  list(
    analysis   = "SDY1276 TIV",
    groups     = "SDY1276-Influenza (IN)",
    ge_times   = c("P+0D", "P+1D"),
    ab_cols    = c("ab_p_0", "ab_p_28"),
    sex_filter = "Female"
  ),
  list(
    analysis   = "PREVAC Ad26/MVA + placebo",
    groups     = c("prevac-Ad26MVA", "prevac-placebo"),
    ge_times   = "P+7D",
    ab_cols    = "ab_p_365",
    sex_filter = NULL
  ),
  list(
    analysis   = "PREVAC rVSV + placebo",
    groups     = c("prevac-rVSV", "prevac-placebo"),
    ge_times   = "P+7D",
    ab_cols    = "ab_p_180",
    sex_filter = NULL
  ),
  list(
    analysis   = "EBOVAC2 Ad26/MVA",
    groups     = "ebovac2-Ad26MVA",
    ge_times   = c("P+0D", "P+7D"),
    ab_cols    = c("ab_p_0", "ab_p_365"),
    sex_filter = NULL
  ),
  list(
    analysis   = "Hamburg rVSV",
    groups     = "hamburg-rVSV",
    ge_times   = c("P+0D", "P+7D"),
    ab_cols    = c("ab_p_0", "ab_p_180"),
    sex_filter = NULL
  )
)

# ---- Build the table: one row per analysis, N = total eligible participants ----

build_plan_row <- function(plan) {
  eligible <- eligible_participants(plan$groups, plan$ge_times, plan$ab_cols, plan$sex_filter)
  tibble(Analysis = plan$analysis, N = nrow(eligible))
}

sample_size_table <- map_dfr(analysis_plans, build_plan_row)

sample_size_table

# ---- Save as a LaTeX table ----

print(
  xtable::xtable(
    sample_size_table,
    caption = "Effective sample sizes by analysis plan: number of eligible participants with complete data for the timepoints required by each planned analysis.",
    label   = "tab:effective_sample_sizes"
  ),
  include.rownames = FALSE,
  booktabs = TRUE,
  file = fs::path(descriptive_tables_folder, "effective_sample_sizes.tex")
)

rm(list = ls())
