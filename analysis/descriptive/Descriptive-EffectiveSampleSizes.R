# =============================================================================
# Effective sample sizes by analysis plan
#
# For each planned RISE analysis, an "eligible" participant is one with
# every piece of data the analysis plan requires (no missing timepoints):
#
#   SDY1276 TIV (paired):                 GE at day 0 & day 1;
#                                          antibody at day 0 & day 28
#   PREVAC Ad26/MVA + placebo (pooled):   GE at day 7;
#                                          antibody at day 365
#   PREVAC rVSV + placebo (pooled):       GE at day 7;
#                                          antibody at day 180
#   EBOVAC2 Ad26/MVA (paired):            GE at day 0 & day 7;
#                                          antibody at day 0 & day 365
#   Hamburg rVSV (paired):                GE at day 0 & day 7;
#                                          antibody at day 0 & day 180
#
# "Pooled" analyses combine two study-vaccine arms into a single
# treated-vs-untreated comparison, so both the per-arm and the pooled
# total are reported. "Paired" analyses are single-arm (baseline vs.
# follow-up within the same participants), so only one count applies.
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

# One row per participant, carrying study_vaccine and the antibody
# columns needed below (values are repeated across a participant's rows
# in df_clinical_all, so distinct() collapses them safely)
ab_df <- df_clinical_all %>%
  distinct(participant_id, study_vaccine, ab_p_0, ab_p_28, ab_p_180, ab_p_365)

# Restrict to the given study-vaccine group(s), then require every listed
# GE timepoint and antibody column to be non-missing
eligible_participants <- function(groups, ge_times = character(0), ab_cols = character(0)) {
  d <- ab_df %>% filter(study_vaccine %in% groups)
  for (t in ge_times) {
    d <- d %>% filter(participant_id %in% ge_participants(t))
  }
  for (col in ab_cols) {
    d <- d %>% filter(!is.na(.data[[col]]))
  }
  d
}

# y-axis-style display labels, matching Descriptive-Ebolavirus.R
group_display_labels <- c(
  "prevac-rVSV"            = "PREVAC rVSV",
  "prevac-Ad26MVA"         = "PREVAC Ad26/MVA",
  "prevac-placebo"         = "PREVAC placebo",
  "ebovac2-Ad26MVA"        = "EBOVAC2 Ad26/MVA",
  "hamburg-rVSV"           = "Hamburg rVSV",
  "SDY1276-Influenza (IN)" = "SDY1276 TIV"
)

# ---- Analysis plans ----

analysis_plans <- list(
  list(
    analysis = "SDY1276 TIV (paired)",
    groups   = "SDY1276-Influenza (IN)",
    ge_times = c("P+0D", "P+1D"),
    ab_cols  = c("ab_p_0", "ab_p_28"),
    pooled   = FALSE
  ),
  list(
    analysis = "PREVAC Ad26/MVA + placebo (pooled)",
    groups   = c("prevac-Ad26MVA", "prevac-placebo"),
    ge_times = "P+7D",
    ab_cols  = "ab_p_365",
    pooled   = TRUE
  ),
  list(
    analysis = "PREVAC rVSV + placebo (pooled)",
    groups   = c("prevac-rVSV", "prevac-placebo"),
    ge_times = "P+7D",
    ab_cols  = "ab_p_180",
    pooled   = TRUE
  ),
  list(
    analysis = "EBOVAC2 Ad26/MVA (paired)",
    groups   = "ebovac2-Ad26MVA",
    ge_times = c("P+0D", "P+7D"),
    ab_cols  = c("ab_p_0", "ab_p_365"),
    pooled   = FALSE
  ),
  list(
    analysis = "Hamburg rVSV (paired)",
    groups   = "hamburg-rVSV",
    ge_times = c("P+0D", "P+7D"),
    ab_cols  = c("ab_p_0", "ab_p_180"),
    pooled   = FALSE
  )
)

# ---- Build the table: one row per group, plus a Total row for pooled analyses ----

build_plan_rows <- function(plan) {
  eligible <- eligible_participants(plan$groups, plan$ge_times, plan$ab_cols)

  group_rows <- tibble(study_vaccine = plan$groups) %>%
    left_join(
      eligible %>% count(study_vaccine, name = "n"),
      by = "study_vaccine"
    ) %>%
    mutate(n = replace_na(n, 0)) %>%
    transmute(
      Analysis = plan$analysis,
      Group    = group_display_labels[study_vaccine],
      N        = n
    )

  if (plan$pooled) {
    total_row <- tibble(
      Analysis = plan$analysis,
      Group    = "Total (pooled)",
      N        = nrow(eligible)
    )
    group_rows <- bind_rows(group_rows, total_row)
  }

  group_rows
}

sample_size_table <- map_dfr(analysis_plans, build_plan_rows)

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
