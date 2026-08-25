# Supplementary. Sensitivity analysis for the Ad26MVA application
# (analysis/application/RISE-Ad26MVA.R): effect of fixing the
# screening-stage non-inferiority margin epsilon at several values
# (instead of choosing it adaptively) on the resulting gamma_Gamma
# composite. Data splitting (PREVAC screening, EBOVAC2 evaluation) is
# identical to the main application; only the screening-stage epsilon
# is varied.

library(SurrogateRank)
library(tidyverse)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "sensitivity_helpers.R"))

set.seed(RISE_SEED)

results_path <- fs::path(output_path, "results", "supplementary")
table_path <- fs::path(output_path, "tables", "supplementary", "epsilon_sensitivity_ad26mva.tex")
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")

df_merged_all <- readRDS(p_load_merged_all)

genelist_prevac <- df_merged_all %>%
  filter(study_accession == "prevac") %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  colnames()

genelist_ebovac2 <- df_merged_all %>%
  filter(study_accession == "ebovac2") %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  colnames()

common_genes <- intersect(genelist_prevac, genelist_ebovac2)

train_df <- df_merged_all %>%
  filter(
    study_accession == "prevac",
    group %in% c("Ad26MVA", "placebo"),
    time == "P+7D",
    !is.na(ab_p_180)
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(participant_id, group, ab_p_180,
                all_of(common_genes))

# Identify participants with entries at BOTH P+0D and P+7D
participants_both_times <- df_merged_all %>%
  filter(time %in% c("P+0D", "P+7D")) %>%
  distinct(participant_id, time) %>%
  count(participant_id) %>%
  filter(n == 2) %>%
  pull(participant_id)

test_df <- df_merged_all %>%
  filter(
    study_accession == "ebovac2",
    group == "Ad26MVA",
    time %in% c("P+0D", "P+7D"),
    !is.na(ab_p_365),
    !is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(participant_id, time, ab_p_365, ab_p_0,
                all_of(common_genes))

yone_train <- train_df %>%
  filter(group == "Ad26MVA") %>%
  pull(ab_p_180)

yzero_train <- train_df %>%
  filter(group == "placebo") %>%
  pull(ab_p_180)

sone_train <- train_df %>%
  filter(group == "Ad26MVA") %>%
  dplyr::select(all_of(common_genes))

szero_train <- train_df %>%
  filter(group == "placebo") %>%
  dplyr::select(all_of(common_genes))

yone_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_365)

yzero_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

sone_test <- test_df %>%
  filter(time == "P+7D") %>%
  dplyr::select(all_of(common_genes))

szero_test <- test_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(all_of(common_genes))

sensitivity_df <- run_rise_epsilon_sensitivity(
  yone_train = yone_train, yzero_train = yzero_train,
  sone_train = sone_train, szero_train = szero_train,
  yone_test = yone_test, yzero_test = yzero_test,
  sone_test = sone_test, szero_test = szero_test,
  screen_paired = FALSE, eval_paired = TRUE, eval_power = 0.9, n_cores = 6
)

saveRDS(sensitivity_df, fs::path(results_path, "epsilon_sensitivity_ad26mva.rds"))

build_epsilon_sensitivity_table(
  sensitivity_df, out_path = table_path,
  caption = "Sensitivity analysis evaluating the effect of varying the screening-stage non-inferiority margin epsilon on the Ad26MVA/placebo (PREVAC, screening) surrogate signature, evaluated in Ad26MVA (EBOVAC2).",
  label = "tab:epsilon-sensitivity-ad26mva"
)

rm(list = ls())
