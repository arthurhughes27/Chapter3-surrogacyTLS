# Supplementary. Sensitivity analysis for the SDY1276 (TIV, Male)
# application
# (analysis/supplementary/supplementary_applications/RISE-SDY1276-Male.R):
# effect of fixing the screening-stage non-inferiority margin epsilon at
# several values (instead of choosing it adaptively) on the resulting
# gamma_Gamma composite. Data loading and the train/test sample split
# are identical to that application; only the screening-stage epsilon
# is varied.

library(SurrogateRank)
library(tidyverse)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "sensitivity_helpers.R"))

set.seed(RISE_SEED)

results_path <- fs::path(output_path, "results", "supplementary")
table_path <- fs::path(output_path, "tables", "supplementary", "epsilon_sensitivity_tiv_male.tex")
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")

df_merged_all <- readRDS(p_load_merged_all)

# Identify participants with entries at BOTH P+0D and P+1D
participants_both_times <- df_merged_all %>%
  filter(time %in% c("P+0D", "P+1D")) %>%
  distinct(participant_id, time) %>%
  count(participant_id) %>%
  filter(n == 2) %>%
  pull(participant_id)

df <- df_merged_all %>%
  filter(
    study_accession == "SDY1276",
    sex == "Male",!is.na(ab_p_28),!is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# Sample splitting
train_ids <- df %>%
  distinct(participant_id) %>%
  slice_sample(prop = 0.75) %>%
  pull(participant_id)

train_df <- df %>%
  filter(participant_id %in% train_ids)

test_df <- df %>%
  filter(!(participant_id %in% train_ids))

yone_train <- train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28)

yzero_train <- train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

sone_train <- train_df %>%
  filter(time == "P+1D") %>%
  dplyr::select(a1cf:zzz3)

szero_train <- train_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

yone_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28)

yzero_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

sone_test <- test_df %>%
  filter(time == "P+1D") %>%
  dplyr::select(a1cf:zzz3)

szero_test <- test_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

sensitivity_df <- run_rise_epsilon_sensitivity(
  yone_train = yone_train, yzero_train = yzero_train,
  sone_train = sone_train, szero_train = szero_train,
  yone_test = yone_test, yzero_test = yzero_test,
  sone_test = sone_test, szero_test = szero_test,
  screen_paired = TRUE, eval_paired = TRUE, eval_power = 0.9, n_cores = 8
)

saveRDS(sensitivity_df, fs::path(results_path, "epsilon_sensitivity_tiv_male.rds"))

build_epsilon_sensitivity_table(
  sensitivity_df, out_path = table_path,
  caption = "Sensitivity analysis evaluating the effect of varying the screening-stage non-inferiority margin epsilon on the SDY1276 (TIV, male) surrogate signature.",
  label = "tab:epsilon-sensitivity-tiv-male"
)

rm(list = ls())
