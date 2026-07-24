library(SurrogateRank)
library(tidyverse)

set.seed(18042025)

processed_data_path <- fs::path("data")
figure_path <- fs::path("output", "figures", "descriptive")
p_load_merged_all = fs::path(file = fs::path(processed_data_path, "df_merged_all.rds"))

df_merged_all = readRDS(p_load_merged_all)

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
    sex == "Female",!is.na(ab_p_28_b_florida_4_2006),!is.na(ab_p_0_b_florida_4_2006),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# Sample splitting
train_ids <- df %>%
  distinct(participant_id) %>%
  slice_sample(prop = 0.75) %>%
  pull(participant_id)

# Step 2: Filter the long dataframe to include only the sampled participants
train_df <- df %>%
  filter(participant_id %in% train_ids)

test_df <- df %>%
  filter(!(participant_id %in% train_ids))

yone_train = train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28_b_florida_4_2006)

yzero_train = train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0_b_florida_4_2006)

sone_train = train_df %>%
  filter(time == "P+1D") %>%
  dplyr::select(a1cf:zzz3)

szero_train = train_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

rise.screen.res = rise.screen(yone = yone_train,
                              yzero = yzero_train,
                              sone = sone_train,
                              szero = szero_train,
                              alpha = 0.05,
                              power.want.s = 0.9,
                              p.correction = "BH",
                              alternative = "two.sided",
                              paired = TRUE,
                              weight.mode = "diff.epsilon",
                              n.cores = 8)

length(rise.screen.res[["significant.markers"]])

markers = rise.screen.res[["significant.markers"]]
weights = rise.screen.res[["screening.weights"]]

yone_test = test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28_b_florida_4_2006)

yzero_test = test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0_b_florida_4_2006)

sone_test = test_df %>%
  filter(time == "P+1D") %>%
  dplyr::select(a1cf:zzz3)

szero_test = test_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

rise.eval.res = rise.evaluate(
  yone = yone_test,
  yzero = yzero_test,
  sone = sone_test,
  szero = szero_test,
  alpha = 0.05,
  power.want.s = 0.9,
  p.correction = "BH",
  alternative = "two.sided",
  paired = TRUE,
  n.cores = 8,
  markers = markers,
  screening.weights = weights
)

rise.eval.res[["gamma.s.evaluate"]]
rise.eval.res[["gamma.s.plot"]]
