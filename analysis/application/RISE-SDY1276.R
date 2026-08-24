# RISE application: SDY1276 (TIV, sample-split into train/test)

library(SurrogateRank)
library(tidyverse)
library(clipr)
library(readr)
library(xtable)
library(egg)
library(grid)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "rise_helpers.R"))

set.seed(RISE_SEED)

figure_path <- fs::path(output_path, "figures", "application")
results_path <- fs::path(output_path, "results", "application")
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
    sex == "Female",!is.na(ab_p_28),!is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# How many input genes
n_inputs <- df %>%
  dplyr::select(a1cf:zzz3) %>%
  ncol()

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

rise_res <- run_rise_pipeline(
  yone_train = yone_train, yzero_train = yzero_train,
  sone_train = sone_train, szero_train = szero_train,
  yone_test = yone_test, yzero_test = yzero_test,
  sone_test = sone_test, szero_test = szero_test,
  screen_label = "A) Screening: top 20 markers ",
  eval_label = "B) Evaluation of 502-gene signature",
  figure_path = fs::path(figure_path, "rise_sdy1276_screening_evaluation.pdf"),
  screen_power = 0.9, screen_paired = TRUE, n_cores = 8
)

rise.screen.res <- rise_res$screen
rise.eval.res <- rise_res$evaluate
markers <- rise_res$markers
weights <- rise_res$weights

length(markers)

saveRDS(markers, fs::path(results_path, "TLS_TIV_SDY1276.rds"))

rise.eval.res[["gamma.s.evaluate"]]

# Check the spearman correlation
gamma_Gamma <- c(rise.eval.res[["gamma.s"]][["gamma.s.one"]], rise.eval.res[["gamma.s"]][["gamma.s.zero"]])
y_all <- c(yone_test, yzero_test)

cor(gamma_Gamma, y_all, method = "spearman")

# Intra-class
cor(rise.eval.res[["gamma.s"]][["gamma.s.one"]], yone_test, method = "spearman")
cor(rise.eval.res[["gamma.s"]][["gamma.s.zero"]], yzero_test, method = "spearman")

# DAVID analysis
genelist <- rise.screen.res[["significant.markers"]]
backgroundlist <- colnames(sone_train)

# Copy to clipboard, one gene per line
# write_clip(genelist)
# write_clip(backgroundlist)

writeLines(as.character(genelist), con = fs::path(results_path, "genelist_SDY1276_DAVID.txt"))
writeLines(as.character(backgroundlist), con = fs::path(results_path, "backgroundlist_SDY1276_DAVID.txt"))

david_results <- read_csv(fs::path(results_path, "DAVID_SDY1276.csv"))

david_table <- david_results %>%
  select(`Full Term`, Count, FDR) %>%
  rename(`Adjusted p-value` = FDR,
         Term = `Full Term`) %>%
  head(n = 10)  %>%
  arrange(`Adjusted p-value`) %>%
  mutate(`Adjusted p-value` = formatC(`Adjusted p-value`, format = "e", digits = 0))

david_table

print(
  xtable(david_table,
         caption = "DAVID over-representation analysis results (SDY1276)",
         label = "tab:david_sdy1276"),
  include.rownames = FALSE,
  booktabs = TRUE
)

rm(list = ls())
