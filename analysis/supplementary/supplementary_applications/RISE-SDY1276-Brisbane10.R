# Supplementary. RISE application: SDY1276 (TIV, sample-split into
# train/test), analysing the A/Brisbane/10/2007 strain-specific
# antibody response instead of the main application's
# mean-across-strains response (ab_p_28/ab_p_0) - otherwise identical
# to analysis/application/RISE-SDY1276.R (same cohort: sex == "Female").

library(SurrogateRank)
library(tidyverse)
library(clipr)
library(readr)
library(xtable)
library(egg)
library(grid)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "rise_helpers.R"))
source(fs::path("R", "sdy1276_strain_helpers.R"))

set.seed(RISE_SEED)

figure_path <- fs::path(output_path, "figures", "supplementary")
results_path <- fs::path(output_path, "results", "supplementary")
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")

df_merged_all <- readRDS(p_load_merged_all)

yone_col <- find_strain_column(colnames(df_merged_all), 28, "brisbane_10")
yzero_col <- find_strain_column(colnames(df_merged_all), 0, "brisbane_10")

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
    sex == "Female",
    !is.na(.data[[yone_col]]),
    !is.na(.data[[yzero_col]]),
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

train_df <- df %>%
  filter(participant_id %in% train_ids)

test_df <- df %>%
  filter(!(participant_id %in% train_ids))

yone_train <- train_df %>%
  filter(time == "P+0D") %>%
  pull(.data[[yone_col]])

yzero_train <- train_df %>%
  filter(time == "P+0D") %>%
  pull(.data[[yzero_col]])

sone_train <- train_df %>%
  filter(time == "P+1D") %>%
  dplyr::select(a1cf:zzz3)

szero_train <- train_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

yone_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(.data[[yone_col]])

yzero_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(.data[[yzero_col]])

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
  eval_label = "B) Evaluation of composite signature",
  figure_path = fs::path(figure_path, "rise_sdy1276_brisbane10_screening_evaluation.pdf"),
  screen_power = 0.9, screen_paired = TRUE, n_cores = 8
)

rise.screen.res <- rise_res$screen
rise.eval.res <- rise_res$evaluate
markers <- rise_res$markers
weights <- rise_res$weights

length(markers)

saveRDS(markers, fs::path(results_path, "TLS_TIV_SDY1276_Brisbane10.rds"))

# Saved for a future significant-markers table, matching
# analysis/supplementary/full_tables/significant_markers_table_sdy1276.R
# for the main (mean-across-strains) application.
saveRDS(
  list(screening.metrics = rise.screen.res[["screening.metrics"]],
       screening.weights = rise.screen.res[["screening.weights"]],
       significant.markers = rise.screen.res[["significant.markers"]]),
  fs::path(results_path, "screening_sdy1276_brisbane10.rds")
)

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

writeLines(as.character(genelist), con = fs::path(results_path, "genelist_SDY1276_Brisbane10_DAVID.txt"))
writeLines(as.character(backgroundlist), con = fs::path(results_path, "backgroundlist_SDY1276_Brisbane10_DAVID.txt"))

# NOTE: DAVID over-representation results are not yet available for
# this strain-specific analysis; the corresponding .csv import/table
# has been skipped pending that analysis.

rm(list = ls())
