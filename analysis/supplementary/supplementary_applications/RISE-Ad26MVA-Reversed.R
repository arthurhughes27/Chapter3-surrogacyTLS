# Supplementary. RISE application: Ad26MVA (EBOVAC2, paired pre/post
# screening) evaluated on Ad26MVA/placebo (PREVAC, independent-groups
# evaluation) - the reverse of analysis/application/RISE-Ad26MVA.R,
# which screens on PREVAC and evaluates on EBOVAC2. Screening/evaluation
# datasets are swapped; screening is paired (matching EBOVAC2's paired
# pre/post design) and evaluation is independent (matching PREVAC's
# two-arm design) - the opposite of the main application. All other
# parameters (screen_power, eval_power, alpha, weighting, correction,
# n_cores) are unchanged from RISE-Ad26MVA.R.

library(SurrogateRank)
library(tidyverse)
library(clipr)
library(purrr)
library(xtable)
library(egg)
library(grid)

source(fs::path("analysis", "config.R"))
source(fs::path("R", "rise_helpers.R"))

set.seed(RISE_SEED)

figure_path <- fs::path(output_path, "figures", "supplementary")
results_path <- fs::path(output_path, "results", "supplementary")
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")
p_load_BTM <- fs::path(processed_data_path, "BTM_processed.rds")

df_merged_all <- readRDS(p_load_merged_all)
BTM <- readRDS(p_load_BTM)

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

# Identify participants with entries at BOTH P+0D and P+7D
participants_both_times <- df_merged_all %>%
  filter(time %in% c("P+0D", "P+7D")) %>%
  distinct(participant_id, time) %>%
  count(participant_id) %>%
  filter(n == 2) %>%
  pull(participant_id)

# Screening (train) dataset: EBOVAC2, paired pre/post within the
# Ad26MVA arm (same participants at P+0D and P+7D) - this was the main
# application's evaluation dataset.
screen_df <- df_merged_all %>%
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

# Evaluation (test) dataset: PREVAC, independent Ad26MVA/placebo groups
# - this was the main application's screening dataset.
eval_df <- df_merged_all %>%
  filter(
    study_accession == "prevac",
    group %in% c("Ad26MVA", "placebo"),
    time == "P+7D",
    !is.na(ab_p_180)
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(participant_id, group, ab_p_180,
                all_of(common_genes))

yone_train <- screen_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_365)

yzero_train <- screen_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

sone_train <- screen_df %>%
  filter(time == "P+7D") %>%
  dplyr::select(all_of(common_genes))

szero_train <- screen_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(all_of(common_genes))

yone_test <- eval_df %>%
  filter(group == "Ad26MVA") %>%
  pull(ab_p_180)

yzero_test <- eval_df %>%
  filter(group == "placebo") %>%
  pull(ab_p_180)

sone_test <- eval_df %>%
  filter(group == "Ad26MVA") %>%
  dplyr::select(all_of(common_genes))

szero_test <- eval_df %>%
  filter(group == "placebo") %>%
  dplyr::select(all_of(common_genes))

rise_res <- run_rise_pipeline(
  yone_train = yone_train, yzero_train = yzero_train,
  sone_train = sone_train, szero_train = szero_train,
  yone_test = yone_test, yzero_test = yzero_test,
  sone_test = sone_test, szero_test = szero_test,
  screen_label = "A) Screening: top 20 markers ",
  eval_label = "B) Evaluation of composite signature",
  figure_path = fs::path(figure_path, "rise_ad26mva_reversed_screening_evaluation.pdf"),
  screen_paired = TRUE, eval_paired = FALSE
)

rise.screen.res <- rise_res$screen
rise.eval.res <- rise_res$evaluate
markers <- rise_res$markers
weights <- rise_res$weights

length(markers)

saveRDS(markers, fs::path(results_path, "TLS_Ad26MVA_ebovac2_reversed.rds"))

# Saved for a future significant-markers table, matching
# analysis/supplementary/full_tables/significant_markers_table_ad26mva.R
# for the main (non-reversed) application.
saveRDS(
  list(screening.metrics = rise.screen.res[["screening.metrics"]],
       screening.weights = rise.screen.res[["screening.weights"]],
       significant.markers = rise.screen.res[["significant.markers"]]),
  fs::path(results_path, "screening_ad26mva_reversed.rds")
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

writeLines(as.character(genelist), con = fs::path(results_path, "genelist_Ad26MVA_Reversed_DAVID.txt"))
writeLines(as.character(backgroundlist), con = fs::path(results_path, "backgroundlist_Ad26MVA_Reversed_DAVID.txt"))

# NOTE: DAVID over-representation results are not yet available for
# this reversed analysis; the corresponding .csv import/table has been
# skipped pending that analysis, matching the main Ad26MVA application.

# For each geneset, find which genes in genelist belong to it
signature_genes <- map_chr(BTM[["genesets"]], function(geneset) {
  matched <- intersect(genelist, geneset)
  paste(matched, collapse = ", ")
})

geneset_table <- data.frame(
  Geneset = BTM[["geneset.names.descriptions"]],
  `Signature genes` = signature_genes,
  check.names = FALSE
)

# Optional: drop genesets with no matching genes from genelist
geneset_table <- geneset_table[geneset_table$`Signature genes` != "", ]

print(
  xtable(geneset_table,
         caption = "Genesets and their matching signature genes (Ad26MVA, reversed screening/evaluation)",
         label = "tab:geneset_signature_genes_ad26mva_reversed"),
  include.rownames = FALSE
)

rm(list = ls())
