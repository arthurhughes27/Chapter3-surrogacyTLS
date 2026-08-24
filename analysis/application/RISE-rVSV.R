# RISE application: rVSV/placebo (PREVAC, screening) evaluated on
# rVSV (Hamburg, evaluation)

library(SurrogateRank)
library(tidyverse)
library(clipr)
library(readr)
library(purrr)
library(xtable)
library(egg)
library(grid)

source(fs::path("analysis", "config.R"))
source(fs::path("analysis", "functions", "rise_helpers.R"))

set.seed(RISE_SEED)

figure_path <- fs::path(output_path, "figures", "application")
results_path <- fs::path(output_path, "results", "application")
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")
p_load_BTM <- fs::path(processed_data_path, "BTM_processed.rds")

df_merged_all <- readRDS(p_load_merged_all)
BTM <- readRDS(p_load_BTM)

genelist_prevac <- df_merged_all %>%
  filter(study_accession == "prevac") %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  colnames()

genelist_hamburg <- df_merged_all %>%
  filter(study_accession == "hamburg") %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  colnames()

common_genes <- intersect(genelist_prevac, genelist_hamburg)

# How many input genes
n_inputs <- length(common_genes)

n_inputs

train_df <- df_merged_all %>%
  filter(
    study_accession == "prevac",
    group %in% c("rVSV", "placebo"),
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
    study_accession == "hamburg",
    group == "rVSV",
    time %in% c("P+0D", "P+7D"),
    !is.na(ab_p_180),
    !is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(participant_id, time, ab_p_180, ab_p_0,
                all_of(common_genes))


# Step 2: Filter the long dataframe to include only the sampled participants
yone_train <- train_df %>%
  filter(group == "rVSV") %>%
  pull(ab_p_180)

yzero_train <- train_df %>%
  filter(group == "placebo") %>%
  pull(ab_p_180)

sone_train <- train_df %>%
  filter(group == "rVSV") %>%
  dplyr::select(all_of(common_genes))

szero_train <- train_df %>%
  filter(group == "placebo") %>%
  dplyr::select(all_of(common_genes))

yone_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_180)

yzero_test <- test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

sone_test <- test_df %>%
  filter(time == "P+7D") %>%
  dplyr::select(all_of(common_genes))

szero_test <- test_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(all_of(common_genes))

rise_res <- run_rise_pipeline(
  yone_train = yone_train, yzero_train = yzero_train,
  sone_train = sone_train, szero_train = szero_train,
  yone_test = yone_test, yzero_test = yzero_test,
  sone_test = sone_test, szero_test = szero_test,
  screen_label = "A) Screening: top 20 markers ",
  eval_label = "B) Evaluation of 41-gene signature",
  figure_path = fs::path(figure_path, "Figure3-13.pdf")
)

rise.screen.res <- rise_res$screen
rise.eval.res <- rise_res$evaluate
markers <- rise_res$markers
weights <- rise_res$weights

length(markers)

saveRDS(markers, fs::path(results_path, "TLS_rVSV_prevac.rds"))

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

writeLines(as.character(genelist), con = fs::path(results_path, "genelist_rVSV_DAVID.txt"))
writeLines(as.character(backgroundlist), con = fs::path(results_path, "backgroundlist_rVSV_DAVID.txt"))

david_results <- read_csv(fs::path(results_path, "DAVID_rVSV.csv"))

david_table <- david_results %>%
  select(`Full Term`, Count, FDR) %>%
  rename(`Adjusted p-value` = FDR,
         Term = `Full Term`) %>%
  arrange(`Adjusted p-value`) %>%
  head(n = 10)  %>%
  mutate(`Adjusted p-value` = formatC(`Adjusted p-value`, format = "e", digits = 0))

david_table

print(
  xtable(david_table,
         caption = "DAVID over-representation analysis results (rVSV)",
         label = "tab:david_rVSV"),
  include.rownames = FALSE,
  booktabs = TRUE
)


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
         caption = "Genesets and their matching signature genes",
         label = "tab:geneset_signature_genes"),
  include.rownames = FALSE
)

# rm(list = ls())
