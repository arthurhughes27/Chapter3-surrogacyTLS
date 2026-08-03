# ---- Libraries ----
library(SurrogateRank)
library(tidyverse)
library(clipr)
library(readr)
library(xtable)
library(egg)
library(grid)
library(MVN)

set.seed(18042025)

# ---- Paths ----
processed_data_path <- fs::path("data")
figure_path         <- fs::path("output", "figures", "application")
results_path        <- fs::path("output", "results", "application")

# ---- Load data ----
p_load_merged_all <- fs::path(processed_data_path, "df_merged_all.rds")
df_merged_all <- readRDS(p_load_merged_all)

# ---- Identify participants with entries at BOTH P+0D and P+1D ----
participants_both_times <- df_merged_all %>%
  filter(time %in% c("P+0D", "P+1D")) %>%
  distinct(participant_id, time) %>%
  count(participant_id) %>%
  filter(n == 2) %>%
  pull(participant_id)

# ---- Filter to analysis cohort ----
df <- df_merged_all %>%
  filter(
    study_accession == "SDY1276",
    sex == "Female",
    !is.na(ab_p_28),
    !is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where(~ !any(is.na(.))))

# Number of candidate gene predictors
n_inputs <- df %>%
  dplyr::select(a1cf:zzz3) %>%
  ncol()

# ---- Sample splitting ----
train_ids <- df %>%
  distinct(participant_id) %>%
  slice_sample(prop = 1) %>%
  pull(participant_id)

train_df <- df %>% filter(participant_id %in% train_ids)
test_df  <- df %>% filter(!(participant_id %in% train_ids))

# ---- Extract placebo-arm outcome and gene expression matrix ----
yzero_train <- train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

szero_train <- train_df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3)

common_genes <- colnames(szero_train)

# ============================================================
# Bivariate normality assessment (Henze-Zirkler test)
# ============================================================

# ---- Randomly subsample predictors for the analysis ----
# NOTE: set n_sample < length(common_genes) to actually subsample;
# as written (nsampled = length(common_genes)) this uses ALL genes.
set.seed(12345)
n_sample <- length(common_genes)
subset_genes <- sample(common_genes, n_sample)

# ---- Function: test bivariate normality of (gene, outcome) pair ----
test_biv_normality <- function(gene, y, s_df, arm_label, test = "hz") {
  s_vals <- s_df[[gene]]
  if (is.null(s_vals)) return(NULL)

  dat <- na.omit(data.frame(S = s_vals, Y = y))
  if (nrow(dat) < 5) return(NULL)

  out <- tryCatch(
    mvn(dat, mvn_test = test, scale = TRUE, bootstrap = F)$multivariate_normality,
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  out %>%
    mutate(
      p.value = as.character(p.value),  # force consistent type before row-binding
      gene = gene,
      arm = arm_label
    )
}

# ---- Run test across sampled genes (placebo arm) ----
results_placebo <- map_dfr(
  subset_genes,
  ~ test_biv_normality(.x, yzero_train, szero_train, "placebo"),
  .progress = TRUE
)

all_results <- results_placebo %>%
  mutate(p_value_num = as.numeric(gsub("<", "", p.value)))

# ---- Summarise rejection rates across sampled candidates ----
summary_table <- all_results %>%
  group_by(arm, Test) %>%
  summarise(
    n_tested = n(),
    prop_rejected = mean(p_value_num < 0.05, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_table)

# ---- Identify the "worst" gene (largest HZ statistic = strongest departure) ----
worst_candidates <- all_results %>%
  group_by(gene, arm) %>%
  summarise(max_stat = max(Statistic, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_stat))
print(head(worst_candidates, 10))

worst_gene <- worst_candidates$gene[1]
worst_arm  <- worst_candidates$arm[1]
cat("Worst gene:", worst_gene, "| Arm:", worst_arm, "\n")

# ---- Plot: scatterplot for the worst gene ----
plot_biv_normal_check <- function(gene, y, s_df, arm_label) {
  dat <- data.frame(S = s_df[[gene]], Y = y)

  ggplot(dat, aes(x = S, y = Y)) +
    geom_point(alpha = 0.5) +
    labs(
      title = paste(gene, "-", arm_label),
      x = "Surrogate expression (scaled)",
      y = "Antibody response (day 180)"
    ) +
    theme_minimal()
}

worst_plot <- plot_biv_normal_check(worst_gene, yzero_train, szero_train, worst_arm)
worst_plot
