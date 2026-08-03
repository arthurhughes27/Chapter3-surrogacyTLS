library(SurrogateRank)
library(tidyverse)
library(clipr)
library(readr)
library(xtable)
library(egg)
library(grid)

set.seed(18042025)

processed_data_path <- fs::path("data")
figure_path <- fs::path("output", "figures", "application")
results_path = fs::path("output", "results", "application")
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
    sex == "Female",!is.na(ab_p_28),!is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# How many input genes
n_inputs = df %>%
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

yone_train = train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28)

yzero_train = train_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

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
                       n.cores = 8,
                       screen.plot.topN = 20)

length(rise.screen.res[["significant.markers"]])

markers = rise.screen.res[["significant.markers"]]
weights = rise.screen.res[["screening.weights"]]

saveRDS(markers, fs::path("output", "results", "application","TLS_TIV_SDY1276.rds"))

p1 = rise.screen.res$plot$screen.plot

p1

yone_test = test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_28)

yzero_test = test_df %>%
  filter(time == "P+0D") %>%
  pull(ab_p_0)

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
p2 = rise.eval.res[["gamma.s.plot"]]

p2

p1_labeled <- p1 +
  labs(title = "A) Screening: top 20 markers ") +
  theme(
    plot.title = element_text(face = "bold", size = 25, hjust = 0),
    plot.title.position = "panel",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

p2_labeled <- p2 +
  labs(title = "B) Evaluation of 502-gene signature") +
  theme(
    plot.title = element_text(face = "bold", size = 25, hjust = 0),
    plot.title.position = "panel",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

# Thin vertical divider between panels
divider <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1),
               colour = "white", linewidth = 0.6) +
  theme_void()

g <- egg::ggarrange(
  p1_labeled, divider, p2_labeled,
  ncol   = 3,
  widths = c(1, 0.2, 1),
  draw   = FALSE
)

grid.newpage()
grid.draw(g)

ggsave(fs::path(figure_path, "Figure3-11.pdf"),
       g, width = 17, height = 6, dpi = 300)










# Check the spearman correlation
gamma_Gamma = c(rise.eval.res[["gamma.s"]][["gamma.s.one"]], rise.eval.res[["gamma.s"]][["gamma.s.zero"]])
y_all = c(yone_test, yzero_test)

cor(gamma_Gamma, y_all, method = "spearman")

# Intra-class
cor(rise.eval.res[["gamma.s"]][["gamma.s.one"]], yone_test, method = "spearman")
cor(rise.eval.res[["gamma.s"]][["gamma.s.zero"]], yzero_test, method = "spearman")

# DAVID analysis
genelist = rise.screen.res[["significant.markers"]]
backgroundlist = colnames(sone_train)

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
