library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

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
    sex == "Female",!is.na(ab_p_28),!is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  filter(time == "P+0D")

# Figure 1: NAb measurements at day 0 and day 28 within strain and mean across strains

# ---- Reshape to long format ----
df_long <- df %>%
  dplyr::select(
    participant_id,
    ab_p_0, ab_p_0_a_brisbane_10_2007, ab_p_0_a_brisbane_59_2007, ab_p_0_b_florida_4_2006,
    ab_p_28, ab_p_28_a_brisbane_10_2007, ab_p_28_a_brisbane_59_2007, ab_p_28_b_florida_4_2006
  ) %>%
  pivot_longer(
    cols = -participant_id,
    names_to = "colname",
    values_to = "value"
  ) %>%
  mutate(
    timepoint = str_extract(colname, "(?<=ab_p_)[0-9]+"),
    timepoint = factor(timepoint, levels = c("0", "28"), labels = c("Day 0", "Day 28")),
    measure = str_remove(colname, "^ab_p_[0-9]+_?"),
    measure = case_when(
      measure == "" ~ "Cross-strain mean",
      str_detect(measure, "brisbane_10")  ~ "A/Brisbane/10/2007",
      str_detect(measure, "brisbane_59")  ~ "A/Brisbane/59/2008",
      str_detect(measure, "florida_4")    ~ "B/Florida/4/2006",
      TRUE ~ measure
    ),
    measure = factor(measure, levels = c("A/Brisbane/10/2007", "A/Brisbane/59/2008",
                                         "B/Florida/4/2006", "Cross-strain mean"))
  ) %>%
  filter(!is.na(value))

# ---- Plot ----
p1 <- ggplot(df_long, aes(x = measure, y = value, fill = measure)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "grey20", fill = "white", alpha = 0.6) +
  geom_jitter(width = 0.06, size = 0.8, alpha = 0.35, colour = "grey20") +
  facet_wrap(~ timepoint) +
  # Default log10 breaks (e.g. 0.3, 1, 3, 10) mix powers of ten with
  # half-decade multiples to fill sparse ranges, which reads as an
  # arbitrary scale. A strict log10 scale also cannot show 0 (log(0) is
  # undefined), so use a pseudo-log transform instead: linear near zero,
  # logarithmic further out, letting 0 sit on the axis alongside clean
  # powers of ten.
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = function(limits) c(0, 10^seq(0, ceiling(log10(limits[2])))),
    labels = scales::label_number()
  ) +
  scale_fill_manual(values = c(
    "A/Brisbane/10/2007" = "#4C72B0",
    "A/Brisbane/59/2008" = "#55A868",
    "B/Florida/4/2006"   = "#C44E52",
    "Cross-strain mean"               = "#8172B2"
  )) +
  labs(
    x = "Response",
    y = "NAb titre",
    title = "Distribution of NAb titres across strains following vaccination with TIV"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    strip.text = element_text(face = "bold", size = 15, colour = "black"),
    strip.background = element_rect(fill = "grey70", colour = NA),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(2.5, "lines"),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "grey98", colour = NA),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 21)
  )

p1

# Save
ggsave(fs::path(figure_path, "sdy1276_nab_response_distribution.pdf"),
       p1, width = 11, height = 5, dpi = 300)

# rm(list = ls())


