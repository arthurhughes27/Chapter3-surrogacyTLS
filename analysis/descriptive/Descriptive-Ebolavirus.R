# =============================================================================
# Script to perform basic study descriptions of PREVAC, Hamburg, EBOVAC2
# =============================================================================

# ---- Libraries ----
library(tidyverse)   # loads dplyr, tidyr, ggplot2, stringr, etc.
library(fs)

# ---- Paths ----
processed_data_path      <- fs::path("data")
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# ---- Load harmonised clinical data ----
df_clinical_all <- readRDS(fs::path(processed_data_path, "df_clinical_all.rds")) %>%
  filter(study_vaccine != "SDY1276-Influenza (IN)")

# ---- Shared ordering / labelling helpers ----

# Study-group display order (used for both antibody and GE figures)
group_order <- c("prevac-rVSV", "prevac-Ad26MVA", "prevac-placebo",
                 "ebovac2-Ad26MVA", "ebovac2-placebo", "hamburg-rVSV")

# y-axis labels: acronyms for PREVAC/EBOVAC2, title case for Hamburg,
# Ad26MVA -> Ad26/MVA for readability
group_label_fun <- function(x) {
  study   <- str_extract(x, "^[^-]+")
  vaccine <- str_remove(x, "^[^-]+-") %>%
    str_replace("Ad26MVA", "Ad26/MVA")
  study_label <- case_when(
    str_to_lower(study) == "prevac"  ~ "PREVAC",
    str_to_lower(study) == "ebovac2" ~ "EBOVAC2",
    TRUE ~ str_to_title(study)               # e.g. "hamburg" -> "Hamburg"
  )
  paste0(study_label, " ", vaccine)
}

# Shared heatmap theme/style, reused for both figures
heatmap_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 13, face = "bold"),
    axis.title.x    = element_text(size = 14, margin = margin(t = 10)),
    plot.title      = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 12, colour = "grey40", hjust = 0.5,
                                   margin = margin(b = 12)),
    panel.grid      = element_blank(),
    plot.margin     = margin(15, 15, 15, 15)
  )

# Build a green heatmap tile + label layer, reused for both figures.
# df must contain n_participants; text colour adapts to tile darkness.
heatmap_layers <- function(df) {
  list(
    geom_tile(color = "white", linewidth = 0.8),
    geom_text(
      aes(label = ifelse(n_participants > 0, n_participants, "")),
      size = 4.5, fontface = "bold",
      color = ifelse(df$n_participants > max(df$n_participants) * 0.55,
                     "white", "grey20")
    ),
    scale_fill_gradient(low = "white", high = "#237A21", guide = "none")
  )
}

# =============================================================================
# Figure 1: Antibody measurement availability by study-group and timepoint
# =============================================================================

# Desired chronological order: prime then boost
time_order <- c(
  "p_0", "p_7", "p_14", "p_28", "p_56", "p_63", "p_84", "p_180", "p_365"
)

df_counts <- df_clinical_all %>%
  pivot_longer(
    cols = starts_with("ab_"),
    names_to = "ab_time",
    values_to = "ab_value"
  ) %>%
  mutate(
    phase         = str_extract(ab_time, "(?<=ab_)[pb]"),
    day           = str_extract(ab_time, "\\d+$"),
    timepoint     = factor(paste0(phase, "_", day), levels = time_order),
    study_vaccine = factor(study_vaccine, levels = rev(group_order))
  ) %>%
  filter(!is.na(ab_value)) %>%
  group_by(study_vaccine, timepoint) %>%
  summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
  complete(study_vaccine, timepoint, fill = list(n_participants = 0)) %>%
  filter(!is.na(timepoint))

# x-axis labels: "Prime + 7D" / "Boost + 7D"
label_fun <- function(x) {
  ifelse(
    str_detect(x, "^p_"),
    paste0("Prime + ", str_remove(x, "^p_"), "D"),
    paste0("Boost + ", str_remove(x, "^b_"), "D")
  )
}

p1 <- ggplot(df_counts, aes(x = timepoint, y = study_vaccine, fill = n_participants)) +
  heatmap_layers(df_counts) +
  scale_x_discrete(labels = label_fun) +
  scale_y_discrete(labels = group_label_fun) +
  labs(
    x = "Timepoint",
    y = NULL,
    title = "Availability of Antibody Measurements for Ebolavirus Vaccines",
    subtitle = "Number of participants with a measurement, by study-group and timepoint"
  ) +
  heatmap_theme

p1

ggsave(
  filename = "Figure3-9.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 26, height = 15, units = "cm"
)

# =============================================================================
# Figure 2: Gene expression sample availability by study-group and timepoint
#           (restricted to participants with a valid immune-response pairing)
# =============================================================================

# Explicit chronological order of GE timepoints
# NOTE: label_fun_ge below handles both "P+" and "B+" prefixes, but this
# vector currently only lists prime timepoints. If boost GE timepoints
# (e.g. "B+0D") exist in `time`, add them here or they'll be dropped by
# factor(time, levels = ge_time_order) turning them into NA.
ge_time_order <- c("P+0D", "P+3H", "P+1D", "P+3D", "P+7D")

group_order <- c("prevac-rVSV", "prevac-Ad26MVA", "prevac-placebo",
                 "ebovac2-Ad26MVA", "hamburg-rVSV")

# Participants with a baseline (P+0D) gene expression sample
participants_baseline_ge <- df_clinical_all %>%
  filter(time == "P+0D") %>%
  distinct(participant_id) %>%
  pull(participant_id)

# Identify participants with a valid immune-response day-0/follow-up pairing
# for their study_vaccine group, AND (for EBOVAC2/Hamburg only) a baseline GE sample:
#   - PREVAC-Ad26MVA:                       day 0 AND day 365 (ab)
#   - EBOVAC2-Ad26MVA, EBOVAC2-placebo:     day 0 AND day 365 (ab) AND baseline GE
#   - PREVAC-rVSV:                          day 0 AND day 180 (ab)
#   - Hamburg-rVSV:                         day 0 AND day 180 (ab) AND baseline GE
#   - PREVAC-placebo:                       day 0 AND (365 OR 180) (ab)
participants_valid_ir <- df_clinical_all %>%
  filter(study_vaccine != "ebovac2-placebo") %>%
  distinct(participant_id, study_vaccine, ab_p_0, ab_p_180, ab_p_365) %>%
  mutate(
    has_valid_ab_pair = case_when(
      study_vaccine %in% c("prevac-Ad26MVA", "ebovac2-Ad26MVA") ~
        !is.na(ab_p_0) & !is.na(ab_p_365),
      study_vaccine %in% c("prevac-rVSV", "hamburg-rVSV") ~
        !is.na(ab_p_0) & !is.na(ab_p_180),
      study_vaccine == "prevac-placebo" ~
        !is.na(ab_p_0) & (!is.na(ab_p_365) | !is.na(ab_p_180)),
      TRUE ~ FALSE
    ),
    requires_baseline_ge = study_vaccine %in% c("ebovac2-Ad26MVA", "hamburg-rVSV"),
    has_baseline_ge = participant_id %in% participants_baseline_ge,
    has_valid_pair = has_valid_ab_pair & (!requires_baseline_ge | has_baseline_ge)
  ) %>%
  filter(has_valid_pair) %>%
  pull(participant_id)

df_counts_ge <- df_clinical_all %>%
  filter(!is.na(time),
         participant_id %in% participants_valid_ir) %>%
  mutate(
    timepoint     = factor(time, levels = ge_time_order),
    study_vaccine = factor(study_vaccine, levels = rev(group_order))
  ) %>%
  group_by(study_vaccine, timepoint) %>%
  summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
  complete(study_vaccine, timepoint, fill = list(n_participants = 0)) %>%
  filter(!is.na(timepoint))

# x-axis labels: "Prime + 3H" / "Boost + 7D"
label_fun_ge <- function(x) {
  x <- as.character(x)
  ifelse(
    str_detect(x, "^P\\+"),
    paste0("Prime + ", str_remove(x, "^P\\+")),
    paste0("Boost + ", str_remove(x, "^B\\+"))
  )
}

p2 <- ggplot(df_counts_ge, aes(x = timepoint, y = study_vaccine, fill = n_participants)) +
  heatmap_layers(df_counts_ge) +
  scale_x_discrete(labels = label_fun_ge) +
  scale_y_discrete(labels = group_label_fun) +
  labs(
    x = "Timepoint",
    y = NULL,
    title = "Availability of Gene Expression Samples for Ebolavirus Vaccines",
    subtitle = "Number of eligible participants with a sample, by study-group and timepoint"
  ) +
  heatmap_theme

p2

ggsave(
  filename = "Figure3-10.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 26, height = 15, units = "cm"
)

rm(list = ls())
