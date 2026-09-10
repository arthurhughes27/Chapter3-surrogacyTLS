# =============================================================================
# Script to perform basic study descriptions of PREVAC, Hamburg, EBOVAC2,
# and SDY1276 (TIV)
# =============================================================================

# ---- Libraries ----
library(tidyverse)   # loads dplyr, tidyr, ggplot2, stringr, etc.
library(fs)
library(patchwork)

# ---- Paths ----
processed_data_path      <- fs::path("data")
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# ---- Load harmonised clinical data ----
df_clinical_all <- readRDS(fs::path(processed_data_path, "df_clinical_all.rds"))

# ---- Shared ordering / labelling helpers ----

# Study-group display order, shared by both figures so the two panels
# line up on an identical, harmonised y-axis. EBOVAC2-placebo is
# excluded from both panels. SDY1276 (TIV) is included alongside the
# Ebolavirus vaccine groups.
group_order <- c("prevac-rVSV", "prevac-Ad26MVA", "prevac-placebo",
                 "ebovac2-Ad26MVA", "hamburg-rVSV", "SDY1276-Influenza (IN)")

# y-axis labels: acronyms for PREVAC/EBOVAC2/SDY1276, title case for
# Hamburg, Ad26MVA -> Ad26/MVA and Influenza (IN) -> TIV for readability
group_label_fun <- function(x) {
  study   <- str_extract(x, "^[^-]+")
  vaccine <- str_remove(x, "^[^-]+-") %>%
    str_replace("Ad26MVA", "Ad26/MVA") %>%
    str_replace(fixed("Influenza (IN)"), "TIV")
  study_label <- case_when(
    str_to_lower(study) == "prevac"   ~ "PREVAC",
    str_to_lower(study) == "ebovac2"  ~ "EBOVAC2",
    str_to_lower(study) == "sdy1276"  ~ "SDY1276",
    TRUE ~ str_to_title(study)               # e.g. "hamburg" -> "Hamburg"
  )
  paste0(study_label, " ", vaccine)
}

# =============================================================================
# Panel A: Antibody measurement availability by study-group and timepoint
#
# Each panel shows all available samples for its own assay, independently
# of whether the participant also has a sample for the other assay -
# i.e. no cross-assay eligibility filtering between the two figures.
# =============================================================================

# Desired chronological order: prime then boost
time_order <- c(
  "p_0", "p_7", "p_14", "p_28", "p_56", "p_63", "p_84", "p_180", "p_365"
)

df_counts <- df_clinical_all %>%
  filter(study_vaccine %in% group_order) %>%
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

# =============================================================================
# Panel B: Gene expression sample availability by study-group and timepoint
# =============================================================================

# Explicit chronological order of GE timepoints
# NOTE: label_fun_ge below handles both "P+" and "B+" prefixes, but this
# vector currently only lists prime timepoints. If boost GE timepoints
# (e.g. "B+0D") exist in `time`, add them here or they'll be dropped by
# factor(time, levels = ge_time_order) turning them into NA.
ge_time_order <- c("P+0D", "P+3H", "P+1D", "P+3D", "P+7D")

df_counts_ge <- df_clinical_all %>%
  filter(study_vaccine %in% group_order, !is.na(time)) %>%
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

# =============================================================================
# Combined, harmonised figure
# =============================================================================

# Both panels share one colour scale (same 0-to-max range across the two
# datasets) so tile darkness is directly comparable, and one legend,
# collected by patchwork below.
shared_fill_max <- max(df_counts$n_participants, df_counts_ge$n_participants)

# Shared heatmap theme/style, reused for both panels
heatmap_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 13, face = "bold"),
    axis.title.x    = element_text(size = 14, margin = margin(t = 10)),
    plot.title      = element_text(size = 17, face = "bold", hjust = 0.5),
    panel.grid      = element_blank(),
    plot.margin     = margin(15, 15, 15, 15)
  )

# Build a green heatmap tile + label layer, reused for both panels, on a
# shared fill scale so the two panels are visually comparable. Empty
# (0-count) cells are recoloured light grey with no numeric label,
# rather than sitting at the pale end of the green gradient.
heatmap_layers <- function(df, fill_max) {
  list(
    geom_tile(aes(fill = ifelse(n_participants == 0, NA_real_, n_participants)),
              color = "white", linewidth = 0.8),
    geom_text(
      aes(label = ifelse(n_participants > 0, n_participants, "")),
      size = 4.5, fontface = "bold",
      color = ifelse(df$n_participants > fill_max * 0.55, "white", "grey20")
    ),
    scale_fill_gradient(
      low = "#E8F5E9", high = "#237A21",
      limits = c(0, fill_max),
      na.value = "grey92",
      name = "Samples"
    )
  )
}

p1 <- ggplot(df_counts, aes(x = timepoint, y = study_vaccine)) +
  heatmap_layers(df_counts, shared_fill_max) +
  scale_x_discrete(labels = label_fun) +
  scale_y_discrete(labels = group_label_fun) +
  labs(x = "Timepoint", y = NULL, title = "Antibody measurements") +
  heatmap_theme

p2 <- ggplot(df_counts_ge, aes(x = timepoint, y = study_vaccine)) +
  heatmap_layers(df_counts_ge, shared_fill_max) +
  scale_x_discrete(labels = label_fun_ge) +
  scale_y_discrete(labels = group_label_fun) +
  labs(x = "Timepoint", y = NULL, title = "Gene expression samples") +
  heatmap_theme

p_combined <- (p1 / p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Availability of Antibody and Gene Expression Measurements",
    theme = theme(
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5)
    )
  )

p_combined

ggsave(
  filename = "vaccine_sample_availability.pdf",
  path = descriptive_figures_folder,
  plot = p_combined,
  width = 26, height = 26, units = "cm"
)

rm(list = ls())
