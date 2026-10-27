# =============================================================================
# Presentation-ready screening + evaluation figures for the TIV (SDY1276)
# application, adapted from analysis/application/RISE-SDY1276.R for use in
# slides. Exports panel A (screening) and panel B (evaluation) as two
# SEPARATE PDFs, so they can be shown one after the other.
#
# This re-runs the exact same data preparation and RISE calls as
# RISE-SDY1276.R (same seed, same train/test split, same markers), since
# the thesis figures' underlying ggplot objects / numeric results aren't
# persisted anywhere else -- it then builds two fully custom ggplot2
# figures from the raw results, rather than restyling the package's own
# screen.plot / gamma.s.plot (which bake in choices, e.g. the p-value
# colour scale and coord_fixed margins, that can't be cleanly overridden
# from the outside).
# =============================================================================

library(SurrogateRank)
library(tidyverse)
library(latex2exp)

source(fs::path("analysis", "config.R"))

set.seed(RISE_SEED)

figure_path <- fs::path(output_path, "figures", "presentation")
fs::dir_create(figure_path)

# =============================================================================
# Deck styling -- EDIT THESE to match your slides
# =============================================================================

# Sans font used by the rest of your deck (plain, non-math text: tick
# labels, legend text, italic gene names). latex2exp::TeX() renders the
# maths below via R's plotmath, whose default serif/italic style already
# closely resembles Computer Modern without installing anything extra;
# DECK_SANS_FONT does not affect those maths-containing axis titles,
# since plotmath renders a whole expression in one family (mixing a sans
# word and a Computer Modern symbol within the SAME title isn't
# supported by base R graphics -- a title with maths in it, e.g. panel
# B's axis titles, will render entirely in the serif/maths look, which
# is the standard convention for this anyway).
DECK_SANS_FONT <- "Latin Modern Sans"  # e.g. "Fira Sans", "Lato", "Open Sans"...

# Rendering fonts reliably (especially non-default ones) needs the Cairo
# PDF device rather than the base pdf() device.
PRESENTATION_DEVICE <- grDevices::cairo_pdf

# Final on-slide size: ~14pt text at this physical canvas size is what
# PRESENTATION_BASE_SIZE / ggsave's width+height below are set to
# achieve directly (ggplot2 text sizes are already in points). The
# LaTeX slide scales the figure to ~70%, landing text at ~10pt there.
PRESENTATION_BASE_SIZE <- 14
PRESENTATION_WIDTH_CM  <- 8
PRESENTATION_HEIGHT_CM <- 8

# Colours
SURROGATE_BLUE  <- "#5461E0"
ENDPOINT_RED    <- "#C8504E"
PRE_VACC_GREY   <- "#BDBDBD"
POST_VACC_TEAL  <- "#2A9D8F"

# How many top markers to show in the screening panel (the thesis uses
# 20; gene names become unreadable at that count on a small slide figure)
PRESENTATION_SCREEN_TOP_N <- 10

# =============================================================================
# Re-run the TIV (SDY1276) RISE application (identical to RISE-SDY1276.R)
# =============================================================================

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
    sex == "Female", !is.na(ab_p_28), !is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# Sample splitting (same seed + call as RISE-SDY1276.R -> identical split)
train_ids <- df %>%
  distinct(participant_id) %>%
  slice_sample(prop = 0.75) %>%
  pull(participant_id)

train_df <- df %>% filter(participant_id %in% train_ids)
test_df  <- df %>% filter(!(participant_id %in% train_ids))

yone_train <- train_df %>% filter(time == "P+0D") %>% pull(ab_p_28)
yzero_train <- train_df %>% filter(time == "P+0D") %>% pull(ab_p_0)
sone_train <- train_df %>% filter(time == "P+1D") %>% dplyr::select(a1cf:zzz3)
szero_train <- train_df %>% filter(time == "P+0D") %>% dplyr::select(a1cf:zzz3)

yone_test <- test_df %>% filter(time == "P+0D") %>% pull(ab_p_28)
yzero_test <- test_df %>% filter(time == "P+0D") %>% pull(ab_p_0)
sone_test <- test_df %>% filter(time == "P+1D") %>% dplyr::select(a1cf:zzz3)
szero_test <- test_df %>% filter(time == "P+0D") %>% dplyr::select(a1cf:zzz3)

# Screening (train), matching run_rise_pipeline()'s call in R/rise_helpers.R
rise.screen.res <- rise.screen(
  yone = yone_train, yzero = yzero_train,
  sone = sone_train, szero = szero_train,
  alpha = 0.05, power.want.s = 0.9,
  p.correction = "BH", alternative = "two.sided",
  paired = TRUE, weight.mode = "diff.epsilon",
  n.cores = 8, return.screen.plot = FALSE
)

markers <- rise.screen.res[["significant.markers"]]
weights <- rise.screen.res[["screening.weights"]]

# Evaluation (test), matching run_rise_pipeline()'s call
rise.eval.res <- rise.evaluate(
  yone = yone_test, yzero = yzero_test,
  sone = sone_test, szero = szero_test,
  alpha = 0.05, power.want.s = 0.9,
  p.correction = "BH", alternative = "two.sided",
  paired = TRUE, n.cores = 8,
  markers = markers, screening.weights = weights
)

# =============================================================================
# Panel A: Screening -- top PRESENTATION_SCREEN_TOP_N markers
# =============================================================================

epsilon.val <- unique(rise.screen.res$screening.metrics$epsilon)  # should be ~ +/-0.289

df_screen <- rise.screen.res$screening.metrics %>%
  arrange(p_unadjusted) %>%
  slice_head(n = PRESENTATION_SCREEN_TOP_N) %>%
  mutate(
    marker   = factor(marker, levels = rev(unique(marker))),
    ci_lower = pmax(ci_lower, -1),
    ci_upper = pmin(ci_upper, 1)
  )

shade_df <- data.frame(
  xmin = -epsilon.val,
  xmax =  epsilon.val,
  ymin = 0.5,
  ymax = nrow(df_screen) + 0.5
)

p_screen <- ggplot(df_screen, aes(x = delta, y = marker)) +
  # Light tint for the equivalence band, instead of the thesis's grey
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = SURROGATE_BLUE, alpha = 0.12
  ) +
  # Small label on the equivalence band, if there's room for it
  annotate(
    "text", x = 0, y = nrow(df_screen) + 0.5, vjust = -0.3,
    label = TeX("$\\pm\\hat{\\epsilon}$"), size = 3.2, color = SURROGATE_BLUE
  ) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_segment(
    aes(x = ci_lower, xend = ci_upper, y = marker, yend = marker),
    color = SURROGATE_BLUE, linewidth = 0.9, lineend = "round"
  ) +
  # One colour for all points -- significance is no longer double-encoded
  # via a separate p-value colour scale (all shown markers are significant)
  geom_point(color = SURROGATE_BLUE, size = 2.3) +
  labs(x = TeX("$\\hat{\\delta}$ (surrogacy)"), y = NULL) +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal(base_size = PRESENTATION_BASE_SIZE) +
  theme(
    text          = element_text(family = DECK_SANS_FONT),
    axis.title.x  = element_text(color = SURROGATE_BLUE, size = PRESENTATION_BASE_SIZE),
    axis.text.y   = element_text(face = "italic", size = 12),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "none",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin   = margin(6, 10, 6, 6)
  )

p_screen

ggsave(
  fs::path(figure_path, "tiv_screening_presentation.pdf"),
  p_screen,
  width = PRESENTATION_WIDTH_CM, height = PRESENTATION_HEIGHT_CM, units = "cm",
  device = PRESENTATION_DEVICE
)

# =============================================================================
# Panel B: Evaluation -- composite surrogate vs clinical endpoint
# =============================================================================

rank_df <- rise.eval.res$gamma.s.plot$data  # columns: treatment, response.rank, gamma.rank

p_eval <- ggplot(rank_df, aes(x = gamma.rank, y = response.rank, color = treatment)) +
  geom_point(size = 2, alpha = 0.9) +
  # Dashed identity line: dark grey, not red/pink
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey40", linewidth = 0.8
  ) +
  scale_color_manual(
    name   = NULL,
    breaks = c("0", "1"),
    values = c("0" = PRE_VACC_GREY, "1" = POST_VACC_TEAL),
    labels = c("0" = "Pre-vaccination", "1" = "Post-vaccination")
  ) +
  labs(
    x = TeX("Composite surrogate $\\hat{\\gamma}_{\\Gamma}$ (rank)"),
    y = TeX("Clinical endpoint $Y$ (rank)")
  ) +
  scale_x_continuous(breaks = c(0, 25, 50)) +
  scale_y_continuous(breaks = c(0, 25, 50)) +
  coord_fixed() +
  theme_minimal(base_size = PRESENTATION_BASE_SIZE) +
  theme(
    text         = element_text(family = DECK_SANS_FONT),
    axis.title.x = element_text(color = SURROGATE_BLUE, size = 12),
    axis.title.y = element_text(color = ENDPOINT_RED,   size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 2),
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin   = margin(6, 10, 6, 6)
  )

# The thesis annotates the printed delta/rho directly on the panel; removed
# here since you said you'd put those numbers on the slide itself. To add
# them back (very large, as discussed), uncomment and adjust:
# delta.val <- rise.eval.res$gamma.s.evaluate[["delta"]]
# p.val     <- rise.eval.res$gamma.s.evaluate[["p_unadjusted"]]
# p_eval <- p_eval +
#   annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
#            label = TeX(sprintf("$\\hat{\\delta} = %.2f,\\; p = %.3f$", delta.val, p.val)),
#            size = 10, color = ENDPOINT_RED)

p_eval

ggsave(
  fs::path(figure_path, "tiv_evaluation_presentation.pdf"),
  p_eval,
  width = PRESENTATION_WIDTH_CM, height = PRESENTATION_HEIGHT_CM, units = "cm",
  device = PRESENTATION_DEVICE
)

# rm(list = ls())
