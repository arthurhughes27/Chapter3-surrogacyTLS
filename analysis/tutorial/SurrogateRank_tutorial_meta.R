library(SurrogateRank)
library(ggplot2)
library(stringr)
library(cowplot)
library(dplyr)

fig_dir <- fs::path("output", "figures", "tutorial")
fs::dir_create(fig_dir)

# Background colour matching the minted code-block background used in the
# thesis LaTeX, so tutorial figures are visually tied to the code that
# produced them.
tutorial_bg <- "#F0F2FF"
tutorial_bg_theme <- theme(
  plot.background = element_rect(fill = tutorial_bg, color = tutorial_bg),
  legend.background = element_rect(fill = tutorial_bg, color = NA)
)

# minted's frame=lines draws a horizontal rule at the top and bottom of
# each code block only (no side rules). plot.background's element_rect
# border can't do sides independently (it's all four or none), so the
# top/bottom rules are drawn separately as full-width lines at the very
# edges of the rendered figure.
#
# Plots built with coord_fixed() (e.g. gamma.s.plot$screen.plot) reserve
# extra "respected" aspect-ratio space around the panel that
# plot.background's fill does not reach, leaving the true canvas edges
# unpainted underneath. Drawing an explicit full-canvas rectangle first,
# before the plot itself, ensures that gap picks up the background colour.
add_top_bottom_border <- function(plot, colour = "black", size = 1, bg = tutorial_bg) {
  ggdraw() +
    draw_grob(grid::rectGrob(gp = grid::gpar(fill = bg, col = NA))) +
    draw_plot(plot) +
    draw_line(x = c(0, 1), y = c(1, 1), color = colour, size = size) +
    draw_line(x = c(0, 1), y = c(0, 0), color = colour, size = size)
}

# forest.plot (from rise.evaluate.meta()) is a cowplot::plot_grid()
# composite of pre-rendered sub-panels (study labels / forest / p-value
# table), each with its own baked-in white panel background. Because
# those sub-panels are already-drawn grobs by the time forest.plot is
# returned, an outer `+ theme(...)` on the composite can only recolour
# its own outermost background, not the white rects nested inside it.
# Instead, walk the rendered grob tree directly and recolour any white
# fills found, wherever in the composite they occur.
recolor_white_backgrounds <- function(plot, to = tutorial_bg) {
  g <- if (inherits(plot, "ggplot")) ggplotGrob(plot) else plot

  is_white <- function(fill) {
    !is.null(fill) && !is.na(fill) &&
      (identical(fill, "white") || identical(toupper(fill), "#FFFFFF"))
  }

  recolor <- function(grob) {
    if (!is.null(grob$gp) && !is.null(grob$gp$fill)) {
      fill <- grob$gp$fill
      fill[vapply(fill, is_white, logical(1))] <- to
      grob$gp$fill <- fill
    }
    if (!is.null(grob$children) && length(grob$children) > 0) {
      grob$children <- do.call(grid::gList, lapply(grob$children, recolor))
    }
    if (!is.null(grob$grobs)) {
      grob$grobs <- lapply(grob$grobs, recolor)
    }
    grob
  }

  recolor(g)
}

# Simulate multi-study, high-dimensional individual participant data:
# 5 studies, 25 treated / 25 untreated per study, 100 candidate
# surrogates, 10% of which are truly valid
ipd_data <- generate.example.data.highdim.multistudy.ipd(
  M = 5,             # number of studies
  n1 = 25,           # treated individuals per study
  n0 = 25,           # untreated individuals per study
  p = 100,           # number of candidate surrogate markers
  prop_valid = 0.10  # proportion of truly valid surrogates
)

# Screening stage: RISE is applied within each study individually,
# then the resulting study-level effects are combined via
# random-effects meta-analysis to identify markers with consistent
# evidence of surrogacy across studies
screen_meta_result <- rise.screen.meta(
  yone      = ipd_data$y1,
  yzero     = ipd_data$y0,
  sone      = ipd_data$s1,
  szero     = ipd_data$s0,
  studyone  = ipd_data$study1,
  studyzero = ipd_data$study0,
  epsilon.study = 0.2, # non-inferiority margin for within-study screening
  epsilon.meta  = 0.2  # non-inferiority margin for the meta-analysis stage
)

sig_markers    <- screen_meta_result$significant.markers
screen_weights <- screen_meta_result$screening.weights

study_level_res = screen_meta_result$screening.metrics.study %>%
  arrange(marker) %>%
  head(n = 5)

study_level_res

meta_level_res = screen_meta_result$screening.metrics.meta %>%
  arrange(p.adjusted) %>%
  head(n = 5)

meta_level_res  %>% as.data.frame()

cat(length(sig_markers), "markers retained after meta-analytic screening\n")

p1 = add_top_bottom_border(screen_meta_result$gamma.s.plot$screen.plot + tutorial_bg_theme)

p1

ggsave(plot = p1, path = fig_dir,
       filename = "tutorial_screen_meta.pdf",
       height = 15, width = 30, units = "cm")

# Evaluation stage: combine the retained markers into a single
# composite surrogate signature and re-evaluate its surrogacy across
# studies via meta-analysis
eval_meta_result <- rise.evaluate.meta(
  yone      = ipd_data$y1,
  yzero     = ipd_data$y0,
  sone      = ipd_data$s1,
  szero     = ipd_data$s0,
  studyone  = ipd_data$study1,
  studyzero = ipd_data$study0,
  epsilon.study = 0.2,
  epsilon.meta  = 0.2,
  markers = screen_meta_result$significant.markers,
  screening.weights = screen_meta_result$screening.weights
)

# Meta-analytic evaluation results for the composite marker
print(eval_meta_result$evaluation.metrics.meta %>% as.data.frame() )

p2 = add_top_bottom_border(recolor_white_backgrounds(eval_meta_result$gamma.s.plot$forest.plot))

p2

ggsave(plot = p2, path = fig_dir,
       filename = "tutorial_evaluate_meta.pdf",
       height = 15, width = 30, units = "cm")
