# ------------------------------------------------------------------
# Simple demonstration of high-dimensional surrogate marker analysis
# with RISE (rise.screen + rise.evaluate), SurrogateRank package
# ------------------------------------------------------------------
library(SurrogateRank)
library(fs)
library(ggplot2)
library(stringr)
library(cowplot)

set.seed(1234)

# Output directory for figures
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
# Plots built with coord_fixed() (e.g. gamma.s.plot) reserve extra
# "respected" aspect-ratio space around the panel that plot.background's
# fill does not reach, leaving the true canvas edges unpainted underneath.
# Drawing an explicit full-canvas rectangle first, before the plot itself,
# ensures that gap picks up the background colour too.
add_top_bottom_border <- function(plot, colour = "black", size = 1, bg = tutorial_bg) {
  ggdraw() +
    draw_grob(grid::rectGrob(gp = grid::gpar(fill = bg, col = NA))) +
    draw_plot(plot) +
    draw_line(x = c(0, 1), y = c(1, 1), color = colour, size = size) +
    draw_line(x = c(0, 1), y = c(0, 0), color = colour, size = size)
}

# ------------------------------------------------------------------
# Simulate a single high-dimensional dataset
# ------------------------------------------------------------------
# generate.example.data.highdim() creates a two-group (treated vs.
# untreated) primary outcome Y and a set of p candidate surrogate
# markers S, with a known proportion of which are constructed to be
# valid surrogates.
full_data <- generate.example.data.highdim(
  n1 = 50,          # Number of treated individuals
  n0 = 50,          # Number of untreated individuals
  p = 200,          # Number of candidate surrogate markers
  prop_valid = 0.10 # Proportion of the 200 candidates that are truly valid surrogates
)

# Structure of the simulated data: y1/y0 (primary outcome by group),
# s1/s0 (n x p matrices of candidate surrogates by group),
# hyp (surrogate validity)
str(full_data)

# ------------------------------------------------------------------
# Sample-split into a screening and evaluation data 50/50
# ------------------------------------------------------------------
n1 <- length(full_data$y1)
n0 <- length(full_data$y0)

# Randomly sample IDs from treated and untreated groups
id_treated_screen   <- sample(seq_len(n1), size = n1 / 2)
id_untreated_screen <- sample(seq_len(n0), size = n0 / 2)
id_treated_evaluate   <- setdiff(seq_len(n1), id_treated_screen)
id_untreated_evaluate <- setdiff(seq_len(n0), id_untreated_screen)

# Build the screening dataset (first half of each group)
screen_data <- list(
  y1 = full_data$y1[id_treated_screen],
  y0 = full_data$y0[id_untreated_screen],
  s1 = full_data$s1[id_treated_screen, ],
  s0 = full_data$s0[id_untreated_screen, ]
)

# Build the evaluation dataset (remaining half of each group)
eval_data <- list(
  y1 = full_data$y1[id_treated_evaluate],
  y0 = full_data$y0[id_untreated_evaluate],
  s1 = full_data$s1[id_treated_evaluate, ],
  s0 = full_data$s0[id_untreated_evaluate, ]
)

# ------------------------------------------------------------------
# Screening stage: rise.screen()
# ------------------------------------------------------------------
screen_result <- rise.screen(
  yone    = screen_data$y1,
  yzero   = screen_data$y0,
  sone    = screen_data$s1,
  szero   = screen_data$s0,
  power.want.s = 0.8,   # desired power to detect a valid surrogate -> used to derive epsilon
  alpha   = 0.05,       # significance level for the screening tests
  p.correction = "BH"   # Benjamini-Hochberg correction for multiple testing
)

# Examine marker-level screening results
head(screen_result[["screening.metrics"]], n= 5)

# Names of the markers whose adjusted p-value is below alpha
sig_markers <- screen_result$significant.markers

# Per-marker weights
screen_weights <- screen_result$screening.weights

cat(length(sig_markers), "markers retained after screening\n")

# Forest-plot-style summary of the top screened candidates
screen_plot <- add_top_bottom_border(screen_result$plot$screen.plot + tutorial_bg_theme)

print(screen_plot)

ggsave(plot = screen_plot, path = fig_dir,
       filename = "tutorial_screen.pdf",
       height = 15, width = 30, units = "cm")

# ------------------------------------------------------------------
# Evaluation stage: rise.evaluate()
# ------------------------------------------------------------------
eval_result <- rise.evaluate(
  yone    = eval_data$y1,
  yzero   = eval_data$y0,
  sone    = eval_data$s1,
  szero   = eval_data$s0,
  power.want.s = 0.8,
  alpha   = 0.05,
  markers = sig_markers,             # markers selected during screening
  screening.weights = screen_weights # their corresponding weights
)

# Surrogacy test results for the composite marker gamma
print(eval_result$gamma.s.evaluate)

# Rank-scale plot of the composite surrogate against the primary
# response, illustrating the strength of the association
eval_plot <- add_top_bottom_border(eval_result$gamma.s.plot + tutorial_bg_theme)

print(eval_plot)

ggsave(plot = eval_plot, path = fig_dir,
       filename = "tutorial_evaluate.pdf",
       height = 15, width = 30, units = "cm")
