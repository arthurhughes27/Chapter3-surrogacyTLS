# Supplementary. Venn diagrams comparing TLS gene- and geneset-level
# signatures across related analyses, mirroring the cross-vaccine Venn
# diagrams in analysis/application/RISE_signatures.R:
#   - TIV (SDY1276): Female (main) vs Male (supplementary)
#   - rVSV: PREVAC->Hamburg (main, unreversed) vs Hamburg->PREVAC
#     (reversed)
#   - Ad26/MVA: PREVAC->EBOVAC2 (main, unreversed) vs EBOVAC2->PREVAC
#     (reversed)
#   - TIV (SDY1276) strains: cross-strain mean vs each of the 3
#     individual strains (4 sets)
# Each comparison is saved as its own figure.

library(ggVennDiagram)
library(ggplot2)
library(patchwork)

# ---- Paths ----
application_results_path   <- fs::path("output", "results", "application")
supplementary_results_path <- fs::path("output", "results", "supplementary")
figure_path                <- fs::path("output", "figures", "supplementary")

BTM <- readRDS(fs::path("data", "BTM_processed.rds"))
genesets <- BTM[["genesets"]]
names(genesets) <- BTM[["geneset.names.descriptions"]]

# a geneset is "hit" if >=1 of its genes is in the signature
genesets_hit_by <- function(sig_genes, genesets) {
  names(genesets)[vapply(genesets, function(gs) any(gs %in% sig_genes), logical(1))]
}

# ---- Shared plotting style for all Venn diagrams (larger text) ----
venn_theme <- theme(
  base_size = 20,
  legend.position = "right",
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 15),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
  plot.margin = margin(10, 20, 10, 20)  # extra right/left margin for spacing
)
venn_style <- function(n_sets) {
  list(
    scale_fill_gradient(low = "#FFFFFF", high = "#4A7FBF"),
    scale_color_manual(values = rep("#2C3E50", n_sets)),
    venn_theme
  )
}

make_venn <- function(lst, subtitle) {
  ggVennDiagram(lst, label = "count", label_alpha = 0, edge_size = 0.8, label_size = 6.5, set_size = 7) +
    venn_style(length(lst)) +
    labs(title = subtitle) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 10, r = 35, b = 10, l = 35, unit = "pt"))
}

# Build the gene-level + geneset-level 2-panel figure for a named list of
# gene-signature vectors and save it
make_overlap_figure <- function(sig_list, title, out_file, width = 15, height = 7) {
  geneset_lists <- lapply(sig_list, genesets_hit_by, genesets = genesets)

  p1 <- make_venn(sig_list, "A) Gene-level")
  p2 <- make_venn(geneset_lists, "B) Geneset-level")

  p_combined <- (p1 | plot_spacer() | p2) +
    plot_layout(widths = c(1, 0.08, 1)) +  # thin spacer column between plots
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25))
    )

  ggsave(fs::path(figure_path, out_file), p_combined, width = width, height = height, dpi = 300)

  p_combined
}

# =============================================================================
# TIV (SDY1276): Female (main) vs Male (supplementary)
# =============================================================================

make_overlap_figure(
  sig_list = list(
    Female = readRDS(fs::path(application_results_path, "TLS_TIV_SDY1276.rds")),
    Male   = readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Male.rds"))
  ),
  title    = "Overlap of TIV (SDY1276) TLS signatures between Females and Males",
  out_file = "rise_signature_overlap_tiv_sex.pdf"
)

# =============================================================================
# rVSV: main (PREVAC screening -> Hamburg evaluation) vs reversed
# (Hamburg screening -> PREVAC evaluation)
# =============================================================================

make_overlap_figure(
  sig_list = list(
    Unreversed = readRDS(fs::path(application_results_path, "TLS_rVSV_prevac.rds")),
    Reversed   = readRDS(fs::path(supplementary_results_path, "TLS_rVSV_hamburg_reversed.rds"))
  ),
  title    = "Overlap of rVSV TLS signatures between the main and reversed-order analyses",
  out_file = "rise_signature_overlap_rvsv_reversed.pdf"
)

# =============================================================================
# Ad26/MVA: main (PREVAC screening -> EBOVAC2 evaluation) vs reversed
# (EBOVAC2 screening -> PREVAC evaluation)
# =============================================================================

make_overlap_figure(
  sig_list = list(
    Unreversed = readRDS(fs::path(application_results_path, "TLS_Ad26MVA_prevac.rds")),
    Reversed   = readRDS(fs::path(supplementary_results_path, "TLS_Ad26MVA_ebovac2_reversed.rds"))
  ),
  title    = "Overlap of Ad26/MVA TLS signatures between the main and reversed-order analyses",
  out_file = "rise_signature_overlap_ad26mva_reversed.pdf"
)

# =============================================================================
# TIV (SDY1276) strains: cross-strain mean vs each of the 3 individual
# strains (4 sets)
# =============================================================================

make_overlap_figure(
  sig_list = list(
    Mean        = readRDS(fs::path(application_results_path, "TLS_TIV_SDY1276.rds")),
    Brisbane10  = readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Brisbane10.rds")),
    Brisbane59  = readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Brisbane59.rds")),
    Florida     = readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Florida.rds"))
  ),
  title    = "Overlap of TIV (SDY1276) TLS signatures across the cross-strain mean and individual strains",
  out_file = "rise_signature_overlap_tiv_strains.pdf",
  # A 4-set Venn diagram needs more horizontal room than a 2-set one
  width    = 17
)

rm(list = ls())
