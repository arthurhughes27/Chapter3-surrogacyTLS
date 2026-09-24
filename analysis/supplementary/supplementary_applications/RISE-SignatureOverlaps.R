# Supplementary. Venn diagrams comparing TLS gene- and geneset-level
# signatures across related analyses, mirroring the cross-vaccine Venn
# diagrams in analysis/application/RISE_signatures.R:
#   - TIV (SDY1276): Female (main) vs Male (supplementary)
#   - rVSV: PREVAC->Hamburg (main) vs Hamburg->PREVAC (reversed)
#   - Ad26/MVA: PREVAC->EBOVAC2 (main) vs EBOVAC2->PREVAC (reversed)
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

# For a single pairwise comparison (exactly 2 sets), list the elements
# that are shared vs. distinct to each set, at both the gene and
# geneset level - i.e. spelling out in citable text exactly what a
# 2-circle Venn diagram from make_overlap_figure() shows. Not defined
# for > 2 sets, since "shared"/"distinct" stops being a simple
# three-way split once there are more than two circles.
#
# Returns (invisibly) a nested list of the six element vectors
# (gene/geneset x only-first/only-second/shared) and, if out_file is
# given, also writes the same summary to a text file.
describe_overlap <- function(sig_list, title, out_file = NULL) {
  if (length(sig_list) != 2) {
    message(
      "describe_overlap() only defines shared/distinct elements for ",
      "exactly 2 sets (got ", length(sig_list), "); skipping."
    )
    return(invisible(NULL))
  }

  set_names <- names(sig_list)
  geneset_list <- lapply(sig_list, genesets_hit_by, genesets = genesets)

  describe_level <- function(lst, level_name) {
    a <- lst[[1]]
    b <- lst[[2]]
    list(
      level    = level_name,
      shared   = sort(intersect(a, b)),
      only_a   = sort(setdiff(a, b)),
      only_b   = sort(setdiff(b, a))
    )
  }

  levels_out <- list(
    gene    = describe_level(sig_list, "Gene-level"),
    geneset = describe_level(geneset_list, "Geneset-level")
  )

  format_block <- function(block) {
    fmt_set <- function(x) if (length(x) == 0) "(none)" else paste(x, collapse = ", ")
    paste0(
      block$level, " overlap: ", set_names[1], " vs ", set_names[2], "\n",
      "  Shared (", length(block$shared), "): ", fmt_set(block$shared), "\n",
      "  Only in ", set_names[1], " (", length(block$only_a), "): ", fmt_set(block$only_a), "\n",
      "  Only in ", set_names[2], " (", length(block$only_b), "): ", fmt_set(block$only_b), "\n"
    )
  }

  text_out <- paste0(
    title, "\n",
    strrep("-", nchar(title)), "\n",
    format_block(levels_out$gene), "\n",
    format_block(levels_out$geneset)
  )

  cat(text_out, "\n")

  if (!is.null(out_file)) {
    writeLines(text_out, fs::path(figure_path, out_file))
  }

  invisible(levels_out)
}

# =============================================================================
# TIV (SDY1276): Female (main) vs Male (supplementary)
# =============================================================================

sig_list_tiv_sex <- list(
  Female = readRDS(fs::path(application_results_path, "TLS_TIV_SDY1276.rds")),
  Male   = readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Male.rds"))
)
title_tiv_sex <- "Overlap of TIV (SDY1276) TLS signatures between Females and Males"

make_overlap_figure(sig_list_tiv_sex, title_tiv_sex, "rise_signature_overlap_tiv_sex.pdf")
describe_overlap(sig_list_tiv_sex, title_tiv_sex, "rise_signature_overlap_tiv_sex.txt")

# =============================================================================
# rVSV: main (PREVAC screening -> Hamburg evaluation) vs reversed
# (Hamburg screening -> PREVAC evaluation)
# =============================================================================

sig_list_rvsv <- list(
  Main     = readRDS(fs::path(application_results_path, "TLS_rVSV_prevac.rds")),
  Reversed = readRDS(fs::path(supplementary_results_path, "TLS_rVSV_hamburg_reversed.rds"))
)
title_rvsv <- "Overlap of rVSV TLS signatures between the main and reversed-order analyses"

make_overlap_figure(sig_list_rvsv, title_rvsv, "rise_signature_overlap_rvsv_reversed.pdf")
describe_overlap(sig_list_rvsv, title_rvsv, "rise_signature_overlap_rvsv_reversed.txt")

# =============================================================================
# Ad26/MVA: main (PREVAC screening -> EBOVAC2 evaluation) vs reversed
# (EBOVAC2 screening -> PREVAC evaluation)
# =============================================================================

sig_list_ad26mva <- list(
  Main     = readRDS(fs::path(application_results_path, "TLS_Ad26MVA_prevac.rds")),
  Reversed = readRDS(fs::path(supplementary_results_path, "TLS_Ad26MVA_ebovac2_reversed.rds"))
)
title_ad26mva <- "Overlap of Ad26/MVA TLS signatures between the main and reversed-order analyses"

make_overlap_figure(sig_list_ad26mva, title_ad26mva, "rise_signature_overlap_ad26mva_reversed.pdf")
describe_overlap(sig_list_ad26mva, title_ad26mva, "rise_signature_overlap_ad26mva_reversed.txt")

# =============================================================================
# TIV (SDY1276) strains: cross-strain mean vs each of the 3 individual
# strains (4 sets). describe_overlap() is not called here since it only
# defines shared/distinct elements for exactly 2 sets.
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
