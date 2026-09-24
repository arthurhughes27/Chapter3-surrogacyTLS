# Supplementary. Overlap of the TIV (SDY1276) TLS gene- and
# geneset-level signatures between the main (Female) application and
# the Male-only supplementary application. Mirrors the cross-vaccine
# Venn diagrams in analysis/application/RISE_signatures.R, but compares
# the two sexes for the same vaccine instead of different vaccines.

library(ggVennDiagram)
library(ggplot2)
library(patchwork)

# ---- Paths & data ----
application_results_path    <- fs::path("output", "results", "application")
supplementary_results_path  <- fs::path("output", "results", "supplementary")
figure_path                 <- fs::path("output", "figures", "supplementary")

sig_female <- readRDS(fs::path(application_results_path, "TLS_TIV_SDY1276.rds"))
sig_male   <- readRDS(fs::path(supplementary_results_path, "TLS_TIV_SDY1276_Male.rds"))
BTM        <- readRDS(fs::path("data", "BTM_processed.rds"))

# ---- Build gene- and geneset-level signature lists ----
gene_lists <- list(Female = sig_female, Male = sig_male)

genesets <- BTM[["genesets"]]
names(genesets) <- BTM[["geneset.names.descriptions"]]

# a geneset is "hit" if >=1 of its genes is in the signature
genesets_hit_by <- function(sig_genes, genesets) {
  names(genesets)[vapply(genesets, function(gs) any(gs %in% sig_genes), logical(1))]
}
geneset_lists <- lapply(gene_lists, genesets_hit_by, genesets = genesets)

# ---- Shared plotting style for both Venn diagrams (larger text) ----
venn_theme <- theme(
  base_size = 20,
  legend.position = "right",
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 15),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
  plot.margin = margin(10, 20, 10, 20)  # extra right/left margin for spacing
)
venn_style <- list(
  scale_fill_gradient(low = "#FFFFFF", high = "#4A7FBF"),
  scale_color_manual(values = rep("#2C3E50", 2)),
  venn_theme
)

make_venn <- function(lst, subtitle) {
  ggVennDiagram(lst, label = "count", label_alpha = 0, edge_size = 0.8, label_size = 6.5, set_size = 7) +
    venn_style +
    labs(title = subtitle) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 10, r = 35, b = 10, l = 35, unit = "pt"))
}

p1 <- make_venn(gene_lists, "A) Gene-level")
p2 <- make_venn(geneset_lists, "B) Geneset-level")

# ---- Combine side by side, with blank spacer + one shared, larger title ----
p_combined <- (p1 | plot_spacer() | p2) +
  plot_layout(widths = c(1, 0.08, 1)) +  # thin spacer column between plots
  plot_annotation(
    title = "Overlap of TIV (SDY1276) TLS signatures between Females and Males",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25))
  )

ggsave(fs::path(figure_path, "rise_signature_overlap_tiv_sex.pdf"),
       p_combined, width = 15, height = 7, dpi = 300)

rm(list = ls())
