library(ggVennDiagram)
library(ggplot2)
library(patchwork)
library(kableExtra)

# ---- Paths & data ----
signature_path <- fs::path("output", "results", "application")
figure_path    <- fs::path("output", "figures", "application")

sig_rVSV    <- readRDS(fs::path(signature_path, "TLS_rVSV_prevac.rds"))
sig_Ad26MVA <- readRDS(fs::path(signature_path, "TLS_Ad26MVA_prevac.rds"))
sig_TIV     <- readRDS(fs::path(signature_path, "TLS_TIV_SDY1276.rds"))
BTM         <- readRDS(fs::path("data", "BTM_processed.rds"))

# ---- Build gene- and geneset-level signature lists ----
gene_lists <- list(rVSV = sig_rVSV, Ad26MVA = sig_Ad26MVA, TIV = sig_TIV)

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
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
  plot.margin = margin(10, 20, 10, 20)  # extra right/left margin for spacing
)
venn_style <- list(
  scale_fill_gradient(low = "#F4F9FF", high = "#4A7FBF"),
  scale_color_manual(values = rep("#2C3E50", 3)),
  venn_theme
)

make_venn <- function(lst, subtitle) {
  ggVennDiagram(lst, label = "count", label_alpha = 0, edge_size = 0.8) +
    venn_style +
    labs(title = subtitle)
}

p1 <- make_venn(gene_lists, "A) Gene-level")
p2 <- make_venn(geneset_lists, "B) Geneset-level")

# ---- Combine side by side, with blank spacer + one shared, larger title ----
p_combined <- (p1 | plot_spacer() | p2) +
  plot_layout(widths = c(1, 0.08, 1)) +  # thin spacer column between plots
  plot_annotation(
    title = "Overlap of TLS Signatures Across Vaccines",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25))
  )

ggsave(fs::path(figure_path, "Figure3-overlap_combined.pdf"),
       p_combined, width = 13, height = 6, dpi = 300)

# ---- Intersections (common + pairwise) at both levels ----
get_intersections <- function(lst) {
  list(
    common_all   = Reduce(intersect, lst),
    AB = intersect(lst[[1]], lst[[2]]),
    AC = intersect(lst[[1]], lst[[3]]),
    BC = intersect(lst[[2]], lst[[3]])
  )
}
gene_intersections    <- get_intersections(gene_lists)
geneset_intersections <- get_intersections(geneset_lists)

# Table to detail all intersections at the geneset level

# ---- Extract exact regions (exclusive intersections) for geneset-level lists ----
venn_data <- process_data(Venn(geneset_lists))
region_df <- venn_data$regionData[, c("name", "item")]

# keep only multi-set intersections (exclude single-vaccine-only regions)
region_df <- region_df[grepl("/", region_df$name), ]

# format: "/"-separated set names -> " ∩ " ; gene/geneset vectors -> comma-separated string
region_df$Intersection <- gsub("/", " $\\\\cap$ ", region_df$name)
region_df$Genesets     <- vapply(region_df$item, paste, collapse = ", ", FUN.VALUE = character(1))

table_out <- region_df[, c("Intersection", "Genesets")]

kbl(table_out, format = "latex", booktabs = TRUE, longtable = TRUE,
    caption = "Genesets shared across vaccine platform signatures",
    col.names = c("Intersection", "Genesets"),
    escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
  column_spec(2, width = "10cm")

# rm(list = ls())
