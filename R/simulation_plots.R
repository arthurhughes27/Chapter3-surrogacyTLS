# Shared (non-analysis) plotting helpers for the RISE simulation figures in
# analysis/simulation/ and analysis/supplementary/. Each function is
# parameterised by data/labels/output path so the same plot code is reused
# across the main-text and supplementary "twin" figures (identical designs,
# differing only by data-generation mode) instead of being copy-pasted.

#' Save a plot as both .pdf and .eps at `out_path` (without extension),
#' creating the destination directory if it does not already exist.
save_plot_pdf_eps <- function(plot, out_path, width, height, dpi = 800) {
  fs::dir_create(fs::path_dir(out_path))
  ggsave(paste0(out_path, ".pdf"), plot, units = "cm", width = width, height = height, dpi = dpi)
  ggsave(paste0(out_path, ".eps"), plot, units = "cm", width = width, height = height, dpi = dpi)
}

#' Boxplot of a false positive rate metric across a grid (e.g. sample size),
#' with the nominal alpha level marked as a dashed reference line.
plot_fpr_boxplot <- function(df, x_var, x_lab, out_path,
                              alpha_line = 0.05, y_max = 0.15,
                              width = 40, height = 25) {
  p <- ggplot(df, aes(x = as.factor(.data[[x_var]]), y = fpr)) +
    geom_boxplot() +
    labs(x = x_lab, y = "False Positive Rate") +
    ylim(0, y_max) +
    geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1.5) +
    theme_minimal(base_size = 35) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          axis.title = element_text(size = 35),
          legend.position = "none")

  save_plot_pdf_eps(p, out_path, width = width, height = height)
  p
}

#' Violin plot of a false positive rate metric across a grid (e.g.
#' inter-predictor correlation), with the nominal alpha level marked.
plot_fpr_violin <- function(df, x_var, x_lab, out_path,
                             alpha_line = 0.05, width = 40, height = 20) {
  p <- ggplot(df, aes(x = as.factor(.data[[x_var]]), y = fpr, fill = .data[[x_var]])) +
    geom_violin() +
    labs(x = x_lab, y = "False Positive Rate") +
    geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6) +
    theme_minimal(base_size = 35) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none",
          axis.title.x = element_text(size = 55)) +
    scale_fill_gradient2(low = "white", mid = "lightcoral", high = "red", midpoint = 0.25)

  save_plot_pdf_eps(p, out_path, width = width, height = height)
  p
}

#' Side-by-side empirical power (TPR) and false discovery proportion (FDR)
#' line plots against average surrogate strength, one line per sample
#' size, sharing a single legend panel.
plot_power_fdr_vs_strength <- function(df, out_path, alpha_line = 0.05,
                                        width = 35, height = 18) {
  fdr_plot <- ggplot(df, aes(x = avg_us, y = avg_fdr, color = as.factor(n), linetype = as.factor(n))) +
    geom_line(linewidth = 1.3) +
    geom_point() +
    ylab("False Discovery Proportion") +
    xlab(expression(bar(U[S]))) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(title = "Sample size"),
           linetype = guide_legend(title = "Sample size")) +
    geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6) +
    theme_minimal(base_size = 24) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  tpr_plot <- ggplot(df, aes(x = avg_us, y = avg_tpr, color = as.factor(n), linetype = as.factor(n))) +
    geom_line(linewidth = 1.3) +
    geom_point() +
    ylab("Empirical power") +
    xlab(expression(bar(U[S]))) +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(title = "Sample size"),
           linetype = guide_legend(title = "Sample size")) +
    theme_minimal(base_size = 24) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  legend <- cowplot::get_legend(fdr_plot + theme(legend.key.width = unit(2, "cm")))

  combined <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(
      tpr_plot + ylim(0, 1) + theme(legend.position = "none"),
      fdr_plot + ylim(0, 1) + theme(legend.position = "none"),
      legend,
      ncol = 3, widths = c(1, 1, 0.3)
    )
  )

  save_plot_pdf_eps(combined, out_path, width = width, height = height)
  combined
}

#' Power/FDR vs surrogate strength, one line per multiple-testing
#' correction method (unadjusted, Bonferroni, BH, BY), for a fixed sample
#' size.
plot_power_fdr_by_correction <- function(df, out_path, alpha_line = 0.05,
                                          width = 41, height = 18) {
  correction_labels <- c(unadjusted = "Unadjusted", bonf = "Bonferroni", bh = "B-H", by = "B-Y")
  correction_linetypes <- c(unadjusted = "solid", bonf = "dashed", bh = "dotted", by = "dotdash")
  palette <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
  correction_colors <- setNames(palette, names(correction_labels))

  fdr_long <- df %>%
    tidyr::pivot_longer(cols = c(avg_fdr_unadjusted, avg_fdr_bonf, avg_fdr_bh, avg_fdr_by),
                         names_to = "correction", values_to = "value") %>%
    dplyr::mutate(correction = sub("^avg_fdr_", "", correction)) %>%
    dplyr::select(avg_us, correction, value)

  tpr_long <- df %>%
    tidyr::pivot_longer(cols = c(avg_tpr_unadjusted, avg_tpr_bonf, avg_tpr_bh, avg_tpr_by),
                         names_to = "correction", values_to = "value") %>%
    dplyr::mutate(correction = sub("^avg_tpr_", "", correction)) %>%
    dplyr::select(avg_us, correction, value)

  build_plot <- function(long_df, y_lab) {
    ggplot(long_df, aes(x = avg_us, y = value, linetype = correction, color = correction)) +
      geom_line(linewidth = 1.3) +
      labs(x = expression(bar(U[S])), y = y_lab, linetype = "Correction", color = "Correction") +
      scale_linetype_manual(values = correction_linetypes, labels = correction_labels) +
      scale_color_manual(values = correction_colors, labels = correction_labels) +
      theme_minimal(base_size = 24) +
      theme(plot.background = element_rect(fill = "white", color = "white"))
  }

  fdr_plot <- build_plot(fdr_long, "False Discovery Proportion") +
    geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6)
  tpr_plot <- build_plot(tpr_long, "Empirical Power")

  legend <- cowplot::get_legend(tpr_plot + theme(legend.key.width = unit(2, "cm")))

  combined <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(
      tpr_plot + ylim(0, 1) + theme(legend.position = "none"),
      fdr_plot + ylim(0, 1) + theme(legend.position = "none"),
      legend,
      ncol = 3, widths = c(1, 1, 0.3)
    )
  )

  save_plot_pdf_eps(combined, out_path, width = width, height = height)
  combined
}

#' Side-by-side empirical power and FDR violin plots against
#' inter-predictor correlation.
plot_power_fdr_violin_by_corr <- function(df, out_path, width = 35, height = 18) {
  tpr_plot <- ggplot(df, aes(x = as.factor(correlation), y = tpr, fill = correlation)) +
    geom_violin() +
    ylim(0, 1) +
    ylab("Empirical Power") +
    xlab(expression(sigma[corr])) +
    theme_minimal(base_size = 24) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none",
          axis.title.x = element_text(size = 35)) +
    scale_fill_gradient2(low = "white", mid = "lightcoral", high = "red", midpoint = 0.25)

  fdr_plot <- ggplot(df, aes(x = as.factor(correlation), y = fdr, fill = correlation)) +
    geom_violin() +
    ylim(0, 1) +
    ylab("False Discovery Proportion") +
    xlab(expression(sigma[corr])) +
    geom_hline(yintercept = 0.05, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6) +
    theme_minimal(base_size = 24) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none",
          axis.title.x = element_text(size = 35)) +
    scale_fill_gradient2(low = "white", mid = "lightcoral", high = "red", midpoint = 0.25)

  combined <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(tpr_plot, fdr_plot, ncol = 2, widths = c(1, 1))
  )

  save_plot_pdf_eps(combined, out_path, width = width, height = height)
  combined
}

#' Boxplot of evaluation-stage p-values against the proportion of invalid
#' surrogates making up gamma_S.
plot_gamma_pvalue_boxplot <- function(df_long, out_path, alpha_line = 0.05,
                                       width = 35, height = 20) {
  p <- ggplot(df_long, aes(x = as.factor(prop_invalid), y = p_value)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "Proportion of invalid surrogates", y = "P-value") +
    geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6) +
    theme_minimal(base_size = 26) +
    theme(plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none")

  save_plot_pdf_eps(p, out_path, width = width, height = height)
  p
}

#' Two gamma_S p-value boxplots (e.g. differing numbers of candidates)
#' side by side, each labelled a)/b).
plot_gamma_pvalue_boxplot_pair <- function(df_long_a, df_long_b, out_path,
                                            alpha_line = 0.05, width = 65, height = 22) {
  build_plot <- function(df_long, hide_y_title = FALSE) {
    p <- ggplot(df_long, aes(x = as.factor(prop_invalid), y = p_value)) +
      geom_boxplot(outlier.shape = NA) +
      labs(x = "Proportion of invalid surrogates", y = "P-value") +
      geom_hline(yintercept = alpha_line, colour = "maroon", linetype = "longdash", linewidth = 1, alpha = 0.6) +
      theme_minimal(base_size = 26) +
      theme(plot.background = element_rect(fill = "white", color = "white"),
            legend.position = "none")
    if (hide_y_title) p <- p + theme(axis.title.y = element_blank())
    p
  }

  panel_a <- build_plot(df_long_a)
  panel_b <- build_plot(df_long_b, hide_y_title = TRUE)

  label_a <- cowplot::ggdraw() + cowplot::draw_label("a)", fontface = "bold", size = 35, x = 0, hjust = 0)
  label_b <- cowplot::ggdraw() + cowplot::draw_label("b)", fontface = "bold", size = 35, x = 0, hjust = 0)
  label_row <- cowplot::plot_grid(label_a, label_b, ncol = 2)

  combined <- cowplot::plot_grid(
    label_row,
    cowplot::plot_grid(panel_a, panel_b, ncol = 2),
    ncol = 1, rel_heights = c(0.1, 1)
  )

  save_plot_pdf_eps(combined, out_path, width = width, height = height)
  combined
}

#' Histogram of p-values under the null hypothesis.
plot_null_pvalue_histogram <- function(df, out_path, width = 25, height = 7) {
  p <- ggplot(df, aes(x = p)) +
    geom_histogram(binwidth = 0.03, fill = "lightblue", color = "black", aes(y = after_stat(density))) +
    labs(x = "P-value", y = "Density") +
    xlim(0, 1) +
    theme_minimal(base_size = 20) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  save_plot_pdf_eps(p, out_path, width = width, height = height)
  p
}
