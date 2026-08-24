# Shared (non-analysis) helper functions used by the RISE application
# scripts in analysis/application/. These wrap the rise.screen() ->
# rise.evaluate() -> combined-plot pipeline that was previously duplicated
# (with only inputs and labels differing) across RISE-Ad26MVA.R,
# RISE-rVSV.R and RISE-SDY1276.R.

#' Run RISE screening followed by evaluation, combine the two diagnostic
#' plots into a single labelled figure, and save it to disk.
#'
#' @return A list with the rise.screen() result, the rise.evaluate() result,
#'   the selected markers/weights, and the combined ggplot grob.
run_rise_pipeline <- function(yone_train, yzero_train, sone_train, szero_train,
                               yone_test, yzero_test, sone_test, szero_test,
                               screen_label, eval_label, figure_path,
                               screen_power = 0.8, screen_paired = FALSE,
                               eval_power = 0.9, eval_paired = TRUE,
                               n_cores = 6, screen_plot_topN = 20,
                               alpha = 0.05, p_correction = "BH",
                               alternative = "two.sided",
                               weight_mode = "diff.epsilon",
                               plot_width = 17, plot_height = 6, plot_dpi = 300) {

  rise.screen.res <- rise.screen(
    yone = yone_train, yzero = yzero_train,
    sone = sone_train, szero = szero_train,
    alpha = alpha, power.want.s = screen_power,
    p.correction = p_correction, alternative = alternative,
    paired = screen_paired, weight.mode = weight_mode,
    n.cores = n_cores, screen.plot.topN = screen_plot_topN
  )

  markers <- rise.screen.res[["significant.markers"]]
  weights <- rise.screen.res[["screening.weights"]]
  p1 <- rise.screen.res$plot$screen.plot

  rise.eval.res <- rise.evaluate(
    yone = yone_test, yzero = yzero_test,
    sone = sone_test, szero = szero_test,
    alpha = alpha, power.want.s = eval_power,
    p.correction = p_correction, alternative = alternative,
    paired = eval_paired, n.cores = n_cores,
    markers = markers, screening.weights = weights
  )

  p2 <- rise.eval.res[["gamma.s.plot"]]

  p1_labeled <- p1 +
    labs(title = screen_label) +
    theme(
      plot.title = element_text(face = "bold", size = 25, hjust = 0),
      plot.title.position = "panel",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )

  p2_labeled <- p2 +
    labs(title = eval_label) +
    theme(
      plot.title = element_text(face = "bold", size = 25, hjust = 0),
      plot.title.position = "panel",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )

  # Thin vertical divider between panels
  divider <- ggplot() +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1),
                 colour = "white", linewidth = 0.6) +
    theme_void()

  g <- egg::ggarrange(
    p1_labeled, divider, p2_labeled,
    ncol = 3,
    widths = c(1, 0.2, 1),
    draw = FALSE
  )

  grid.newpage()
  grid.draw(g)

  ggsave(figure_path, g, width = plot_width, height = plot_height, dpi = plot_dpi)

  list(
    screen = rise.screen.res,
    evaluate = rise.eval.res,
    markers = markers,
    weights = weights,
    plot = g
  )
}
