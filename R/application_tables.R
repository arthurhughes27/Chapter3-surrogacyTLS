# Shared (non-analysis) helper for building the "significant markers"
# supplementary LaTeX tables from a completed RISE screening run
# (analysis/application/RISE-*.R), used by
# analysis/supplementary/significant_markers_table_*.R.

#' Build a LaTeX table (.tex) of significant screening-stage markers:
#' marker name (capitalised), delta with a 90% Wald confidence interval,
#' adjusted p-value, and normalised screening weight, sorted by adjusted
#' p-value (most to least significant).
#'
#' The 90% CI is computed here as delta +/- z * sd (z = qnorm(0.95)),
#' since rise.screen()'s own ci_lower/ci_upper columns are the one-sided
#' bound used internally for the non-inferiority test at whatever alpha
#' screening was run with, not a descriptive two-sided interval.
#'
#' @param screening_metrics data.frame - rise.screen.res[["screening.metrics"]]
#'   (must have columns marker, delta, sd, p_adjusted).
#' @param screening_weights data.frame - rise.screen.res[["screening.weights"]]
#'   (columns marker, weight; weight is assumed already normalised, i.e.
#'   from rise.screen()'s default normalise.weights = TRUE).
#' @param significant_markers character vector of significant marker
#'   names, e.g. rise.screen.res[["significant.markers"]].
#' @param out_path .tex output file path.
#' @param caption,label passed to xtable::xtable().
#' @param ci_level confidence level for the reported interval (default 0.90).
build_significant_markers_table <- function(screening_metrics, screening_weights,
                                             significant_markers, out_path,
                                             caption, label, ci_level = 0.90) {
  z <- qnorm(1 - (1 - ci_level) / 2)

  table_df <- screening_metrics %>%
    dplyr::filter(marker %in% significant_markers) %>%
    dplyr::left_join(screening_weights %>% dplyr::select(marker, weight), by = "marker") %>%
    dplyr::arrange(p_adjusted) %>%
    dplyr::transmute(
      Marker = toupper(marker),
      `Delta (90% CI)` = sprintf("%.2f (%.2f, %.2f)", delta, delta - z * sd, delta + z * sd),
      `Adjusted p-value` = formatC(p_adjusted, format = "e", digits = 2),
      `Normalised weight` = sprintf("%.2f", weight)
    )

  fs::dir_create(fs::path_dir(out_path))
  print(
    xtable::xtable(table_df, caption = caption, label = label),
    include.rownames = FALSE,
    booktabs = TRUE,
    file = out_path
  )

  table_df
}
