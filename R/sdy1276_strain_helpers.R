# Shared (non-analysis) helper for the per-strain SDY1276 (TIV)
# supplementary applications in
# analysis/supplementary/supplementary_applications/RISE-SDY1276-*.R
# (Brisbane10, Brisbane59, Florida), which analyse each vaccine strain's
# antibody response separately instead of the mean-across-strains
# response (ab_p_28/ab_p_0) used by the main SDY1276 application.

#' Find the exact df_merged_all column name for a specific SDY1276 (TIV)
#' antibody strain and timepoint.
#'
#' Strain-level columns are produced during preprocessing
#' (preprocessing_is2_immuneresponse.R) as
#' ab_p_<timepoint>_<clean strain name>, where the strain name suffix
#' comes from janitor::make_clean_names() applied to the raw virus name
#' in the source data - its exact text (e.g. a trailing "2007" vs
#' "2008") isn't independently confirmed from code alone. Matching on a
#' short, stable keyword substring (e.g. "florida", "brisbane_10")
#' avoids depending on that exact suffix, and this errors loudly if
#' zero or more than one column matches, rather than silently selecting
#' the wrong column or an empty one.
find_strain_column <- function(cols, timepoint, keyword) {
  pattern <- paste0("^ab_p_", timepoint, "_.*", keyword)
  matches <- grep(pattern, cols, value = TRUE, ignore.case = TRUE)
  if (length(matches) != 1) {
    stop(sprintf(
      "Expected exactly one column matching keyword '%s' at timepoint %s, found %d: %s",
      keyword, timepoint, length(matches),
      if (length(matches) == 0) "(none)" else paste(matches, collapse = ", ")
    ))
  }
  matches
}
