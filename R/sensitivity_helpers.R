# Shared (non-analysis) helpers for the epsilon sensitivity analyses in
# analysis/supplementary/sensitivity_analyses/. These evaluate the effect
# of fixing the screening-stage non-inferiority margin (epsilon) at
# several values, instead of choosing it adaptively via power.want.s as
# analysis/application/RISE-*.R (via R/rise_helpers.R's
# run_rise_pipeline()) do. Evaluation-stage epsilon is left adaptive
# (via eval_power) in every case - only the screening-stage epsilon
# varies - and every other setting (p-value correction, weighting mode,
# alternative, paired/independent) matches run_rise_pipeline()'s
# defaults, so this is a minimal, targeted variant of the main
# application pipeline rather than a separately-tuned analysis.
#
# The epsilon grid (0.05, 0.1, ..., 0.35) matches the original RISE
# thesis/paper's own sensitivity analysis for SDY1276 (females) exactly
# (see RISE_application_ver2.Rmd, Table S2), applied uniformly across
# all four applications here for consistency.

#' For a single application's train/test split, sweep a grid of fixed
#' screening-stage epsilon values and, for each, report the resulting
#' evaluation-stage composite gamma_Gamma's characteristics: number of
#' significant markers, delta with a 90% CI, its standard deviation, and
#' the (unadjusted) evaluation p-value. Grid points yielding zero
#' significant markers (screening epsilon too strict) are reported with
#' NA evaluation columns rather than silently dropped.
run_rise_epsilon_sensitivity <- function(yone_train, yzero_train, sone_train, szero_train,
                                          yone_test, yzero_test, sone_test, szero_test,
                                          epsilon_grid = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                                          screen_paired = FALSE, eval_power = 0.9, eval_paired = TRUE,
                                          n_cores = 6, alpha = 0.05, p_correction = "BH",
                                          alternative = "two.sided", weight_mode = "diff.epsilon",
                                          ci_level = 0.90) {
  z <- qnorm(1 - (1 - ci_level) / 2)

  rows <- lapply(epsilon_grid, function(eps) {
    screen_res <- rise.screen(
      yone = yone_train, yzero = yzero_train, sone = sone_train, szero = szero_train,
      alpha = alpha, epsilon = eps, p.correction = p_correction, alternative = alternative,
      paired = screen_paired, weight.mode = weight_mode, n.cores = n_cores
    )

    markers <- screen_res[["significant.markers"]]
    n_markers <- length(markers)

    if (n_markers == 0) {
      return(data.frame(
        epsilon_screen = eps, n_predictors = 0L,
        delta_ci = NA_character_, sigma_delta = NA_real_, p_value = NA_character_
      ))
    }

    # return.all.evaluate/return.plot.evaluate must stay at their
    # defaults (TRUE) - the installed SurrogateRank::rise.evaluate()
    # references individual.metrics/gamma.s.plot unconditionally in its
    # return() list, so passing FALSE for either throws (see
    # R/simulation_helpers.R for the same issue elsewhere).
    eval_res <- rise.evaluate(
      yone = yone_test, yzero = yzero_test, sone = sone_test, szero = szero_test,
      alpha = alpha, power.want.s = eval_power, p.correction = p_correction,
      alternative = alternative, paired = eval_paired, n.cores = n_cores,
      markers = markers, screening.weights = screen_res[["screening.weights"]]
    )

    ev <- eval_res[["gamma.s.evaluate"]]
    delta <- unname(ev["delta"])
    sigma_delta <- unname(ev["sd"])
    p_unadjusted <- unname(ev["p_unadjusted"])

    data.frame(
      epsilon_screen = eps, n_predictors = n_markers,
      delta_ci = sprintf("%.3f (%.3f, %.3f)", delta, delta - z * sigma_delta, delta + z * sigma_delta),
      sigma_delta = round(sigma_delta, 3),
      p_value = formatC(p_unadjusted, format = "e", digits = 1)
    )
  })

  do.call(rbind, rows)
}

#' Write an epsilon sensitivity table (from run_rise_epsilon_sensitivity())
#' to a .tex file, with the exact LaTeX column headers from the thesis
#' (surrogate notation S -> Gamma): epsilon (screening), No. of genes in
#' gamma_Gamma, delta_gammaGamma (90% CI), p-value. The standard
#' deviation column from sensitivity_df is dropped (deemed irrelevant
#' for the table, though it remains in the saved sensitivity_df/.rds).
build_epsilon_sensitivity_table <- function(sensitivity_df, out_path, caption, label) {
  table_df <- sensitivity_df[, c("epsilon_screen", "n_predictors", "delta_ci", "p_value")]
  colnames(table_df) <- c(
    "$\\boldsymbol{\\epsilon}$ \\textbf{(screening)}",
    "\\textbf{No. of genes in} $\\boldsymbol{\\gamma_{\\Gamma}}$",
    "$\\boldsymbol{\\delta_{\\gamma_{\\Gamma}}}$ \\textbf{(90\\% C.I.)}",
    "\\textbf{p-value}"
  )

  fs::dir_create(fs::path_dir(out_path))
  print(
    xtable::xtable(table_df, caption = caption, label = label),
    include.rownames = FALSE,
    booktabs = TRUE,
    sanitize.colnames.function = function(x) x,
    file = out_path
  )

  table_df
}
