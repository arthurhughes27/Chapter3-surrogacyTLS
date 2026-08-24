# Shared (non-analysis) helper functions for the RISE simulation studies in
# analysis/simulation/ and analysis/supplementary/.
#
# These replace the bespoke data-generation, ground-truth and
# screening/evaluation logic that used to live inline, and repeated, across
# every chunk of the old simulation_combined.Rmd:
#   - every input is an explicit function argument (no reliance on globals
#     set by an earlier, unrelated chunk);
#   - the gamma_S screen/evaluate pipeline (simulate_gamma_pvalues) now
#     calls SurrogateRank::rise.screen()/rise.evaluate() directly, rather
#     than the old bespoke rise_screen()/rise_evaluate() functions with a
#     different (now unavailable) API;
#   - the n_sim replicate loop (independent draws) is parallelised via
#     pbmcapply instead of a plain sequential for loop.

#' Stop with a clear, single error message if any parallel replicate
#' failed, instead of letting the failure surface many calls later as a
#' confusing "index out of bounds" or similar error when the (silently
#' missing) result is used downstream.
check_replicates <- function(reps) {
  failed <- vapply(reps, function(r) inherits(r, "try-error"), logical(1))
  if (any(failed)) {
    first_error <- reps[[which(failed)[1]]]
    stop(sprintf(
      "%d of %d parallel simulation replicates failed. First error:\n%s",
      sum(failed), length(reps), conditionMessage(attr(first_error, "condition"))
    ))
  }
  reps
}

#' Simulate primary outcomes and candidate surrogates for one replicate.
#'
#' @param mode "simple": valid surrogates are the primary response plus
#'   Gaussian noise, invalid surrogates are Gaussian. "complex": valid
#'   surrogates are the cubed primary response plus Gaussian noise, invalid
#'   surrogates are Exponential.
gen_data <- function(n1, n0, p, prop_valid, valid_sigma, corr,
                      y1_mean, y1_sd, y0_mean, y0_sd,
                      mode = c("simple", "complex")) {
  mode <- match.arg(mode)

  p_valid <- as.numeric(prop_valid * p)
  p_invalid <- as.numeric((1 - prop_valid) * p)

  y1 <- rnorm(n1, y1_mean, y1_sd)
  y0 <- rnorm(n0, y0_mean, y0_sd)

  if (mode == "simple") {
    mm <- runif(p_invalid, min = 0.5, max = 2.5)
    ss <- runif(p_invalid, min = 0.5, max = 2)
    Sigma_invalid <- matrix(corr, nrow = p_invalid, ncol = p_invalid)
    diag(Sigma_invalid) <- ss

    if (prop_valid != 0) {
      if (p_valid == 1) {
        Sigma_valid <- c(valid_sigma)
        s1.valid <- y1 + MASS::mvrnorm(n = n1, mu = 0, Sigma = Sigma_valid)
        s0.valid <- y0 + MASS::mvrnorm(n = n0, mu = 0, Sigma = Sigma_valid)
      } else {
        Sigma_valid <- matrix(corr * valid_sigma, nrow = p_valid, ncol = p_valid)
        diag(Sigma_valid) <- rep(valid_sigma, p_valid)
        s1.valid <- matrix(y1, nrow = n1, ncol = p_valid, byrow = TRUE) +
          MASS::mvrnorm(n = n1, mu = rep(0, p_valid), Sigma = Sigma_valid)
        s0.valid <- matrix(y0, nrow = n0, ncol = p_valid, byrow = TRUE) +
          MASS::mvrnorm(n = n0, mu = rep(0, p_valid), Sigma = Sigma_valid)
      }

      if (prop_valid != 1) {
        if (p_invalid > 1) {
          s1.invalid <- MASS::mvrnorm(n = n1, mu = mm, Sigma = Sigma_invalid)
          s0.invalid <- MASS::mvrnorm(n = n0, mu = mm, Sigma = Sigma_invalid)
        } else {
          mm1 <- runif(1, min = 0.5, max = 2.5)
          ss1 <- runif(1, min = 0.5, max = 2)
          s1.invalid <- rnorm(n = n1, mean = mm1, sd = ss1)
          s0.invalid <- rnorm(n = n0, mean = mm1, sd = ss1)
        }
        s1 <- cbind(s1.valid, s1.invalid)
        s0 <- cbind(s0.valid, s0.invalid)
      } else {
        s1 <- s1.valid
        s0 <- s0.valid
      }
    } else {
      s1 <- MASS::mvrnorm(n = n1, mu = mm, Sigma = Sigma_invalid)
      s0 <- MASS::mvrnorm(n = n0, mu = mm, Sigma = Sigma_invalid)
    }
  } else { # complex
    lambda <- runif(p_invalid, min = 0.5, max = 2.5)

    if (prop_valid != 0) {
      Sigma_valid <- matrix(corr * valid_sigma, nrow = p_valid, ncol = p_valid)
      diag(Sigma_valid) <- rep(valid_sigma, p_valid)
      s1.valid <- matrix(y1^3, nrow = n1, ncol = p_valid, byrow = TRUE) +
        MASS::mvrnorm(n = n1, mu = rep(0, p_valid), Sigma = Sigma_valid)
      s0.valid <- matrix(y0^3, nrow = n0, ncol = p_valid, byrow = TRUE) +
        MASS::mvrnorm(n = n0, mu = rep(0, p_valid), Sigma = Sigma_valid)

      if (prop_valid != 1) {
        s0.invalid <- sapply(lambda, function(rate) rexp(n0, rate))
        s1.invalid <- sapply(lambda, function(rate) rexp(n1, rate))
        s1 <- cbind(s1.valid, s1.invalid)
        s0 <- cbind(s0.valid, s0.invalid)
      } else {
        s1 <- s1.valid
        s0 <- s0.valid
      }
    } else {
      s0 <- sapply(lambda, function(rate) rexp(n0, rate))
      s1 <- sapply(lambda, function(rate) rexp(n1, rate))
    }
  }

  hyp <- c(rep("null false", p_valid), rep("null true", p_invalid))
  list(y1 = y1, y0 = y0, s1 = s1, s0 = s0, hyp = hyp)
}

#' True U_Y and per-surrogate U_S values via a large simulated sample, used
#' to translate a target average surrogate strength into a valid_sigma
#' value for the smaller-sample simulations.
calc_truth <- function(p, prop_valid, valid_sigma, corr,
                        y1_mean, y1_sd, y0_mean, y0_sd,
                        mode = c("simple", "complex"), n_truth = 10000) {
  mode <- match.arg(mode)
  dd <- gen_data(n1 = n_truth, n0 = n_truth, p = p, prop_valid = prop_valid,
                 valid_sigma = valid_sigma, corr = corr,
                 y1_mean = y1_mean, y1_sd = y1_sd, y0_mean = y0_mean, y0_sd = y0_sd,
                 mode = mode)

  uy <- (n_truth * n_truth)^(-1) * wilcox.test(dd$y1, dd$y0)$statistic

  valid_idx <- which(dd$hyp == "null false")
  us <- numeric(length(valid_idx))
  for (j in valid_idx) {
    us[j] <- (n_truth * n_truth)^(-1) * wilcox.test(dd$s1[, j], dd$s0[, j])$statistic
  }

  list(uy_true = uy, us_true = us,
       delta_true = mean(uy - us), delta_true_sd = sd(uy - us))
}

#' Run n_sim replicates of per-candidate univariate surrogate testing (via
#' SurrogateRank::test.surrogate), and summarise FPR/FDR/TPR/PPV under no
#' correction, Bonferroni, Benjamini-Hochberg and Benjamini-Yekutieli
#' correction. Replicates are independent, so are parallelised via
#' pbmcapply when n_cores > 1.
simulate_screening_metrics <- function(n1, n0, p, prop_valid, n_sim, valid_sigma, corr,
                                        y1_mean, y1_sd, y0_mean, y0_sd,
                                        mode = c("simple", "complex"),
                                        alpha = 0.05, n_cores = 1) {
  mode <- match.arg(mode)

  one_replicate <- function(k) {
    data <- gen_data(n1, n0, p, prop_valid, valid_sigma, corr,
                      y1_mean, y1_sd, y0_mean, y0_sd, mode = mode)

    u_y_estimated <- SurrogateRank::test.surrogate(
      yone = data$y1, yzero = data$y0, sone = data$y1, szero = data$y0, epsilon = 0.1
    )$u.y
    eps <- max(0, u_y_estimated - 0.5)

    p_unadjusted <- vapply(seq_len(p), function(j) {
      ss.test <- SurrogateRank::test.surrogate(
        yone = data$y1, yzero = data$y0,
        sone = data$s1[, j], szero = data$s0[, j], epsilon = eps
      )
      pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta)
    }, numeric(1))

    p_bonf <- p.adjust(p_unadjusted, method = "bonferroni")
    p_bh <- p.adjust(p_unadjusted, method = "BH")
    p_by <- p.adjust(p_unadjusted, method = "BY")

    classify <- function(p_values) {
      TP <- sum(data$hyp == "null false" & p_values < alpha)
      FP <- sum(data$hyp == "null true" & p_values < alpha)
      TN <- sum(data$hyp == "null true" & p_values >= alpha)
      FN <- sum(data$hyp == "null false" & p_values >= alpha)

      if (prop_valid == 0) {
        list(fpr = FP / (FP + TN), fdr = 0, tpr = 0, ppv = 0)
      } else {
        list(fpr = FP / (FP + TN), fdr = FP / (TP + FP),
             tpr = TP / (TP + FN), ppv = TP / (TP + FP))
      }
    }

    list(
      p_unadjusted = p_unadjusted, p_bonf = p_bonf, p_bh = p_bh, p_by = p_by,
      unadjusted = classify(p_unadjusted), bonf = classify(p_bonf),
      bh = classify(p_bh), by = classify(p_by)
    )
  }

  reps <- check_replicates(pbmcapply::pbmclapply(seq_len(n_sim), one_replicate, mc.cores = n_cores))

  summarise_method <- function(name) {
    metric_mean <- function(metric) {
      v <- vapply(reps, function(r) r[[name]][[metric]], numeric(1))
      v[is.nan(v)] <- 0
      mean(v)
    }
    list(avg_fpr = metric_mean("fpr"), avg_fdr = metric_mean("fdr"),
         avg_tpr = metric_mean("tpr"), avg_ppv = metric_mean("ppv"))
  }

  p_values <- list(
    p_unadjusted = do.call(rbind, lapply(reps, `[[`, "p_unadjusted")),
    p_bonf = do.call(rbind, lapply(reps, `[[`, "p_bonf")),
    p_bh = do.call(rbind, lapply(reps, `[[`, "p_bh")),
    p_by = do.call(rbind, lapply(reps, `[[`, "p_by"))
  )

  metrics <- list(
    metrics_unadjusted = summarise_method("unadjusted"),
    metrics_bonf = summarise_method("bonf"),
    metrics_bh = summarise_method("bh"),
    metrics_by = summarise_method("by")
  )

  # Per-replicate metrics (one row per simulation x correction method), for
  # figures that show the distribution across replicates (boxplots/violins)
  # rather than only the mean.
  per_replicate <- do.call(rbind, lapply(c("unadjusted", "bonf", "bh", "by"), function(name) {
    data.frame(
      replicate = seq_len(n_sim),
      correction = name,
      fpr = vapply(reps, function(r) r[[name]]$fpr, numeric(1)),
      fdr = vapply(reps, function(r) r[[name]]$fdr, numeric(1)),
      tpr = vapply(reps, function(r) r[[name]]$tpr, numeric(1)),
      ppv = vapply(reps, function(r) r[[name]]$ppv, numeric(1))
    )
  }))
  per_replicate[is.nan(per_replicate$fpr), "fpr"] <- 0
  per_replicate[is.nan(per_replicate$fdr), "fdr"] <- 0
  per_replicate[is.nan(per_replicate$tpr), "tpr"] <- 0
  per_replicate[is.nan(per_replicate$ppv), "ppv"] <- 0

  list(p_values = p_values, metrics = metrics, per_replicate = per_replicate)
}

#' For a single (prop_invalid) design point, generate n_sim replicate
#' gamma_S composite markers built from ALL p candidates at a fixed
#' proportion of invalid surrogates, and return the evaluation-stage
#' unadjusted p-value from SurrogateRank::rise.evaluate() for each
#' replicate.
#'
#' Mirrors the original bespoke pipeline (screen for per-marker weights,
#' then evaluate the composite built from every candidate, not just the
#' significant ones) using the current package: rise.screen() is called
#' with weight.mode = "inverse.delta" (== the original's
#' abs(1/delta) weight), normalise.weights = FALSE and
#' return.all.weights = TRUE so every candidate (not only significant
#' ones) gets a usable weight, and that full weight table is passed to
#' rise.evaluate() together with markers = all candidate names.
simulate_gamma_pvalues <- function(n1, n0, p, prop_invalid, valid_sigma, corr,
                                    y1_mean, y1_sd, y0_mean, y0_sd,
                                    mode = c("simple", "complex"),
                                    n_sim = 500, n_cores = 1) {
  mode <- match.arg(mode)
  prop_valid <- 1 - prop_invalid

  one_replicate <- function(k) {
    data <- gen_data(n1, n0, p, prop_valid, valid_sigma, corr,
                      y1_mean, y1_sd, y0_mean, y0_sd, mode = mode)

    marker_names <- paste0("marker", seq_len(p))
    sone <- data$s1
    szero <- data$s0
    colnames(sone) <- marker_names
    colnames(szero) <- marker_names

    u_y_estimated <- SurrogateRank::test.surrogate(
      yone = data$y1, yzero = data$y0, sone = data$y1, szero = data$y0, epsilon = 0.1
    )$u.y
    eps <- max(0, u_y_estimated - 0.5)

    screen_res <- rise.screen(
      yone = data$y1, yzero = data$y0, sone = sone, szero = szero,
      alpha = 0.05, epsilon = eps, p.correction = "none", n.cores = 1,
      weight.mode = "inverse.delta", normalise.weights = FALSE,
      return.all.weights = TRUE
    )

    eval_res <- rise.evaluate(
      yone = data$y1, yzero = data$y0, sone = sone, szero = szero,
      alpha = 0.05, epsilon = eps, p.correction = "none", n.cores = 1,
      markers = marker_names, screening.weights = screen_res[["screening.weights"]],
      return.plot.evaluate = FALSE, return.all.evaluate = FALSE
    )

    unname(eval_res[["gamma.s.evaluate"]]["p_unadjusted"])
  }

  reps <- check_replicates(pbmcapply::pbmclapply(seq_len(n_sim), one_replicate, mc.cores = n_cores))
  unlist(reps)
}

#' Read a cached .rds file if it exists, otherwise evaluate `expr`, save
#' the result to `path`, and return it. Set `force = TRUE` to ignore an
#' existing cache and recompute (e.g. after changing simulation
#' parameters).
cache_rds <- function(path, expr, force = FALSE) {
  if (!force && file.exists(path)) {
    return(readRDS(path))
  }
  result <- expr
  fs::dir_create(fs::path_dir(path))
  saveRDS(result, path)
  result
}
