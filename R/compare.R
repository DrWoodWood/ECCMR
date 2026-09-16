## Benchmarking ECC-MR against mainstream MR methods on simulated data.

#' Compare ECC-MR with mainstream MR methods on one data set
#'
#' Applies ECC-MR and a set of conventional summary-data MR methods
#' (IVW, MR-Egger, weighted median, MR-Lasso) to the same data set. By
#' default the conventional methods are applied to one LD-pruned
#' representative instrument per LD block (standard clumping practice),
#' whereas ECC-MR uses all instruments and models the LD structure directly
#' -- this contrast is precisely where ECC-MR gains its information
#' advantage.
#'
#' @param beta_X,beta_Y Numeric vectors of exposure and outcome association
#'   estimates.
#' @param se_Y Numeric vector of standard errors of `beta_Y`.
#' @param blocks Optional list of integer vectors of LD blocks; derived from
#'   `ld_mat` when omitted.
#' @param ld_mat Optional p x p matrix of pairwise LD correlations.
#' @param r2_thresh Numeric; r^2 threshold used when deriving blocks from
#'   `ld_mat`.
#' @param methods Character vector; any of `"ECC-MR"`, `"IVW"`, `"MR-Egger"`,
#'   `"Weighted median"`, `"MR-Lasso"`.
#' @param prune Logical; apply conventional methods to one LD-pruned
#'   representative per block (default `TRUE`). If `FALSE`, they use all
#'   instruments, violating their independence assumption.
#' @param n_boot Integer; bootstrap replicates for the ECC-MR standard error.
#' @param median_boot,lasso_boot Integer; bootstrap replicates for the
#'   weighted-median and MR-Lasso standard errors.
#' @param nlambda,cv_folds Tuning controls passed to [eccmr()].
#' @param seed Optional seed for reproducibility.
#' @param verbose Logical; print progress messages.
#'
#' @return A data frame with one row per method: `method`, `theta`, `se`,
#'   `ci_lo`, `ci_hi`, and the number of instruments used (`n_instruments`).
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
#'                       seed = 42)
#' compare_mr_methods(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
#'                    n_boot = 50, seed = 1)
#'
#' @export
compare_mr_methods <- function(beta_X, beta_Y, se_Y, blocks = NULL,
                               ld_mat = NULL, r2_thresh = 0.3,
                               methods = c("ECC-MR", "IVW", "MR-Egger",
                                           "Weighted median", "MR-Lasso"),
                               prune = TRUE, n_boot = 100L,
                               median_boot = 200L, lasso_boot = 200L,
                               nlambda = 10L,
                               cv_folds = 5L, seed = NULL, cores = 1L,
                               verbose = FALSE) {
  if (is.null(blocks)) {
    if (is.null(ld_mat)) {
      stop("Provide either `blocks` or `ld_mat` to define the LD structure.")
    }
    blocks <- partition_ld_blocks(ld_mat, r2_thresh = r2_thresh)
  }
  if (!is.null(seed)) set.seed(seed)
  p_all <- length(beta_X)
  idx <- seq_len(p_all)
  if (prune) idx <- sort(.prune_ld(beta_X, se_Y, blocks))

  rows <- list()
  for (meth in methods) {
    if (verbose) message("Fitting ", meth, "...")
    res <- switch(
      meth,
      "ECC-MR" = {
        f <- eccmr(beta_X, beta_Y, se_Y, blocks = blocks, n_boot = n_boot,
                   nlambda = nlambda,
                   cv_folds = cv_folds, cores = cores, verbose = verbose)
        list(theta = f$theta, se = f$se, ci_lo = f$ci[1L], ci_hi = f$ci[2L],
             n_instruments = p_all)
      },
      "IVW" = {
        f <- mr_ivw(beta_X[idx], beta_Y[idx], se_Y[idx])
        c(f[c("theta", "se", "ci_lo", "ci_hi")], n_instruments = length(idx))
      },
      "MR-Egger" = {
        f <- mr_egger(beta_X[idx], beta_Y[idx], se_Y[idx])
        c(f[c("theta", "se", "ci_lo", "ci_hi")], n_instruments = length(idx))
      },
      "Weighted median" = {
        f <- mr_weighted_median(beta_X[idx], beta_Y[idx], se_Y[idx],
                                n_boot = median_boot)
        c(f[c("theta", "se", "ci_lo", "ci_hi")], n_instruments = length(idx))
      },
      "MR-Lasso" = {
        f <- mr_lasso(beta_X[idx], beta_Y[idx], se_Y[idx],
                      n_boot = lasso_boot)
        c(f[c("theta", "se", "ci_lo", "ci_hi")], n_instruments = length(idx))
      },
      stop("Unknown method: '", meth, "'. Available: 'ECC-MR', 'IVW', ",
           "'MR-Egger', 'Weighted median', 'MR-Lasso'.")
    )
    rows[[meth]] <- data.frame(
      method = meth, theta = res$theta, se = res$se,
      ci_lo = res$ci_lo, ci_hi = res$ci_hi,
      n_instruments = res$n_instruments,
      row.names = NULL
    )
  }
  do.call(rbind, rows)
}

#' Simulation benchmark of ECC-MR against mainstream MR methods
#'
#' Repeatedly simulates GWAS summary statistics with [simulate_eccmr()],
#' applies ECC-MR and conventional MR methods to each replicate via
#' [compare_mr_methods()], and aggregates bias, Monte-Carlo variability,
#' RMSE, confidence-interval coverage and rejection rates. This reproduces
#' the simulation study of the ECC-MR manuscript and demonstrates where
#' ECC-MR outperforms the alternatives.
#'
#' @param nsim Integer; number of simulation replicates (default 20; use 50+
#'   for stable Monte-Carlo error).
#' @param scenario Character; pleiotropy scenario passed to
#'   [simulate_eccmr()] (`"A"`, `"B"` or `"C"`).
#' @param methods Character vector of methods, see [compare_mr_methods()].
#' @param prune Logical; LD-prune instruments for conventional methods
#'   (default `TRUE`).
#' @param n_snps,block_size,theta,... Simulation parameters passed to
#'   [simulate_eccmr()] (e.g. `pi_pleio`, `directional`, `delta`).
#' @param n_boot,nlambda,cv_folds ECC-MR controls for each
#'   replicate; defaults are deliberately light to keep the benchmark fast.
#' @param cores Integer; number of parallel workers (Unix only, via
#'   `parallel::mclapply`); `1` (default) runs serially.
#' @param seed Integer base seed; replicate `k` uses `seed + k - 1`.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `eccmr_sim`: a list with
#'   \describe{
#'     \item{summary}{Data frame with one row per method: mean estimate,
#'       bias, empirical SD, RMSE, mean reported SE, 95% CI coverage and
#'       rejection rate.}
#'     \item{raw}{Data frame with one row per method per replicate.}
#'     \item{theta, scenario, nsim}{Benchmark meta-information.}
#'   }
#'
#' @examples
#' \dontrun{
#' bench <- run_simulation(nsim = 50, scenario = "A", directional = 0.05,
#'                         seed = 1)
#' print(bench)
#' plot(bench)
#' }
#'
#' @export
run_simulation <- function(nsim = 20L, scenario = "A",
                           methods = c("ECC-MR", "IVW", "MR-Egger",
                                       "Weighted median", "MR-Lasso"),
                           prune = TRUE, n_snps = 300L, block_size = 10L,
                           theta = 0.3, ..., n_boot = 100L, nlambda = 10L,
                           cv_folds = 5L,
                           cores = 1L, seed = 1L, verbose = TRUE) {
  sim_args <- list(...)
  one_rep <- function(k) {
    sim <- do.call(simulate_eccmr,
                   c(list(n_snps = n_snps, block_size = block_size,
                          theta = theta, scenario = scenario,
                          seed = seed + k - 1), sim_args))
    res <- compare_mr_methods(sim$beta_X, sim$beta_Y, sim$se_Y,
                              blocks = sim$blocks, methods = methods,
                              prune = prune, n_boot = n_boot,
                              nlambda = nlambda,
                              cv_folds = cv_folds, verbose = FALSE)
    res$rep <- k
    res
  }
  if (cores > 1L) {
    reps <- parallel::mclapply(seq_len(nsim), one_rep, mc.cores = cores)
  } else {
    reps <- lapply(seq_len(nsim), function(k) {
      if (verbose) message("Replicate ", k, " / ", nsim)
      one_rep(k)
    })
  }
  raw <- do.call(rbind, reps)

  summ <- do.call(rbind, lapply(split(raw, raw$method), function(d) {
    has_se <- is.finite(d$se)
    data.frame(
      method = d$method[1L],
      estimate = mean(d$theta),
      bias = mean(d$theta) - theta,
      mc_sd = stats::sd(d$theta),
      rmse = sqrt(mean((d$theta - theta)^2)),
      mean_se = if (any(has_se)) mean(d$se[has_se]) else NA_real_,
      coverage = if (any(has_se))
        mean(d$ci_lo[has_se] <= theta & d$ci_hi[has_se] >= theta) else NA_real_,
      rejection = if (any(has_se))
        mean(abs(d$theta[has_se] / d$se[has_se]) > stats::qnorm(0.975)) else NA_real_,
      row.names = NULL
    )
  }))
  ## keep the user-specified method order
  summ <- summ[match(methods, summ$method), ]
  summ <- summ[!is.na(summ$method), ]

  out <- list(summary = summ, raw = raw, theta = theta,
              scenario = scenario, nsim = nsim, prune = prune)
  class(out) <- "eccmr_sim"
  out
}

#' @export
print.eccmr_sim <- function(x, digits = 4, ...) {
  cat("ECC-MR simulation benchmark\n")
  cat(sprintf("  scenario %s | true theta = %g | %d replicates | %s\n",
              x$scenario, x$theta, x$nsim,
              if (x$prune) "conventional methods LD-pruned" else "no pruning"))
  tab <- x$summary
  disp <- data.frame(
    method = tab$method,
    estimate = round(tab$estimate, digits),
    bias = round(tab$bias, digits),
    sd = round(tab$mc_sd, digits),
    rmse = round(tab$rmse, digits),
    coverage = round(tab$coverage, 3L),
    rejection = round(tab$rejection, 3L)
  )
  best <- which.min(tab$rmse)
  disp$method <- as.character(disp$method)
  disp$method[best] <- paste0(disp$method[best], " *")
  print(disp, row.names = FALSE)
  cat("(* lowest RMSE)\n")
  invisible(x)
}

#' @export
plot.eccmr_sim <- function(x, ...) {
  s <- x$summary
  cols <- ifelse(s$method == "ECC-MR", "#B2182B", "grey65")
  old <- graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
  on.exit(graphics::par(old))
  ## bias panel with Monte-Carlo error bars
  mc_err <- s$mc_sd / sqrt(x$nsim)
  rng <- range(c(s$bias - mc_err, s$bias + mc_err, 0))
  pad <- 0.15 * diff(rng)
  bp <- graphics::barplot(s$bias, col = cols, las = 2, names.arg = s$method,
                          ylim = rng + c(-pad, pad),
                          main = sprintf("Bias (scenario %s)", x$scenario),
                          ylab = "Bias")
  graphics::arrows(bp, s$bias - mc_err, bp, s$bias + mc_err,
                   angle = 90, code = 3, length = 0.05)
  graphics::abline(h = 0, col = "grey30")
  ## RMSE panel
  graphics::barplot(s$rmse, col = cols, las = 2, names.arg = s$method,
                    main = sprintf("RMSE (scenario %s)", x$scenario),
                    ylab = "RMSE")
  invisible(x)
}
