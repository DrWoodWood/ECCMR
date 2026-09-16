## Operating-boundary scan: when does ECC-MR gain over conventional MR?

#' Scan the operating boundary of ECC-MR across pleiotropy regimes
#'
#' Runs a grid of simulation regimes varying (i) the fraction of pleiotropic
#' instruments (`pis`), (ii) the mean direct effect (`directionals`; 0 gives
#' balanced pleiotropy) and (iii) the per-variant pleiotropy scale
#' (`tau_alphas`), and reports bias / Monte-Carlo SD / RMSE of ECC-MR and
#' conventional summary-data MR methods in every regime. This answers the
#' question "when does modelling LD with parity checks actually help":
#' ECC-MR is expected to dominate whenever pleiotropic effects are
#' individually detectable (roughly |alpha| >~ 2-3 se_Y) or balanced, and to
#' fall back to IVW-like behaviour when pleiotropy is both frequent and
#' weaker than the noise floor.
#'
#' Only point estimates are computed (no bootstrap), so the scan is fast.
#'
#' @param pis Numeric vector; fractions of pleiotropic instruments.
#' @param directionals Numeric vector; mean direct effects. Use 0 for
#'   balanced pleiotropy and a positive value for directional pleiotropy.
#' @param tau_alphas Numeric vector; SD of the direct effects.
#' @param nsim Integer; Monte-Carlo replicates per grid cell.
#' @param n_snps,block_size,theta,se_y Passed to [simulate_eccmr()].
#' @param methods Character vector; subset of `"ECC-MR"`, `"IVW"`,
#'   `"MR-Egger"`, `"Weighted median"`, `"MR-Lasso"`.
#' @param prune Logical; apply conventional methods to one LD-pruned
#'   representative per block (default `TRUE`), while ECC-MR always uses all
#'   instruments.
#' @param nlambda Tuning controls passed to [eccmr()].
#' @param cores Integer; parallel workers over replicates (Unix only).
#' @param seed Integer; base seed (each cell/replicate gets a derived seed).
#' @param verbose Logical; print per-cell progress.
#'
#' @return An object of class `eccmr_boundary`: a data frame with one row per
#'   regime x method, columns `directional`, `pi`, `tau_alpha`, `method`,
#'   `bias`, `mc_sd`, `rmse`, `n_reps`.
#'
#' @examples
#' bd <- run_boundary_scan(pis = c(0.1, 0.3), directionals = c(0.05, 0),
#'                         tau_alphas = c(0.02, 0.08), nsim = 4)
#' print(bd)
#'
#' @export
run_boundary_scan <- function(pis = c(0.1, 0.2, 0.3),
                              directionals = c(0.05, 0),
                              tau_alphas = c(0.02, 0.04, 0.08),
                              nsim = 10L, n_snps = 300L, block_size = 10L,
                              theta = 0.3, se_y = 0.015,
                              methods = c("ECC-MR", "IVW", "MR-Egger",
                                          "Weighted median", "MR-Lasso"),
                              prune = TRUE, nlambda = 12L,
                              cores = 1L, seed = 1L, verbose = TRUE) {
  methods <- match.arg(methods, several.ok = TRUE,
                       choices = c("ECC-MR", "IVW", "MR-Egger",
                                   "Weighted median", "MR-Lasso"))

  one_rep <- function(rep_seed, pi, direc, tau) {
    sim <- simulate_eccmr(n_snps = n_snps, block_size = block_size,
                          theta = theta, se_y = se_y, scenario = "A",
                          pi_pleio = pi, directional = direc,
                          tau_alpha = tau, seed = rep_seed)
    bx <- sim$beta_X; by <- sim$beta_Y; sy <- sim$se_Y
    if (prune) {
      kp <- .prune_ld(bx, sy, sim$blocks)
      bx_p <- bx[kp]; by_p <- by[kp]; sy_p <- sy[kp]
    } else {
      bx_p <- bx; by_p <- by; sy_p <- sy
    }
    out <- numeric(length(methods)); names(out) <- methods
    for (m in methods) {
      out[m] <- switch(m,
        "ECC-MR" = eccmr(bx, by, sy, blocks = sim$blocks,
                         nlambda = nlambda,
                         n_boot = 0L, seed = rep_seed)$theta,
        "IVW" = mr_ivw(bx_p, by_p, sy_p)$theta,
        "MR-Egger" = mr_egger(bx_p, by_p, sy_p)$theta,
        "Weighted median" = mr_weighted_median(bx_p, by_p, sy_p,
                                               n_boot = 0L)$theta,
        "MR-Lasso" = mr_lasso(bx_p, by_p, sy_p, n_boot = 0L)$theta)
    }
    out
  }

  apply_fun <- if (cores > 1L) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply
  }

  grid <- expand.grid(directional = directionals, pi = pis,
                      tau_alpha = tau_alphas)
  rows <- list()
  for (g in seq_len(nrow(grid))) {
    d <- grid$directional[g]; pi <- grid$pi[g]; tau <- grid$tau_alpha[g]
    cell_seed <- seed + g * 10000L
    ests <- apply_fun(seq_len(nsim), function(r)
      one_rep(cell_seed + r, pi, d, tau))
    ests <- do.call(rbind, ests)
    for (m in methods) {
      v <- ests[, m]
      rows[[length(rows) + 1L]] <- data.frame(
        directional = d, pi = pi, tau_alpha = tau, method = m,
        bias = mean(v) - theta, mc_sd = stats::sd(v),
        rmse = sqrt(mean((v - theta)^2)), n_reps = nsim,
        stringsAsFactors = FALSE)
    }
    if (verbose) {
      ecc_rmse <- sqrt(mean((ests[, "ECC-MR"] - theta)^2))
      message(sprintf("directional=%.2f pi=%.1f tau=%.2f done: ECC-MR rmse=%.4f",
                      d, pi, tau, ecc_rmse))
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("eccmr_boundary", "data.frame")
  out
}

#' @export
print.eccmr_boundary <- function(x, ...) {
  df <- as.data.frame(x)
  cells <- unique(df[, c("directional", "pi", "tau_alpha")])
  cat("ECC-MR operating-boundary scan\n")
  cat("==============================\n")
  for (i in seq_len(nrow(cells))) {
    sub <- merge(df, cells[i, , drop = FALSE])
    cat(sprintf("\ndirectional=%.2f  pi=%.1f  tau_alpha=%.2f  (n=%d)\n",
                cells$directional[i], cells$pi[i], cells$tau_alpha[i],
                sub$n_reps[1]))
    tab <- sub[, c("method", "bias", "mc_sd", "rmse")]
    best_other <- min(tab$rmse[tab$method != "ECC-MR"])
    tab$rmse <- sprintf("%.4f%s", tab$rmse,
                        ifelse(tab$method == "ECC-MR" &
                               tab$rmse <= best_other, " *", ""))
    tab$bias <- sprintf("%+.4f", tab$bias)
    tab$mc_sd <- sprintf("%.4f", tab$mc_sd)
    print(tab, row.names = FALSE)
    ecc <- sub$rmse[sub$method == "ECC-MR"]
    cat(sprintf("ECC-MR / best conventional RMSE ratio: %.2f\n",
                ecc / best_other))
  }
  cat("\n(* = lowest RMSE in the regime)\n")
  invisible(x)
}
