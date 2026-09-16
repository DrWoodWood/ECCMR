#' ECC-MR: Error-correcting code inspired robust Mendelian randomization
#'
#' Fits the ECC-MR model to GWAS summary statistics. Horizontal pleiotropy is
#' modelled as sparse corruption of SNP-level causal signals; a parity-check
#' matrix built from the LD structure annihilates the causal signal
#' (`H beta_X = 0`), so the syndrome `H beta_Y` isolates pleiotropic errors.
#' The causal effect theta and the pleiotropic effects alpha are estimated
#' jointly by convex optimization; tuning parameters are selected by a BIC
#' criterion along the penalty path (or optionally by LD-block
#' cross-validation), the selected support is refitted without the L1
#' penalty to remove shrinkage bias, and standard errors are obtained from
#' an LD-block bootstrap. Unlike MR-PRESSO-style approaches, pleiotropic
#' instruments are corrected and retained rather than discarded, and
#' correlated instruments can be used directly without LD pruning.
#'
#' @param beta_X,beta_Y Numeric vectors of length p; exposure and outcome
#'   association estimates for the genetic instruments.
#' @param se_Y Numeric vector of length p; standard errors of `beta_Y`.
#' @param blocks Optional list of integer vectors giving the LD block of each
#'   instrument (e.g. external annotations). Required unless `ld_mat` is
#'   given.
#' @param ld_mat Optional p x p matrix of pairwise LD correlations, used to
#'   derive blocks with [partition_ld_blocks()] and to threshold the checks.
#' @param r2_thresh Numeric; r^2 threshold for edges/checks (default 0.3).
#' @param lambda,gamma Optional tuning parameters. `lambda` = L1 decoding
#'   penalty (BIC/CV-selected when `NULL`). `gamma` is deprecated since
#'   v0.3.0: the syndrome-consistency penalty is no longer part of the core
#'   objective (it defaults to 0; a positive value is still honoured as an
#'   optional regulariser). When `NULL` (the
#'   default), they are selected according to `tuning`.
#' @param tuning Character; tuning-parameter selection strategy. `"bic"`
#'   (default) minimizes a BIC criterion along the lambda path,
#'   which rewards actually fitting the pleiotropic effects; `"cv"` uses
#'   K-fold cross-validation over LD blocks with the one-standard-error
#'   rule.
#' @param nlambda Integer; size of the lambda grid (default 20).
#' @param lambda_min_ratio Numeric; smallest lambda as a fraction of the
#'   largest (default 0.01).
#'   `tr(W) / tr(H' Omega^+ H)` searched during tuning.
#' @param cv_folds Integer; number of cross-validation folds (only used when
#'   `tuning = "cv"`; default 10).
#' @param n_boot Integer; number of LD-block bootstrap replicates (default
#'   500; set to 0 to skip the bootstrap).
#' @param snp_names Optional character vector of instrument names (e.g.
#'   rsIDs), used to label the detected pleiotropic variants.
#' @param cores Integer; number of parallel workers for the LD-block
#'   bootstrap and the tuning path (default 1 = sequential). Uses a
#'   PSOCK cluster, so it works
#'   on Windows; results are identical for any `cores` at fixed `seed`.
#' @param seed Optional integer seed for reproducibility of the
#'   cross-validation fold assignment and the bootstrap.
#' @param max_iter,tol Solver controls passed to [ecc_fit()].
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `eccmr`: a list with components
#'   \describe{
#'     \item{theta}{Estimated causal effect from the relaxed (debiased)
#'       refit; `theta_l1` gives the L1-penalized estimate before debiasing.}
#'     \item{se, ci}{Bootstrap standard error and 95% percentile confidence
#'       interval (when `n_boot > 0`).}
#'     \item{alpha}{Estimated pleiotropic effect of each instrument (relaxed
#'       refit; `alpha_l1` gives the penalized version used for support
#'       selection).}
#'     \item{pleiotropic}{Indices (and names, when supplied) of instruments
#'       with nonzero estimated pleiotropic effect.}
#'     \item{syndrome}{Result of [syndrome_test()]: global pleiotropy test
#'       and standardized check statistics.}
#'     \item{lambda, gamma}{Selected tuning parameters.}
#'     \item{tune_table}{Data frame of BIC (or cross-validation) scores
#'       along the tuning path.}
#'     \item{boot_theta}{Bootstrap replicates of the causal estimate.}
#'     \item{n_instruments, n_checks}{Problem dimensions.}
#'     \item{fit}{The final relaxed [ecc_fit()] result (iterations,
#'       convergence).}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
#' Framework for Robust Mendelian Randomization. Manuscript.
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
#'                       seed = 42)
#' fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
#'              nlambda = 10, n_boot = 100, seed = 1)
#' print(fit)
#' summary(fit)
#'
#' @export
eccmr <- function(beta_X, beta_Y, se_Y, blocks = NULL, ld_mat = NULL,
                  r2_thresh = 0.3, lambda = NULL, gamma = NULL,
                  tuning = c("bic", "cv"),
                  nlambda = 20L, lambda_min_ratio = 0.01,
                  cv_folds = 10L, n_boot = 500L, snp_names = NULL,
                  seed = NULL, max_iter = 200L, tol = 1e-8,
                  cores = 1L, verbose = FALSE) {
  tuning <- match.arg(tuning)
  cl <- match.call()
  p <- length(beta_X)
  stopifnot(length(beta_Y) == p, length(se_Y) == p,
            all(se_Y > 0), all(is.finite(beta_X)))
  if (!is.null(seed)) set.seed(seed)

  ## ---- 1. LD blocks and parity-check matrix ---------------------------
  if (is.null(blocks)) {
    if (is.null(ld_mat)) {
      stop("Provide either `blocks` (LD-block membership of each instrument) ",
           "or `ld_mat` (pairwise LD correlations) so that the LD structure ",
           "can be derived.")
    }
    blocks <- partition_ld_blocks(ld_mat, r2_thresh = r2_thresh)
  }
  pc <- build_parity_check(beta_X, blocks,
                           r2_mat = if (is.null(ld_mat)) NULL else as.matrix(ld_mat)^2,
                           tau = r2_thresh)
  H <- pc$H
  m <- nrow(H)
  if (verbose) message("ECC-MR: ", p, " instruments, ", length(blocks),
                       " LD blocks, ", m, " parity checks.")

  ## ---- 2. Syndrome and global pleiotropy test --------------------------
  syn <- syndrome_test(beta_Y, se_Y, H, check_block = pc$check_block)

  ## ---- 3. Tuning-parameter selection -----------------------------------
  tune_table <- NULL
  ## v0.3.0: the syndrome-consistency penalty is removed from the core
  ## objective (extensive simulation showed it brings no systematic gain;
  ## the syndrome is retained as a *diagnostic* via syndrome_test()).
  ## gamma = 0 reduces the objective to the weighted sparse-pleiotropy
  ## decoder; a user-supplied gamma > 0 is still honoured as an optional
  ## regulariser.
  if (is.null(gamma)) gamma <- 0
  c_gammas <- 0
  if (is.null(lambda)) {
    g0 <- if (gamma > 0) .gamma0(se_Y, H, pc$check_block) else 0
    lambdas <- .lambda_grid(beta_X, beta_Y, se_Y, H, pc$check_block,
                            nlambda = nlambda,
                            lambda_min_ratio = lambda_min_ratio)
    if (tuning == "bic") {
      bic_res <- .bic_eccmr(beta_X, beta_Y, se_Y, H, pc$check_block,
                            lambdas, c_gammas,
                            max_iter = max_iter, tol = tol, cores = cores)
      tune_table <- bic_res$table
      best <- which.min(tune_table$bic)
      lambda <- tune_table$lambda[best]
    } else if (length(blocks) >= 5L) {
      tune_table <- .cv_eccmr(beta_X, beta_Y, se_Y, H, pc$edges,
                              pc$check_block, blocks, lambdas, c_gammas,
                              cv_folds = cv_folds, max_iter = max_iter,
                              tol = tol, cores = cores)
      best <- which.min(tune_table$cv_error)
      cg_sel <- tune_table$c_gamma[best]
      if (is.null(lambda)) {
        ## one-standard-error rule within the selected c_gamma: the largest
        ## (sparsest) lambda whose CV error is within one SE of the minimum
        sub <- tune_table[tune_table$c_gamma == cg_sel, ]
        thr <- min(sub$cv_error) + sub$cv_se[which.min(sub$cv_error)]
        ok <- sub$cv_error <= thr
        lambda <- max(sub$lambda[ok])
      }
    } else {
      ## too few blocks for cross-validation: universal-threshold default
      lambda <- 2 * sqrt(log(p)) * stats::median(1 / se_Y)
      if (verbose) {
        message("Too few LD blocks for cross-validation; ",
                "using universal-threshold defaults.")
      }
    }
    if (verbose) message("Selected lambda = ", signif(lambda, 4),
                         " (", tuning, "); gamma = 0 (syndrome penalty ",
                         "removed in v0.3.0).")
  }

  ## ---- 4. Final fit: L1 decoding + relaxed (debiased) refit --------------
  fit <- ecc_fit(beta_X, beta_Y, se_Y, H, lambda = lambda, gamma = gamma,
                 check_block = pc$check_block, max_iter = max_iter, tol = tol)
  fit_relaxed <- .relaxed_refit(beta_X, beta_Y, se_Y, H, pc$check_block,
                                alpha_hat = fit$alpha, gamma = gamma,
                                max_iter = max_iter, tol = tol)

  ## ---- 5. LD-block bootstrap --------------------------------------------
  boot_theta <- NULL
  se <- NA_real_
  ci <- c(NA_real_, NA_real_)
  if (n_boot > 0L) {
    if (verbose) message("Running ", n_boot, " LD-block bootstrap replicates...")
    boot_theta <- .boot_eccmr(beta_X, beta_Y, se_Y, blocks,
                              lambda = lambda, gamma = gamma,
                              B = n_boot, max_iter = max_iter, tol = tol,
                              cores = cores)
    se <- stats::sd(boot_theta)
    ci <- as.numeric(stats::quantile(boot_theta, c(0.025, 0.975),
                                     names = FALSE, na.rm = TRUE))
  }

  ## ---- 6. Collect results (relaxed refit is the primary estimate) --------
  pleio_idx <- which(abs(fit$alpha) > 0)
  pleio <- if (!is.null(snp_names) && length(pleio_idx) > 0L) {
    stats::setNames(fit_relaxed$alpha[pleio_idx], snp_names[pleio_idx])
  } else {
    fit_relaxed$alpha[pleio_idx]
  }
  out <- list(
    theta = fit_relaxed$theta, theta_l1 = fit$theta,
    se = se, ci = ci,
    alpha = fit_relaxed$alpha, alpha_l1 = fit$alpha,
    pleiotropic_index = pleio_idx,
    pleiotropic = pleio,
    syndrome = syn,
    lambda = lambda, gamma = gamma, tuning = tuning, tune_table = tune_table,
    boot_theta = boot_theta,
    n_instruments = p, n_checks = m,
    fit = list(iter = fit_relaxed$iter, converged = fit_relaxed$converged,
               objective = fit_relaxed$objective),
    call = cl
  )
  class(out) <- "eccmr"
  out
}
