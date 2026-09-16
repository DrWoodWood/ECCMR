#' eccmr: Error-Correcting Code Inspired Robust Mendelian Randomization
#'
#' ECC-MR reframes horizontal pleiotropy as sparse, correctable corruption of
#' SNP-level causal signals and leverages linkage disequilibrium (LD) among
#' instruments as structured redundancy. A parity-check matrix built from the
#' LD graph annihilates the causal signal, so the syndrome of the observed
#' outcome associations isolates pleiotropic errors; sparse syndrome decoding
#' within a convex optimization framework jointly estimates the causal effect
#' and the pleiotropic effects without discarding instruments.
#'
#' The main entry point is [eccmr()]. Lower-level building blocks are
#' [partition_ld_blocks()], [build_parity_check()], [syndrome_test()] and
#' [ecc_fit()]. [simulate_eccmr()] generates synthetic summary statistics for
#' examples and testing; [compare_mr_methods()] and [run_simulation()]
#' benchmark ECC-MR against mainstream MR methods (IVW, MR-Egger, weighted
#' median, MR-Lasso) on simulated data.
#'
#' @references
#' Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
#' Framework for Robust Mendelian Randomization. Manuscript.
#'
#' @keywords internal
#' @name eccmr-package
#' @aliases eccmr-package
"_PACKAGE"
