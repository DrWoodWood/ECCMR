#' Syndrome and global test of pleiotropy
#'
#' Computes the syndrome `s = H %*% beta_Y` of the ECC-MR framework and the
#' associated global test of pleiotropy. Under the null hypothesis of no
#' horizontal pleiotropy, the syndrome equals `H %*% epsilon` and the
#' quadratic form `T_syn = s' Omega^+ s` (Moore-Penrose pseudoinverse) follows
#' a chi-squared distribution with `rank(H)` degrees of freedom. Standardized
#' check statistics localize the corrupted regions, providing an ECC-native
#' analogue of the MR-PRESSO global and outlier tests that remains valid for
#' correlated instruments.
#'
#' @param beta_Y Numeric vector of length p; outcome association estimates.
#' @param se_Y Numeric vector of length p; standard errors of `beta_Y`.
#' @param H Sparse m x p parity-check matrix, from [build_parity_check()].
#' @param check_block Optional integer vector of length m assigning each check
#'   to an LD block; derived from `H` when omitted.
#' @param sigma Optional p x p error covariance matrix (use when the exposure
#'   and outcome GWAS share participants); defaults to `diag(se_Y^2)`.
#' @param tol_eig Numeric; relative eigenvalue tolerance for the pseudoinverse.
#'
#' @return A list with components
#'   \describe{
#'     \item{statistic}{The syndrome statistic `T_syn`.}
#'     \item{df}{Degrees of freedom, equal to the numerical rank of H.}
#'     \item{p_value}{P-value of the global pleiotropy test.}
#'     \item{syndrome}{The syndrome vector s.}
#'     \item{z}{Standardized check statistics `s_k / sqrt(Omega_kk)`.}
#'   }
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' pc <- build_parity_check(sim$beta_X, sim$blocks)
#' syndrome_test(sim$beta_Y, sim$se_Y, pc$H, pc$check_block)
#'
#' @export
syndrome_test <- function(beta_Y, se_Y, H, check_block = NULL,
                          sigma = NULL, tol_eig = 1e-10) {
  m <- nrow(H)
  if (m == 0L) {
    return(list(statistic = 0, df = 0L, p_value = NA_real_,
                syndrome = numeric(0L), z = numeric(0L)))
  }
  if (is.null(check_block)) check_block <- .check_blocks(H)
  Sigma <- if (is.null(sigma)) Matrix::Diagonal(x = se_Y^2) else sigma
  s <- as.numeric(H %*% beta_Y)
  Omega <- H %*% Matrix::tcrossprod(Sigma, H)

  stat <- 0
  df <- 0L
  idx_split <- split(seq_len(m), check_block)
  for (idx in idx_split) {
    Ob <- as.matrix(Omega[idx, idx, drop = FALSE])
    ee <- eigen(Ob, symmetric = TRUE)
    d <- ee$values
    keep <- d > tol_eig * max(d, 0)
    if (!any(keep)) next
    proj <- as.numeric(crossprod(ee$vectors[, keep, drop = FALSE], s[idx]))
    stat <- stat + sum(proj^2 / d[keep])
    df <- df + sum(keep)
  }
  p_value <- if (df > 0L) stats::pchisq(stat, df = df, lower.tail = FALSE) else NA_real_
  z <- s / sqrt(diag(as.matrix(Omega)))
  list(statistic = stat, df = df, p_value = p_value, syndrome = s, z = z)
}
