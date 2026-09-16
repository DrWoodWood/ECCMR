#' Partition instruments into LD blocks
#'
#' Builds approximately independent LD blocks from a pairwise LD (correlation)
#' matrix by thresholding squared correlations and taking connected components
#' of the resulting graph.
#'
#' @param ld_mat A p x p matrix (or sparse `Matrix`) of pairwise LD
#'   correlations r between instruments, estimated from a reference panel of
#'   matched ancestry.
#' @param r2_thresh Numeric; edges connect instruments with r^2 > `r2_thresh`.
#'
#' @return A list of integer vectors; each vector contains the instrument
#'   indices of one LD block (singletons allowed).
#'
#' @examples
#' ld <- diag(6); ld[1, 2] <- ld[2, 1] <- 0.9; ld[3, 4] <- ld[4, 3] <- 0.8
#' partition_ld_blocks(ld, r2_thresh = 0.3)
#'
#' @export
partition_ld_blocks <- function(ld_mat, r2_thresh = 0.3) {
  r2 <- as.matrix(ld_mat)^2
  p <- nrow(r2)
  edges <- which(r2 > r2_thresh & upper.tri(r2), arr.ind = TRUE)
  lab <- .union_find_components(edges, p)
  blocks <- split(seq_len(p), lab)
  blocks <- lapply(blocks, sort)
  names(blocks) <- NULL
  blocks
}

#' Construct the LD parity-check matrix
#'
#' Builds the parity-check matrix H of the ECC-MR framework. Each retained
#' within-block edge (i, j) defines one parity check (one row of H) with
#' entries `H[k, i] = beta_X[j]` and `H[k, j] = -beta_X[i]`, so that
#' `H %*% beta_X = 0` by construction: the causal signal lies in the null
#' space (code space) of H, whereas pleiotropic effects generate a nonzero
#' syndrome.
#'
#' @param beta_X Numeric vector of length p; exposure association estimates.
#' @param blocks A list of integer vectors giving the LD blocks, e.g. from
#'   [partition_ld_blocks()] or from external LD-block annotations
#'   (Berisa & Pickrell, 2016).
#' @param r2_mat Optional p x p matrix of squared LD correlations; when given,
#'   only within-block pairs with r^2 > `tau` generate a check. When `NULL`,
#'   all within-block pairs are used.
#' @param tau Numeric; r^2 threshold applied when `r2_mat` is supplied.
#'
#' @return A list with components
#'   \describe{
#'     \item{H}{Sparse (`dgCMatrix`) m x p parity-check matrix.}
#'     \item{edges}{m x 2 integer matrix of instrument index pairs (checks).}
#'     \item{check_block}{Integer vector of length m with the LD-block label
#'       of each check.}
#'   }
#'
#' @references
#' Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
#' Framework for Robust Mendelian Randomization. Manuscript.
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, seed = 1)
#' pc <- build_parity_check(sim$beta_X, sim$blocks)
#' max(abs(as.numeric(pc$H %*% sim$beta_X)))  # exactly 0
#'
#' @export
build_parity_check <- function(beta_X, blocks, r2_mat = NULL, tau = 0.3) {
  p <- length(beta_X)
  edge_list <- list()
  block_id <- list()
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    if (length(blk) < 2L) next
    pairs <- t(utils::combn(blk, 2L))
    if (!is.null(r2_mat)) {
      keep <- r2_mat[pairs] > tau
      pairs <- pairs[keep, , drop = FALSE]
    }
    if (nrow(pairs) == 0L) next
    edge_list[[length(edge_list) + 1L]] <- pairs
    block_id[[length(block_id) + 1L]] <- rep(b, nrow(pairs))
  }
  if (length(edge_list) == 0L) {
    return(list(
      H = Matrix::sparseMatrix(i = integer(0L), j = integer(0L),
                               x = numeric(0L), dims = c(0L, p)),
      edges = matrix(integer(0L), ncol = 2L),
      check_block = integer(0L)
    ))
  }
  edges <- do.call(rbind, edge_list)
  check_block <- unlist(block_id, use.names = FALSE)
  m <- nrow(edges)
  H <- Matrix::sparseMatrix(
    i = c(seq_len(m), seq_len(m)),
    j = c(edges[, 1L], edges[, 2L]),
    x = c(beta_X[edges[, 2L]], -beta_X[edges[, 1L]]),
    dims = c(m, p)
  )
  list(H = H, edges = edges, check_block = check_block)
}
