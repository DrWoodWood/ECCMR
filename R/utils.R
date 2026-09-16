## Internal utilities for the eccmr package.

## Multiply a (sparse) matrix by a vector and return a plain numeric vector.
.spmv <- function(A, x) {
  as(A %*% x, "matrix")[, 1L]
}

## Soft-thresholding operator S(z, kappa) = sign(z) * max(|z| - kappa, 0).
.soft_threshold <- function(z, kappa) {
  sign(z) * max(abs(z) - kappa, 0)
}

## Connected components of a graph given by an edge list (two-column matrix
## of vertex indices in 1..n). Returns an integer vector of component labels.
## NOTE: only plain local assignments are used (no `<<-`): inside the
## function body `parent[...] <<- ...` would modify the *global* environment
## instead of the local `parent`, and it can fail outright once the package
## is byte-compiled at installation.
.union_find_components <- function(edges, n) {
  parent <- seq_len(n)
  find_root <- function(x) {
    while (parent[x] != x) x <- parent[x]
    x
  }
  if (length(edges) > 0L) {
    for (k in seq_len(nrow(edges))) {
      ra <- find_root(edges[k, 1L])
      rb <- find_root(edges[k, 2L])
      if (ra != rb) parent[ra] <- rb
    }
  }
  roots <- vapply(seq_len(n), find_root, integer(1L))
  match(roots, unique(roots))
}

## Derive a block label for every check (row) of a parity-check matrix:
## two checks belong to the same group if they share an instrument. Because
## H is built from disjoint LD blocks, this recovers the block structure.
.check_blocks <- function(H) {
  m <- nrow(H)
  if (m == 0L) return(integer(0L))
  Ht <- methods::as(H, "TsparseMatrix")  # triplets: 0-based @i (row), @j (col)
  cols <- split(Ht@j + 1L, Ht@i + 1L)    # two instrument indices per check
  pairs <- do.call(rbind, lapply(cols, sort))
  storage.mode(pairs) <- "integer"
  .union_find_components(pairs, m)
}

## Moore-Penrose pseudoinverse of a symmetric block-diagonal matrix, computed
## block by block through eigendecomposition. `check_block` gives the block
## label of each row. Returns a sparse block-diagonal Matrix whose rows and
## columns keep the original ordering.
.pinv_blockdiag <- function(Omega, check_block, tol_eig = 1e-10) {
  m <- nrow(Omega)
  idx_split <- split(seq_len(m), check_block)
  trips <- lapply(idx_split, function(idx) {
    Ob <- as.matrix(Omega[idx, idx, drop = FALSE])
    ee <- eigen(Ob, symmetric = TRUE)
    d <- ee$values
    keep <- d > tol_eig * max(d, 0)
    if (!any(keep)) return(NULL)
    ## Moore-Penrose pseudoinverse: U diag(1/d) U'. Scale the eigenvectors
    ## by 1/sqrt(d) BEFORE the outer product -- scaling by 1/d instead would
    ## erroneously yield U diag(1/d^2) U'.
    U <- ee$vectors[, keep, drop = FALSE]
    U <- sweep(U, 2L, 1 / sqrt(d[keep]), `*`)
    Mb <- tcrossprod(U)
    gg <- expand.grid(seq_along(idx), seq_along(idx))
    data.frame(i = idx[gg[, 1L]], j = idx[gg[, 2L]], x = as.numeric(Mb))
  })
  trips <- do.call(rbind, trips)
  Matrix::sparseMatrix(i = trips$i, j = trips$j, x = trips$x,
                       dims = c(m, m))
}

## Inverse-variance weighted (IVW) estimator, used for initialisation.
.ivw <- function(beta_X, beta_Y, se_Y) {
  w <- 1 / se_Y^2
  sum(w * beta_X * beta_Y) / sum(w * beta_X^2)
}
