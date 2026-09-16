#' Fit the ECC-MR model for fixed tuning parameters
#'
#' Minimizes the convex ECC-MR objective
#'
#'   0.5 * sum_i w_i (beta_Yi - theta beta_Xi - alpha_i)^2
#'     + lambda * ||alpha||_1
#'     + (gamma / 2) * ||s - H alpha||^2_{Omega^+}
#'
#' by block coordinate descent: a closed-form inverse-variance weighted update
#' for the causal effect theta alternates with cyclic soft-thresholding
#' coordinate-descent sweeps for the pleiotropic effects alpha. Joint
#' convexity guarantees convergence to the global minimum. One coordinate
#' sweep costs O(number of checks), so decoding remains feasible for large
#' instrument sets.
#'
#' @param beta_X,beta_Y Numeric vectors of length p; exposure and outcome
#'   association estimates.
#' @param se_Y Numeric vector of length p; standard errors of `beta_Y`.
#' @param H Sparse m x p parity-check matrix, from [build_parity_check()].
#' @param lambda Numeric; sparsity (L1) tuning parameter. May be a vector of
#'   length p giving a per-instrument penalty (use a very large value such
#'   as `1e300` to force `alpha_i = 0`); a scalar is recycled.
#' @param gamma Numeric; syndrome-consistency tuning parameter.
#' @param check_block Optional integer vector assigning each check to an LD
#'   block; derived from `H` when omitted.
#' @param sigma Optional p x p error covariance matrix (use when the exposure
#'   and outcome GWAS overlap); defaults to `diag(se_Y^2)`.
#' @param theta_init,alpha_init Optional starting values; default to the IVW
#'   estimate and zero.
#' @param max_iter Integer; maximum number of outer iterations.
#' @param tol Numeric; convergence tolerance on the relative objective change.
#' @param tol_eig Numeric; relative eigenvalue tolerance for the pseudoinverse.
#'
#' @return A list with components
#'   \describe{
#'     \item{theta}{Estimated causal effect.}
#'     \item{alpha}{Estimated pleiotropic effects (length p).}
#'     \item{iter}{Number of outer iterations used.}
#'     \item{converged}{Logical; whether the tolerance was reached.}
#'     \item{objective}{Final objective value.}
#'   }
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' pc <- build_parity_check(sim$beta_X, sim$blocks)
#' fit <- ecc_fit(sim$beta_X, sim$beta_Y, sim$se_Y, pc$H,
#'                lambda = 5, gamma = 1)
#' fit$theta
#'
#' @export
ecc_fit <- function(beta_X, beta_Y, se_Y, H, lambda, gamma,
                    check_block = NULL, sigma = NULL,
                    theta_init = NULL, alpha_init = NULL,
                    max_iter = 200L, tol = 1e-8, tol_eig = 1e-10) {
  p <- length(beta_X)
  m <- nrow(H)
  w <- 1 / se_Y^2
  if (length(lambda) == 1L) lambda <- rep(lambda, p)
  stopifnot(length(lambda) == p)

  ## v0.3.0: when gamma = 0 the syndrome penalty vanishes and NONE of the
  ## H-based matrices (Omega, Minv, HtMH -- O(m^2) memory) are needed; the
  ## objective reduces to a diagonal quadratic. This matters at biobank
  ## scale: m can exceed 7e4, where densifying Omega would need >30 GiB.
  if (m > 0L && gamma > 0) {
    if (is.null(check_block)) check_block <- .check_blocks(H)
    Sigma <- if (is.null(sigma)) Matrix::Diagonal(x = se_Y^2) else sigma
    s <- .spmv(H, beta_Y)
    Omega <- H %*% Matrix::tcrossprod(Sigma, H)
    Minv <- .pinv_blockdiag(Omega, check_block, tol_eig = tol_eig)
    MH <- Minv %*% H                                # m x p
    HtMH <- Matrix::crossprod(H, MH)                # p x p, block diagonal
    ## NOTE: use Matrix::crossprod / Matrix::colSums rather than base t() /
    ## colSums() on sparse matrices -- the base S3 generics do not reliably
    ## dispatch to Matrix's S4 methods from inside an installed package.
    HtMs <- as.numeric(Matrix::crossprod(H, Minv %*% s))
    Gd <- w + gamma * Matrix::colSums(H * MH)       # diag(W + gamma H' M H)
    ## crossprod() returns a symmetric-storage dsCMatrix (lower triangle
    ## only); expand to a general dgCMatrix so that the column triplets
    ## below contain every entry of each column.
    G <- methods::as(Matrix::Diagonal(x = w) + gamma * HtMH, "dgCMatrix")
    b0 <- w * beta_Y + gamma * HtMs
  } else {
    s <- numeric(0L)
    Minv <- NULL
    Gd <- w
    G <- Matrix::Diagonal(x = w)
    b0 <- w * beta_Y
  }

  ## Precompute the column structure of G for fast coordinate updates.
  ## The diagonal of G is strictly positive (w_i > 0), so every column is
  ## present; factor levels keep the columns aligned with 1..p.
  Gt <- methods::as(G, "TsparseMatrix")
  col_ix <- split(seq_along(Gt@x), factor(Gt@j + 1L, levels = seq_len(p)))
  col_rows <- lapply(col_ix, function(ix) Gt@i[ix] + 1L)
  col_vals <- lapply(col_ix, function(ix) Gt@x[ix])

  theta <- if (is.null(theta_init)) .ivw(beta_X, beta_Y, se_Y) else theta_init
  alpha <- if (is.null(alpha_init)) rep(0, p) else alpha_init

  objective <- function(theta, alpha) {
    res <- beta_Y - theta * beta_X - alpha
    quad <- 0
    if (!is.null(Minv)) {          # gamma > 0 only (v0.3.1)
      rs <- s - .spmv(H, alpha)
      quad <- sum(rs * .spmv(Minv, rs))
    }
    0.5 * sum(w * res^2) + sum(lambda * abs(alpha)) + 0.5 * gamma * quad
  }

  obj_prev <- Inf
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    ## theta update: closed-form IVW on corrected outcome associations
    theta <- sum(w * beta_X * (beta_Y - alpha)) / sum(w * beta_X^2)
    ## alpha update: cyclic coordinate descent with soft thresholding,
    ## maintaining Ga = G %*% alpha incrementally
    b <- b0 - w * (theta * beta_X)
    Ga <- .spmv(G, alpha)
    for (i in seq_len(p)) {
      ri <- b[i] - Ga[i] + Gd[i] * alpha[i]
      a_new <- .soft_threshold(ri, lambda[i]) / Gd[i]
      dlt <- a_new - alpha[i]
      if (dlt != 0) {
        Ga[col_rows[[i]]] <- Ga[col_rows[[i]]] + col_vals[[i]] * dlt
        alpha[i] <- a_new
      }
    }
    obj <- objective(theta, alpha)
    if (abs(obj_prev - obj) < tol * (1 + abs(obj_prev))) {
      converged <- TRUE
      obj_prev <- obj
      break
    }
    obj_prev <- obj
  }
  list(theta = theta, alpha = alpha, iter = iter,
       converged = converged, objective = obj_prev)
}
