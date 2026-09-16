## Internal: LD-block bootstrap for ECC-MR.
##
## Approximately independent LD blocks are resampled with replacement; the
## full decoding procedure is re-run on each replicate with the tuning
## parameters fixed at their selected values. Resampling whole blocks
## preserves the LD structure within blocks and propagates the uncertainty of
## error identification into the causal estimate.

## Precompute, for each block, the within-block check pairs in block-local
## indices.
.block_pairs <- function(blocks) {
  lapply(blocks, function(blk) {
    if (length(blk) < 2L) return(matrix(integer(0L), ncol = 2L))
    t(utils::combn(length(blk), 2L))
  })
}

## Build a parity-check matrix for one bootstrap replicate: `bs` is the
## vector of sampled block indices; beta vectors are concatenated block by
## block. Returns the replicate data set.
.boot_replicate <- function(beta_X, beta_Y, se_Y, blocks, pairs_list, bs) {
  n_b <- lengths(blocks)[bs]
  offset <- c(0L, cumsum(n_b))
  bx <- unlist(lapply(bs, function(k) beta_X[blocks[[k]]]), use.names = FALSE)
  by <- unlist(lapply(bs, function(k) beta_Y[blocks[[k]]]), use.names = FALSE)
  sy <- unlist(lapply(bs, function(k) se_Y[blocks[[k]]]), use.names = FALSE)

  ii <- jj <- xx <- numeric(0L)
  cb <- integer(0L)
  m_cum <- 0L
  for (q in seq_along(bs)) {
    prs <- pairs_list[[bs[q]]]
    if (nrow(prs) == 0L) next
    gi <- offset[q] + prs[, 1L]
    gj <- offset[q] + prs[, 2L]
    m_q <- nrow(prs)
    ii <- c(ii, seq_len(m_q) + m_cum, seq_len(m_q) + m_cum)
    jj <- c(jj, gi, gj)
    xx <- c(xx, bx[gj], -bx[gi])
    cb <- c(cb, rep(q, m_q))
    m_cum <- m_cum + m_q
  }
  p <- offset[length(offset)]
  H <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(m_cum, p))
  list(beta_X = bx, beta_Y = by, se_Y = sy, H = H, check_block = cb)
}

.boot_eccmr <- function(beta_X, beta_Y, se_Y, blocks, lambda, gamma,
                        B = 500L, max_iter = 200L, tol = 1e-8, cores = 1L) {
  pairs_list <- .block_pairs(blocks)
  ## per-replicate seeds generated on the master: reproducible regardless of
  ## the number of workers, and independent across replicates
  rep_seeds <- sample.int(.Machine$integer.max, B)
  one_rep <- function(sd) {
    set.seed(sd)
    bs <- sample(length(blocks), replace = TRUE)
    rep_data <- .boot_replicate(beta_X, beta_Y, se_Y, blocks, pairs_list, bs)
    fit <- ecc_fit(rep_data$beta_X, rep_data$beta_Y, rep_data$se_Y,
                   rep_data$H, lambda = lambda, gamma = gamma,
                   check_block = rep_data$check_block,
                   max_iter = max_iter, tol = tol)
    fit <- .relaxed_refit(rep_data$beta_X, rep_data$beta_Y, rep_data$se_Y,
                          rep_data$H, rep_data$check_block,
                          alpha_hat = fit$alpha, gamma = gamma,
                          max_iter = max_iter, tol = tol)
    fit$theta
  }
  unlist(.par_apply(rep_seeds, one_rep, cores), use.names = FALSE)
}
