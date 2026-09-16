## Internal: tuning-parameter selection by K-fold cross-validation over LD
## blocks. In each fold the model is fitted without the held-out blocks and
## the weighted prediction error of beta_Y on the held-out blocks is
## recorded; instruments in held-out blocks contribute their pleiotropy-free
## prediction theta * beta_X.

## Data-driven lambda grid: from lambda_max down geometrically. lambda_max
## is the largest weighted IVW-residual score, a conservative upper bound on
## the smallest lambda that sets every alpha to zero. The grid is derived
## from the data score only; the syndrome term must not inflate it.
.lambda_grid <- function(beta_X, beta_Y, se_Y, H, check_block,
                         nlambda = 20L, lambda_min_ratio = 0.01) {
  w <- 1 / se_Y^2
  theta0 <- .ivw(beta_X, beta_Y, se_Y)
  b0 <- w * (beta_Y - theta0 * beta_X)
  lambda_max <- max(abs(b0))
  lambda_max * lambda_min_ratio^seq(0, 1, length.out = nlambda)
}

## Natural scale for gamma: tr(W) / tr(H' Omega^+ H). Returns 0 when there
## are no checks (ECC-MR then reduces to a sparse-pleiotropy model).
.gamma0 <- function(se_Y, H, check_block) {
  m <- nrow(H)
  if (m == 0L) return(0)
  w <- 1 / se_Y^2
  Sigma <- Matrix::Diagonal(x = se_Y^2)
  Omega <- H %*% Matrix::tcrossprod(Sigma, H)
  Minv <- .pinv_blockdiag(Omega, check_block)
  MH <- Minv %*% H
  tr_HtMH <- sum(Matrix::colSums(H * MH))
  sum(w) / max(tr_HtMH, .Machine$double.eps)
}

## Restrict the problem to a subset of instruments, keeping only checks whose
## two endpoints both lie in the subset; column indices are remapped to the
## subset positions and H is rebuilt from the edge list.
.subset_problem <- function(beta_X, beta_Y, se_Y, edges, keep_snps) {
  if (nrow(edges) > 0L) {
    keep_row <- edges[, 1L] %in% keep_snps & edges[, 2L] %in% keep_snps
    edges_sub <- edges[keep_row, , drop = FALSE]
    remap <- integer(length(beta_X))
    remap[keep_snps] <- seq_along(keep_snps)
    ## rebuild H directly from the kept edges
    m_sub <- nrow(edges_sub)
    gi <- remap[edges_sub[, 1L]]
    gj <- remap[edges_sub[, 2L]]
    bx <- beta_X[keep_snps]
    H_sub <- Matrix::sparseMatrix(
      i = c(seq_len(m_sub), seq_len(m_sub)),
      j = c(gi, gj),
      x = c(bx[gj], -bx[gi]),
      dims = c(m_sub, length(keep_snps))
    )
    cb_sub <- .check_blocks(H_sub)
  } else {
    H_sub <- Matrix::sparseMatrix(i = integer(0L), j = integer(0L),
                                  x = numeric(0L),
                                  dims = c(0L, length(keep_snps)))
    cb_sub <- integer(0L)
  }
  list(beta_X = beta_X[keep_snps], beta_Y = beta_Y[keep_snps],
       se_Y = se_Y[keep_snps], H = H_sub, check_block = cb_sub)
}


## Internal: parallel lapply over a PSOCK cluster (works on Windows).
## Package internals resolve via the loaded namespace on the workers.
.par_apply <- function(X, FUN, cores) {
  cores <- max(1L, min(as.integer(cores), parallel::detectCores(TRUE)))
  if (cores <= 1L || length(X) < 8L) return(lapply(X, FUN))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(eccmr))
  ## dynamic load balancing: workers fetch the next chunk as soon as they
  ## finish, so slow replicates never leave the cluster half-idle
  parallel::parLapplyLB(cl, X, FUN,
                        chunk.size = max(1L, length(X) %/% (cores * 4L)))
}

## K-fold CV over LD blocks. Returns a data.frame with columns c_gamma,
## lambda, cv_error, cv_se.
.cv_eccmr <- function(beta_X, beta_Y, se_Y, H, edges, check_block, blocks,
                      lambdas, c_gammas, cv_folds = 10L,
                      max_iter = 200L, tol = 1e-8, cores = 1L) {
  w <- 1 / se_Y^2
  n_blocks <- length(blocks)
  cv_folds <- min(cv_folds, n_blocks)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = n_blocks))
  all_snps <- unlist(blocks, use.names = FALSE)

  grid <- expand.grid(c_gamma = c_gammas, lambda = lambdas)
  err <- matrix(NA_real_, nrow = nrow(grid), ncol = cv_folds)

  ## precompute per-fold sub-problems on the master (fold split uses the
  ## master RNG, so results are identical for any number of cores)
  subs <- lapply(seq_len(cv_folds), function(f) {
    train_snps <- setdiff(all_snps,
                          unlist(blocks[fold_id == f], use.names = FALSE))
    sub <- .subset_problem(beta_X, beta_Y, se_Y, edges, train_snps)
    sub$g0 <- .gamma0(sub$se_Y, sub$H, sub$check_block)
    sub
  })
  tasks <- expand.grid(f = seq_len(cv_folds), g = seq_len(nrow(grid)))
  vals <- .par_apply(seq_len(nrow(tasks)), function(i) {
    f <- tasks$f[i]; g <- tasks$g[i]
    sub <- subs[[f]]
    fit <- ecc_fit(sub$beta_X, sub$beta_Y, sub$se_Y, sub$H,
                   lambda = grid$lambda[g], gamma = grid$c_gamma[g] * sub$g0,
                   check_block = sub$check_block,
                   max_iter = max_iter, tol = tol)
    test_snps <- unlist(blocks[fold_id == f], use.names = FALSE)
    pred <- fit$theta * beta_X[test_snps]
    sum(w[test_snps] * (beta_Y[test_snps] - pred)^2)
  }, cores)
  err[cbind(tasks$g, tasks$f)] <- unlist(vals, use.names = FALSE)
  grid$cv_error <- rowSums(err)
  grid$cv_se <- apply(err, 1L, stats::sd) / sqrt(cv_folds)
  grid
}

## Tuning by BIC along the (c_gamma, lambda) path. For each candidate the
## model is refitted on the full data and scored by
##   BIC = weighted SSE + gamma * syndrome quadratic + log(p) * ||alpha||_0 .
## Unlike held-out prediction (CV over LD blocks), the in-sample fit rewards
## modelling pleiotropy, which makes the path strongly informative even when
## pleiotropic effects are unpredictable across blocks. Returns a data.frame
## with columns c_gamma, lambda, bic, df, theta.
.bic_eccmr <- function(beta_X, beta_Y, se_Y, H, check_block,
                       lambdas, c_gammas, max_iter = 200L, tol = 1e-8,
                       cores = 1L) {
  w <- 1 / se_Y^2
  p <- length(beta_X)
  m <- nrow(H)
  need_syn <- m > 0L && any(c_gammas != 0)
  if (need_syn) {
    Sigma <- Matrix::Diagonal(x = se_Y^2)
    s <- .spmv(H, beta_Y)
    Omega <- H %*% Matrix::tcrossprod(Sigma, H)
    Minv <- .pinv_blockdiag(Omega, check_block)
  } else {
    s <- numeric(0L)
    Minv <- NULL
  }
  g0 <- if (need_syn) .gamma0(se_Y, H, check_block) else 0

  grid <- expand.grid(c_gamma = c_gammas, lambda = lambdas)
  cells <- .par_apply(seq_len(nrow(grid)), function(g) {
    fit <- ecc_fit(beta_X, beta_Y, se_Y, H,
                   lambda = grid$lambda[g],
                   gamma = grid$c_gamma[g] * g0,
                   check_block = check_block,
                   max_iter = max_iter, tol = tol)
    res <- beta_Y - fit$theta * beta_X - fit$alpha
    sse <- sum(w * res^2)
    if (!is.null(Minv) && grid$c_gamma[g] != 0) {
      rs <- s - .spmv(H, fit$alpha)
      sse <- sse + grid$c_gamma[g] * g0 * sum(rs * .spmv(Minv, rs))
    }
    c(bic = sse + log(p) * sum(abs(fit$alpha) > 0),
      df = sum(abs(fit$alpha) > 0), theta = fit$theta)
  }, cores)
  grid$bic   <- vapply(cells, "[[", numeric(1), 1L)
  grid$df    <- vapply(cells, "[[", numeric(1), 2L)
  grid$theta <- vapply(cells, "[[", numeric(1), 3L)
  list(table = grid, g0 = g0)
}

## Syndrome-guided gamma selection.
##
## Motivation: BIC penalises the number of non-zero alpha, so it always
## prefers the weakest syndrome-consistency penalty even when correlated
## noise (LD / sample overlap) calls for strong correction. The ECC-native
## rule instead increases gamma until the *residual* syndrome
## (s - H alpha_hat) is consistent with pure noise, i.e. the global syndrome
## test is no longer rejected: "raise the correction strength until no
## decodable error remains".
##
## lambda is fixed to the BIC-selected value; only gamma is scanned.
## Returns the smallest c_gamma whose residual syndrome statistic is below
## the (1 - alpha_level) chi-squared quantile; falls back to the largest
## c_gamma when none qualifies, and to the smallest when even the weakest
## penalty already passes (iid-noise case -> gamma stays small).
.select_gamma_syndrome <- function(beta_X, beta_Y, se_Y, H, check_block,
                                   lambda, c_gammas = c(0.1, 0.2, 0.5, 1, 2,
                                                        5, 10, 50, 100),
                                   alpha_level = 0.05,
                                   max_iter = 200L, tol = 1e-8) {
  g0 <- .gamma0(se_Y, H, check_block)
  m <- nrow(H)
  if (m == 0L) return(list(c_gamma = c_gammas[1], gamma = c_gammas[1] * g0,
                           table = NULL))
  Sigma <- Matrix::Diagonal(x = se_Y^2)
  s <- .spmv(H, beta_Y)
  Omega <- H %*% Matrix::tcrossprod(Sigma, H)
  Minv <- .pinv_blockdiag(Omega, check_block)
  df <- sum(eigen(as.matrix(Omega), symmetric = TRUE,
                  only.values = TRUE)$values > 1e-10 * max(diag(as.matrix(Omega))))
  crit <- stats::qchisq(1 - alpha_level, df = df)

  rows <- lapply(c_gammas, function(cg) {
    fit <- ecc_fit(beta_X, beta_Y, se_Y, H, lambda = lambda,
                   gamma = cg * g0, check_block = check_block,
                   max_iter = max_iter, tol = tol)
    rs <- s - .spmv(H, fit$alpha)
    t_res <- sum(rs * .spmv(Minv, rs))
    data.frame(c_gamma = cg, theta = fit$theta,
               n_pleio = sum(abs(fit$alpha) > 0),
               syndrome_stat = t_res, pass = t_res <= crit)
  })
  tab <- do.call(rbind, rows)
  ok <- which(tab$pass)
  cg_sel <- if (length(ok)) tab$c_gamma[ok[1]] else max(c_gammas)
  list(c_gamma = cg_sel, gamma = cg_sel * g0, table = tab,
       df = df, crit = crit)
}

## Relaxed (debiased) refit: keep the support selected by the L1 fit and
## re-estimate theta and the supported alpha without the L1 penalty (the
## syndrome-consistency penalty is retained). Off-support coordinates are
## forced to zero via an infinite per-coordinate penalty. This removes the
## L1 shrinkage bias that would otherwise leak into the causal estimate.
.relaxed_refit <- function(beta_X, beta_Y, se_Y, H, check_block,
                           alpha_hat, gamma, max_iter = 200L, tol = 1e-8) {
  ## large finite penalty (not Inf, which would make the objective NaN)
  lambda_vec <- ifelse(abs(alpha_hat) > 0, 0, 1e300)
  ecc_fit(beta_X, beta_Y, se_Y, H, lambda = lambda_vec, gamma = gamma,
          check_block = check_block, alpha_init = alpha_hat,
          max_iter = max_iter, tol = tol)
}
