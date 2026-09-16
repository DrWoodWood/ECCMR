test_that("block pseudoinverse satisfies M Omega M = M", {
  sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 5)
  pc <- build_parity_check(sim$beta_X, sim$blocks)
  Sigma <- Matrix::Diagonal(x = sim$se_Y^2)
  Omega <- pc$H %*% Matrix::tcrossprod(Sigma, pc$H)
  Minv <- .pinv_blockdiag(Omega, pc$check_block)
  resid <- as.matrix(Minv %*% Omega %*% Minv - Minv)
  expect_lt(max(abs(resid)), 1e-6 * max(abs(as.matrix(Minv))))
})

test_that("gamma natural scale is of order 1 (whitened parametrization)", {
  sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A", seed = 6)
  pc <- build_parity_check(sim$beta_X, sim$blocks)
  g0 <- .gamma0(sim$se_Y, pc$H, pc$check_block)
  ## tr(H' Omega^+ H) = rank(H) / se_y^2 for homoscedastic errors, so
  ## g0 = sum(w) / tr(H'Omega^+H) = p / rank(H) ~ 1
  expect_gt(g0, 0.5)
  expect_lt(g0, 5)
})

test_that("BIC tuning detects pleiotropic instruments and debiased refit beats IVW", {
  sim <- simulate_eccmr(n_snps = 300, block_size = 10, theta = 0.3,
                        scenario = "A", directional = 0.05, seed = 42)
  fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
               n_boot = 20, seed = 1)
  expect_gt(length(fit$pleiotropic_index), 30)
  expect_lt(abs(fit$theta - 0.3), abs(fit$theta_l1 - 0.3) + 1e-8)
  ivw <- mr_ivw(sim$beta_X, sim$beta_Y, sim$se_Y)
  expect_lt(abs(fit$theta - 0.3), abs(ivw$theta - 0.3))
})
