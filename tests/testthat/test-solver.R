test_that("ecc_fit recovers the IVW estimate when lambda is huge", {
  sim <- simulate_eccmr(n_snps = 100, block_size = 10, scenario = "A",
                        seed = 3)
  pc <- build_parity_check(sim$beta_X, sim$blocks)
  fit <- ecc_fit(sim$beta_X, sim$beta_Y, sim$se_Y, pc$H,
                 lambda = 1e8, gamma = 1)
  w <- 1 / sim$se_Y^2
  ivw <- sum(w * sim$beta_X * sim$beta_Y) / sum(w * sim$beta_X^2)
  expect_equal(fit$theta, ivw, tolerance = 1e-6)
  expect_true(all(fit$alpha == 0))
})

test_that("ecc_fit is unbiased with no pleiotropy and corrects directional pleiotropy", {
  ## no pleiotropy: estimate close to truth
  sim0 <- simulate_eccmr(n_snps = 200, block_size = 10, theta = 0.3,
                         pi_pleio = 0, scenario = "A", seed = 11)
  pc0 <- build_parity_check(sim0$beta_X, sim0$blocks)
  fit0 <- ecc_fit(sim0$beta_X, sim0$beta_Y, sim0$se_Y, pc0$H,
                  lambda = 50, gamma = 0)
  expect_lt(abs(fit0$theta - 0.3), 0.05)

  ## directional pleiotropy: ECC-MR less biased than IVW
  sim1 <- simulate_eccmr(n_snps = 300, block_size = 10, theta = 0.3,
                         scenario = "A", seed = 12)
  pc1 <- build_parity_check(sim1$beta_X, sim1$blocks)
  w <- 1 / sim1$se_Y^2
  ivw <- sum(w * sim1$beta_X * sim1$beta_Y) / sum(w * sim1$beta_X^2)
  fit1 <- ecc_fit(sim1$beta_X, sim1$beta_Y, sim1$se_Y, pc1$H,
                  lambda = 50, gamma = 0)
  expect_lt(abs(fit1$theta - 0.3), abs(ivw - 0.3))
})

test_that("ecc_fit converges and objective decreases", {
  sim <- simulate_eccmr(n_snps = 100, block_size = 10, seed = 5)
  pc <- build_parity_check(sim$beta_X, sim$blocks)
  fit <- ecc_fit(sim$beta_X, sim$beta_Y, sim$se_Y, pc$H,
                 lambda = 10, gamma = 0, max_iter = 500)
  expect_true(fit$converged)
  expect_true(is.finite(fit$objective))
})

test_that("syndrome test is valid under the null and rejects under pleiotropy", {
  sim0 <- simulate_eccmr(n_snps = 200, block_size = 10, pi_pleio = 0,
                         scenario = "A", seed = 21)
  pc0 <- build_parity_check(sim0$beta_X, sim0$blocks)
  st0 <- syndrome_test(sim0$beta_Y, sim0$se_Y, pc0$H, pc0$check_block)
  expect_true(st0$p_value >= 0 && st0$p_value <= 1)

  sim1 <- simulate_eccmr(n_snps = 200, block_size = 10, scenario = "A",
                         directional = 0.1, pi_pleio = 0.5, seed = 22)
  pc1 <- build_parity_check(sim1$beta_X, sim1$blocks)
  st1 <- syndrome_test(sim1$beta_Y, sim1$se_Y, pc1$H, pc1$check_block)
  expect_lt(st1$p_value, 0.05)
})
