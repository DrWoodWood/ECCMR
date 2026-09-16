test_that("eccmr runs end-to-end and returns a valid object", {
  sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                        seed = 42)
  fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
               nlambda = 10, c_gammas = c(0.5, 1, 2),
               n_boot = 30, seed = 1)
  expect_s3_class(fit, "eccmr")
  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$se) && fit$se > 0)
  expect_length(fit$alpha, 120)
  expect_length(fit$boot_theta, 30)
  expect_output(print(fit), "ECC-MR")
  expect_output(summary(fit), "Causal estimate")
})

test_that("eccmr is less biased than IVW under directional pleiotropy", {
  sim <- simulate_eccmr(n_snps = 300, block_size = 10, theta = 0.3,
                        scenario = "A", directional = 0.05, seed = 100)
  fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
               nlambda = 10, c_gammas = c(0.5, 1, 2), n_boot = 0, seed = 1)
  w <- 1 / sim$se_Y^2
  ivw <- sum(w * sim$beta_X * sim$beta_Y) / sum(w * sim$beta_X^2)
  expect_lt(abs(fit$theta - 0.3), abs(ivw - 0.3))
})

test_that("eccmr works with ld_mat input and without checks", {
  sim <- simulate_eccmr(n_snps = 60, block_size = 6, seed = 7)
  ld <- matrix(0.9, 60, 60); diag(ld) <- 1
  fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, ld_mat = ld,
               n_boot = 0, seed = 2)
  expect_s3_class(fit, "eccmr")

  ## no LD at all -> no checks -> reduces to a sparse-pleiotropy model
  ld0 <- diag(60)
  fit0 <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, ld_mat = ld0,
                n_boot = 0, seed = 2)
  expect_equal(fit0$n_checks, 0L)
  expect_true(is.finite(fit0$theta))
})
