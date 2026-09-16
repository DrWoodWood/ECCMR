test_that("comparator methods are approximately unbiased without pleiotropy", {
  sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                        pi_pleio = 0, theta = 0.3, seed = 7)
  ivw <- mr_ivw(sim$beta_X, sim$beta_Y, sim$se_Y)
  med <- mr_weighted_median(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 100)
  egg <- mr_egger(sim$beta_X, sim$beta_Y, sim$se_Y)
  expect_lt(abs(ivw$theta - 0.3), 3 * ivw$se)
  expect_lt(abs(med$theta - 0.3), 0.05)
  expect_lt(abs(egg$theta - 0.3), 0.1)
})

test_that("mr_lasso flags pleiotropic instruments", {
  sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                        pi_pleio = 0.3, directional = 0.06, seed = 3)
  fit <- mr_lasso(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 50)
  expect_gt(fit$n_pleiotropic, 0)
  expect_true(is.finite(fit$theta))
})

test_that("compare_mr_methods returns one row per method and prunes LD", {
  sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                        seed = 11)
  cmp <- compare_mr_methods(sim$beta_X, sim$beta_Y, sim$se_Y,
                            blocks = sim$blocks, n_boot = 10,
                            median_boot = 50, lasso_boot = 50,
                            nlambda = 6, cv_folds = 3, seed = 1)
  expect_s3_class(cmp, "data.frame")
  expect_equal(nrow(cmp), 5L)
  expect_true(all(c("method", "theta", "se", "ci_lo", "ci_hi",
                    "n_instruments") %in% names(cmp)))
  ## ECC-MR keeps all instruments; conventional methods are LD-pruned
  expect_equal(cmp$n_instruments[cmp$method == "ECC-MR"], 120L)
  expect_equal(cmp$n_instruments[cmp$method == "IVW"], length(sim$blocks))
  expect_true(all(is.finite(cmp$theta)))
})

test_that("run_simulation aggregates and ECC-MR beats IVW under directional pleiotropy", {
  bench <- run_simulation(nsim = 5, scenario = "A", directional = 0.06,
                          n_snps = 120, block_size = 10, theta = 0.3,
                          n_boot = 10, nlambda = 6, c_gammas = c(0.5, 1, 2),
                          cv_folds = 3, seed = 99, verbose = FALSE)
  expect_s3_class(bench, "eccmr_sim")
  expect_equal(nrow(bench$summary), 5L)
  expect_equal(nrow(bench$raw), 25L)
  expect_true(all(c("bias", "rmse", "coverage") %in% names(bench$summary)))
  b_ecc <- abs(bench$summary$bias[bench$summary$method == "ECC-MR"])
  b_ivw <- abs(bench$summary$bias[bench$summary$method == "IVW"])
  expect_lt(b_ecc, b_ivw)
  expect_output(print(bench), "ECC-MR")
})
