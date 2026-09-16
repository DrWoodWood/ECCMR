test_that("run_boundary_scan returns one row per regime x method", {
  bd <- run_boundary_scan(pis = 0.1, directionals = c(0.05, 0),
                          tau_alphas = 0.08, nsim = 2L,
                          methods = c("ECC-MR", "IVW"),
                          nlambda = 8L, c_gammas = c(1, 2),
                          verbose = FALSE)
  expect_s3_class(bd, "eccmr_boundary")
  expect_equal(nrow(bd), 4L)
  expect_true(all(c("directional", "pi", "tau_alpha", "method",
                    "bias", "mc_sd", "rmse") %in% names(bd)))
})

test_that("ECC-MR dominates in the detectable-pleiotropy regime", {
  bd <- run_boundary_scan(pis = 0.1, directionals = 0.05,
                          tau_alphas = 0.08, nsim = 4L,
                          methods = c("ECC-MR", "IVW", "MR-Lasso"),
                          nlambda = 8L, c_gammas = c(1, 2),
                          seed = 7L, verbose = FALSE)
  ecc <- bd$rmse[bd$method == "ECC-MR"]
  others <- bd$rmse[bd$method != "ECC-MR"]
  expect_lt(ecc, min(others))
})
