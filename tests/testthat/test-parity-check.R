test_that("parity-check matrix annihilates the exposure vector", {
  sim <- simulate_eccmr(n_snps = 60, block_size = 6, seed = 1)
  pc <- build_parity_check(sim$beta_X, sim$blocks)
  expect_s4_class(pc$H, "sparseMatrix")
  expect_equal(nrow(pc$H), nrow(pc$edges))
  expect_equal(length(pc$check_block), nrow(pc$H))
  ## key property: H %*% beta_X == 0 exactly
  expect_equal(as.numeric(pc$H %*% sim$beta_X),
               rep(0, nrow(pc$H)),
               tolerance = 1e-12)
})

test_that("r2 thresholding removes weakly correlated checks", {
  sim <- simulate_eccmr(n_snps = 40, block_size = 4, seed = 2)
  r2 <- matrix(0.9, 40, 40)
  diag(r2) <- 1
  r2[1, 2] <- r2[2, 1] <- 0.1  # below tau = 0.3
  pc_all <- build_parity_check(sim$beta_X, sim$blocks)
  pc_thr <- build_parity_check(sim$beta_X, sim$blocks,
                               r2_mat = r2, tau = 0.3)
  expect_lt(nrow(pc_thr$H), nrow(pc_all$H))
})

test_that("partition_ld_blocks finds connected components", {
  ld <- diag(6)
  ld[1, 2] <- ld[2, 1] <- 0.9
  ld[2, 3] <- ld[3, 2] <- 0.8   # chain: 1-2-3 one block
  ld[4, 5] <- ld[5, 4] <- 0.7
  blocks <- partition_ld_blocks(ld, r2_thresh = 0.3)
  sizes <- sort(lengths(blocks), decreasing = TRUE)
  expect_equal(sizes, c(3L, 2L, 1L))
})
