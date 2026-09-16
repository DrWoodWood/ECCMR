#' Simulate GWAS summary statistics with LD-block structure
#'
#' Generates synthetic summary statistics mimicking the simulation scenarios
#' of the ECC-MR manuscript. Instruments are organised in equally sized LD
#' blocks.
#'
#' @param n_snps Integer; number of instruments.
#' @param block_size Integer; instruments per LD block.
#' @param theta Numeric; true causal effect.
#' @param se_y Numeric; standard error of the outcome associations.
#' @param scenario Character; one of
#'   \describe{
#'     \item{"A"}{Independent directional pleiotropy: a fraction
#'       `pi_pleio` of instruments receive direct effects with mean
#'       `directional` and SD `tau_alpha`. IVW is biased.}
#'     \item{"B"}{LD-block clustered pleiotropy with random sign: balanced
#'       across blocks, IVW approximately unbiased.}
#'     \item{"C"}{LD-block clustered, directional, individually weak
#'       pleiotropy: within a pleiotropic block, instrument i receives
#'       alpha_i = r_iu * delta with r_iu in [0.4, 1].}
#'   }
#' @param pi_pleio,pi_block,directional,tau_alpha,delta Scenario parameters.
#' @param seed Optional random seed.
#'
#' @return A list with `beta_X`, `beta_Y`, `se_Y`, `alpha` (true pleiotropic
#'   effects), `blocks` (LD-block membership) and the true `theta`.
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 100, block_size = 10, scenario = "A", seed = 1)
#' str(sim)
#'
#' @export
simulate_eccmr <- function(n_snps = 300L, block_size = 10L, theta = 0.3,
                           se_y = 0.015, scenario = c("A", "B", "C"),
                           pi_pleio = 0.3, pi_block = 0.2,
                           directional = 0.03, tau_alpha = 0.04,
                           delta = 0.02, seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  if (n_snps %% block_size != 0L) {
    n_snps <- (n_snps %/% block_size) * block_size
    warning("`n_snps` rounded down to a multiple of `block_size` (",
            n_snps, ").", call. = FALSE)
  }
  n_blocks <- n_snps %/% block_size
  blocks <- split(seq_len(n_blocks * block_size),
                  rep(seq_len(n_blocks), each = block_size))
  names(blocks) <- NULL

  beta_X <- stats::rnorm(n_snps, mean = 0.08, sd = 0.02)
  beta_X[abs(beta_X) < 0.02] <- 0.03

  alpha <- numeric(n_snps)
  if (scenario == "A") {
    idx <- sample(n_snps, size = floor(pi_pleio * n_snps))
    alpha[idx] <- stats::rnorm(length(idx), mean = directional, sd = tau_alpha)
  } else if (scenario == "B") {
    for (blk in blocks) {
      if (stats::runif(1L) < pi_block) {
        r_iu <- stats::runif(length(blk), 0.5, 1) *
          sample(c(-1, 1), length(blk), replace = TRUE)
        alpha[blk] <- r_iu * delta
      }
    }
  } else {
    for (blk in blocks) {
      if (stats::runif(1L) < pi_block) {
        r_iu <- stats::runif(length(blk), 0.4, 1)
        alpha[blk] <- r_iu * delta
      }
    }
  }

  se_Y <- rep(se_y, n_snps)
  beta_Y <- theta * beta_X + alpha + stats::rnorm(n_snps, sd = se_y)
  list(beta_X = beta_X, beta_Y = beta_Y, se_Y = se_Y, alpha = alpha,
       blocks = blocks, theta = theta)
}
