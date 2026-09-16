## Reference implementations of mainstream summary-data MR methods, used as
## comparators in simulation benchmarks. All methods use only (beta_X, beta_Y,
## se_Y) and base R, so the package stays dependency-free. For production
## analyses we recommend the dedicated packages (TwoSampleMR,
## MendelianRandomization, MRPRESSO); the versions here follow the published
## estimators closely and are intended for benchmarking ECC-MR.

## Standard normal confidence interval.
.ci_normal <- function(theta, se) {
  theta + c(-1, 1) * stats::qnorm(0.975) * se
}

#' Inverse-variance weighted (IVW) Mendelian randomization
#'
#' Multiplicative random-effects IVW estimator: the standard error is inflated
#' by the square root of the residual heterogeneity statistic when
#' Cochran's Q exceeds its degrees of freedom. Instruments are assumed
#' independent.
#'
#' @param beta_X,beta_Y Numeric vectors of exposure and outcome association
#'   estimates.
#' @param se_Y Numeric vector of standard errors of `beta_Y`.
#'
#' @return A list with `theta`, `se`, `p_value`, confidence interval
#'   (`ci_lo`, `ci_hi`), Cochran's `Q` and the residual inflation factor
#'   `phi`.
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' mr_ivw(sim$beta_X, sim$beta_Y, sim$se_Y)$theta
#'
#' @export
mr_ivw <- function(beta_X, beta_Y, se_Y) {
  w <- 1 / se_Y^2
  p <- length(beta_X)
  theta <- sum(w * beta_X * beta_Y) / sum(w * beta_X^2)
  q <- sum(w * (beta_Y - theta * beta_X)^2)
  phi <- max(1, q / (p - 1))
  se <- sqrt(phi / sum(w * beta_X^2))
  ci <- .ci_normal(theta, se)
  list(theta = theta, se = se, p_value = 2 * stats::pnorm(-abs(theta / se)),
       ci_lo = ci[1L], ci_hi = ci[2L], Q = q, phi = phi)
}

#' MR-Egger regression
#'
#' Weighted regression of `beta_Y` on `beta_X` with an unconstrained
#' intercept. The slope is a consistent causal estimate under the InSIDE
#' assumption; the intercept tests for directional pleiotropy. Instruments
#' are assumed independent.
#'
#' @inheritParams mr_ivw
#'
#' @return A list with the causal estimate `theta` (slope), its `se` and
#'   `p_value`, confidence interval (`ci_lo`, `ci_hi`), and the pleiotropy
#'   intercept with its standard error and p-value.
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' mr_egger(sim$beta_X, sim$beta_Y, sim$se_Y)$theta
#'
#' @export
mr_egger <- function(beta_X, beta_Y, se_Y) {
  w <- 1 / se_Y^2
  p <- length(beta_X)
  xbar <- sum(w * beta_X) / sum(w)
  ybar <- sum(w * beta_Y) / sum(w)
  sxx <- sum(w * (beta_X - xbar)^2)
  slope <- sum(w * (beta_X - xbar) * (beta_Y - ybar)) / sxx
  intercept <- ybar - slope * xbar
  rss <- sum(w * (beta_Y - intercept - slope * beta_X)^2)
  phi <- max(1, rss / (p - 2))
  se_slope <- sqrt(phi / sxx)
  se_int <- sqrt(phi * (1 / sum(w) + xbar^2 / sxx))
  ci <- .ci_normal(slope, se_slope)
  list(theta = slope, se = se_slope,
       p_value = 2 * stats::pnorm(-abs(slope / se_slope)),
       ci_lo = ci[1L], ci_hi = ci[2L],
       intercept = intercept, se_intercept = se_int,
       p_intercept = 2 * stats::pnorm(-abs(intercept / se_int)))
}

## Weighted median of x with weights w (linear interpolation).
.wtd_median <- function(x, w) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  stats::approx(cw, x[o], xout = 0.5, ties = "ordered")$y
}

#' Weighted median Mendelian randomization
#'
#' Weighted median of the SNP-specific ratio estimates with weights
#' `beta_X^2 / se_Y^2`. Consistent when up to half of the weight comes from
#' invalid instruments. The standard error is obtained by bootstrapping
#' instruments.
#'
#' @inheritParams mr_ivw
#' @param n_boot Integer; number of bootstrap replicates for the standard
#'   error (default 200).
#'
#' @return A list with `theta`, `se`, `p_value` and confidence interval
#'   (`ci_lo`, `ci_hi`).
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' mr_weighted_median(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 50)$theta
#'
#' @export
mr_weighted_median <- function(beta_X, beta_Y, se_Y, n_boot = 200L) {
  ratios <- beta_Y / beta_X
  w <- beta_X^2 / se_Y^2
  theta <- .wtd_median(ratios, w)
  se <- NA_real_
  if (n_boot > 0L) {
    p <- length(beta_X)
    boots <- replicate(n_boot, {
      idx <- sample(p, replace = TRUE)
      .wtd_median(ratios[idx], w[idx])
    })
    se <- stats::sd(boots)
  }
  ci <- .ci_normal(theta, se)
  list(theta = theta, se = se,
       p_value = 2 * stats::pnorm(-abs(theta / se)),
       ci_lo = ci[1L], ci_hi = ci[2L])
}

#' MR-Lasso
#'
#' Lasso-penalized IVW: an L1 penalty on the pleiotropic effects identifies
#' and removes invalid instruments. The penalty parameter is chosen by the
#' heterogeneity stopping rule of the original publication: decrease lambda
#' from its maximal value until Cochran's Q computed on the corrected
#' associations is no longer significant at the 5% level. The standard error
#' is obtained by bootstrapping instruments.
#'
#' @inheritParams mr_ivw
#' @param nlambda Integer; size of the lambda grid (default 25).
#' @param max_iter Integer; maximum coordinate-descent iterations per lambda.
#' @param n_boot Integer; number of bootstrap replicates for the standard
#'   error (default 200).
#'
#' @return A list with `theta`, `se`, `p_value`, confidence interval
#'   (`ci_lo`, `ci_hi`), the selected `lambda` and the number of instruments
#'   flagged as pleiotropic (`n_pleiotropic`).
#'
#' @examples
#' sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
#' mr_lasso(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 50)$theta
#'
#' @export
mr_lasso <- function(beta_X, beta_Y, se_Y, nlambda = 25L, max_iter = 100L,
                     n_boot = 200L) {
  fit_one <- function(bx, by, sy) {
    w <- 1 / sy^2
    p <- length(bx)
    theta0 <- sum(w * bx * by) / sum(w * bx^2)
    lam_max <- max(abs(w * (by - theta0 * bx)))
    lambdas <- lam_max * 0.01^seq(0, 1, length.out = nlambda)
    alpha <- numeric(p)
    theta <- theta0
    chosen <- list(theta = theta, alpha = alpha, lambda = lam_max)
    for (lam in lambdas) {
      for (it in seq_len(max_iter)) {
        theta <- sum(w * bx * (by - alpha)) / sum(w * bx^2)
        z <- w * (by - theta * bx)
        a_new <- sign(z) * pmax(abs(z) - lam, 0) / w
        if (max(abs(a_new - alpha)) < 1e-10) {
          alpha <- a_new
          break
        }
        alpha <- a_new
      }
      q <- sum(w * (by - theta * bx - alpha)^2)
      df <- max(p - 1 - sum(alpha != 0), 1)
      chosen <- list(theta = theta, alpha = alpha, lambda = lam)
      ## heterogeneity stopping rule: accept once Q is no longer significant
      if (stats::pchisq(q, df = df, lower.tail = FALSE) > 0.05) break
    }
    chosen
  }
  fit <- fit_one(beta_X, beta_Y, se_Y)
  se <- NA_real_
  if (n_boot > 0L) {
    p <- length(beta_X)
    boots <- replicate(n_boot, {
      idx <- sample(p, replace = TRUE)
      fit_one(beta_X[idx], beta_Y[idx], se_Y[idx])$theta
    })
    se <- stats::sd(boots)
  }
  ci <- .ci_normal(fit$theta, se)
  list(theta = fit$theta, se = se,
       p_value = 2 * stats::pnorm(-abs(fit$theta / se)),
       ci_lo = ci[1L], ci_hi = ci[2L],
       lambda = fit$lambda, n_pleiotropic = sum(fit$alpha != 0))
}

## Pick one representative instrument per LD block (the strongest, by
## |beta_X| / se_Y), mimicking standard LD clumping applied before
## conventional MR analyses.
.prune_ld <- function(beta_X, se_Y, blocks) {
  vapply(blocks, function(blk) blk[which.max(abs(beta_X[blk]) / se_Y[blk])],
         integer(1L))
}
