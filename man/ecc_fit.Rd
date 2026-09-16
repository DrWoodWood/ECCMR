\name{ecc_fit}
\alias{ecc_fit}
\title{Fit the ECC-MR model for fixed tuning parameters}
\usage{
ecc_fit(beta_X, beta_Y, se_Y, H, lambda, gamma, check_block = NULL,
        sigma = NULL, theta_init = NULL, alpha_init = NULL,
        max_iter = 200L, tol = 1e-08, tol_eig = 1e-10)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of length p; exposure and outcome
association estimates.}
\item{se_Y}{Numeric vector of length p; standard errors of \code{beta_Y}.}
\item{H}{Sparse m x p parity-check matrix, from
\code{\link{build_parity_check}}.}
\item{lambda}{Numeric; sparsity (L1) tuning parameter.}
\item{gamma}{Numeric; syndrome-consistency tuning parameter.}
\item{check_block}{Optional integer vector assigning each check to an LD
block; derived from \code{H} when omitted.}
\item{sigma}{Optional p x p error covariance matrix (use when the exposure
and outcome GWAS overlap); defaults to diag(se_Y^2).}
\item{theta_init, alpha_init}{Optional starting values; default to the IVW
estimate and zero.}
\item{max_iter}{Integer; maximum number of outer iterations.}
\item{tol}{Numeric; convergence tolerance on the relative objective change.}
\item{tol_eig}{Numeric; relative eigenvalue tolerance for the
pseudoinverse.}
}
\value{
A list with components
\describe{
  \item{theta}{Estimated causal effect.}
  \item{alpha}{Estimated pleiotropic effects (length p).}
  \item{iter}{Number of outer iterations used.}
  \item{converged}{Logical; whether the tolerance was reached.}
  \item{objective}{Final objective value.}
}
}
\description{
Minimizes the convex ECC-MR objective
0.5 * sum_i w_i (beta_Yi - theta beta_Xi - alpha_i)^2
+ lambda * ||alpha||_1 + (gamma / 2) * ||s - H alpha||^2_{Omega^+}
by block coordinate descent: a closed-form inverse-variance weighted update
for the causal effect theta alternates with cyclic soft-thresholding
coordinate-descent sweeps for the pleiotropic effects alpha. Joint
convexity guarantees convergence to the global minimum. One coordinate
sweep costs O(number of checks), so decoding remains feasible for large
instrument sets.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
pc <- build_parity_check(sim$beta_X, sim$blocks)
fit <- ecc_fit(sim$beta_X, sim$beta_Y, sim$se_Y, pc$H,
               lambda = 5, gamma = 1)
fit$theta
}
