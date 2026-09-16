\name{mr_egger}
\alias{mr_egger}
\title{MR-Egger regression}
\usage{
mr_egger(beta_X, beta_Y, se_Y)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of exposure and outcome association
estimates.}
\item{se_Y}{Numeric vector of standard errors of \code{beta_Y}.}
}
\value{
A list with the causal estimate \code{theta} (slope), its \code{se} and
\code{p_value}, confidence interval (\code{ci_lo}, \code{ci_hi}), and the
pleiotropy \code{intercept} with its standard error and p-value.
}
\description{
Weighted regression of \code{beta_Y} on \code{beta_X} with an unconstrained
intercept. The slope is a consistent causal estimate under the InSIDE
assumption; the intercept tests for directional pleiotropy. Instruments are
assumed independent. Provided as a comparator for benchmarking ECC-MR.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
mr_egger(sim$beta_X, sim$beta_Y, sim$se_Y)$theta
}
