\name{mr_ivw}
\alias{mr_ivw}
\title{Inverse-variance weighted (IVW) Mendelian randomization}
\usage{
mr_ivw(beta_X, beta_Y, se_Y)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of exposure and outcome association
estimates.}
\item{se_Y}{Numeric vector of standard errors of \code{beta_Y}.}
}
\value{
A list with \code{theta}, \code{se}, \code{p_value}, confidence interval
(\code{ci_lo}, \code{ci_hi}), Cochran's \code{Q} and the residual inflation
factor \code{phi}.
}
\description{
Multiplicative random-effects IVW estimator: the standard error is inflated
by the square root of the residual heterogeneity statistic when Cochran's Q
exceeds its degrees of freedom. Instruments are assumed independent.
Provided as a comparator for benchmarking ECC-MR.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
mr_ivw(sim$beta_X, sim$beta_Y, sim$se_Y)$theta
}
