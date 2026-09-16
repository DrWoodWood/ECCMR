\name{mr_weighted_median}
\alias{mr_weighted_median}
\title{Weighted median Mendelian randomization}
\usage{
mr_weighted_median(beta_X, beta_Y, se_Y, n_boot = 200L)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of exposure and outcome association
estimates.}
\item{se_Y}{Numeric vector of standard errors of \code{beta_Y}.}
\item{n_boot}{Integer; number of bootstrap replicates for the standard
error (default 200).}
}
\value{
A list with \code{theta}, \code{se}, \code{p_value} and confidence interval
(\code{ci_lo}, \code{ci_hi}).
}
\description{
Weighted median of the SNP-specific ratio estimates with weights
\code{beta_X^2 / se_Y^2}. Consistent when up to half of the weight comes
from invalid instruments. The standard error is obtained by bootstrapping
instruments. Provided as a comparator for benchmarking ECC-MR.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
mr_weighted_median(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 50)$theta
}
