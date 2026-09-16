\name{mr_lasso}
\alias{mr_lasso}
\title{MR-Lasso}
\usage{
mr_lasso(beta_X, beta_Y, se_Y, nlambda = 25L, max_iter = 100L,
         n_boot = 200L)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of exposure and outcome association
estimates.}
\item{se_Y}{Numeric vector of standard errors of \code{beta_Y}.}
\item{nlambda}{Integer; size of the lambda grid (default 25).}
\item{max_iter}{Integer; maximum coordinate-descent iterations per lambda.}
\item{n_boot}{Integer; number of bootstrap replicates for the standard
error (default 200).}
}
\value{
A list with \code{theta}, \code{se}, \code{p_value}, confidence interval
(\code{ci_lo}, \code{ci_hi}), the selected \code{lambda} and the number of
instruments flagged as pleiotropic (\code{n_pleiotropic}).
}
\description{
Lasso-penalized IVW: an L1 penalty on the pleiotropic effects identifies
and removes invalid instruments. The penalty parameter is chosen by the
heterogeneity stopping rule of the original publication: decrease lambda
from its maximal value until Cochran's Q computed on the corrected
associations is no longer significant at the 5\% level. The standard error
is obtained by bootstrapping instruments. Provided as a comparator for
benchmarking ECC-MR.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
mr_lasso(sim$beta_X, sim$beta_Y, sim$se_Y, n_boot = 50)$theta
}
