\name{syndrome_test}
\alias{syndrome_test}
\title{Syndrome and global test of pleiotropy}
\usage{
syndrome_test(beta_Y, se_Y, H, check_block = NULL, sigma = NULL,
              tol_eig = 1e-10)
}
\arguments{
\item{beta_Y}{Numeric vector of length p; outcome association estimates.}
\item{se_Y}{Numeric vector of length p; standard errors of \code{beta_Y}.}
\item{H}{Sparse m x p parity-check matrix, from
\code{\link{build_parity_check}}.}
\item{check_block}{Optional integer vector of length m assigning each check
to an LD block; derived from \code{H} when omitted.}
\item{sigma}{Optional p x p error covariance matrix (use when the exposure
and outcome GWAS share participants); defaults to diag(se_Y^2).}
\item{tol_eig}{Numeric; relative eigenvalue tolerance for the
pseudoinverse.}
}
\value{
A list with components
\describe{
  \item{statistic}{The syndrome statistic T_syn.}
  \item{df}{Degrees of freedom, equal to the numerical rank of H.}
  \item{p_value}{P-value of the global pleiotropy test.}
  \item{syndrome}{The syndrome vector s.}
  \item{z}{Standardized check statistics s_k / sqrt(Omega_kk).}
}
}
\description{
Computes the syndrome s = H \%*\% beta_Y of the ECC-MR framework and the
associated global test of pleiotropy. Under the null hypothesis of no
horizontal pleiotropy, the syndrome equals H \%*\% epsilon and the quadratic
form T_syn = s' Omega^+ s (Moore-Penrose pseudoinverse) follows a
chi-squared distribution with rank(H) degrees of freedom. Standardized
check statistics localize the corrupted regions, providing an ECC-native
analogue of the MR-PRESSO global and outlier tests that remains valid for
correlated instruments.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, scenario = "A", seed = 1)
pc <- build_parity_check(sim$beta_X, sim$blocks)
syndrome_test(sim$beta_Y, sim$se_Y, pc$H, pc$check_block)
}
