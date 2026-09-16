\name{build_parity_check}
\alias{build_parity_check}
\title{Construct the LD parity-check matrix}
\usage{
build_parity_check(beta_X, blocks, r2_mat = NULL, tau = 0.3)
}
\arguments{
\item{beta_X}{Numeric vector of length p; exposure association estimates.}
\item{blocks}{A list of integer vectors giving the LD blocks, e.g. from
\code{\link{partition_ld_blocks}} or from external LD-block annotations
(Berisa & Pickrell, 2016).}
\item{r2_mat}{Optional p x p matrix of squared LD correlations; when given,
only within-block pairs with r^2 > \code{tau} generate a check. When
\code{NULL}, all within-block pairs are used.}
\item{tau}{Numeric; r^2 threshold applied when \code{r2_mat} is supplied.}
}
\value{
A list with components
\describe{
  \item{H}{Sparse (\code{dgCMatrix}) m x p parity-check matrix.}
  \item{edges}{m x 2 integer matrix of instrument index pairs (checks).}
  \item{check_block}{Integer vector of length m with the LD-block label of
  each check.}
}
}
\description{
Builds the parity-check matrix H of the ECC-MR framework. Each retained
within-block edge (i, j) defines one parity check (one row of H) with
entries H[k, i] = beta_X[j] and H[k, j] = -beta_X[i], so that
H \%*\% beta_X = 0 by construction: the causal signal lies in the null
space (code space) of H, whereas pleiotropic effects generate a nonzero
syndrome.
}
\examples{
sim <- simulate_eccmr(n_snps = 60, block_size = 6, seed = 1)
pc <- build_parity_check(sim$beta_X, sim$blocks)
max(abs(as.numeric(pc$H \%*\% sim$beta_X)))  # exactly 0
}
\references{
Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
Framework for Robust Mendelian Randomization. Manuscript.
}
