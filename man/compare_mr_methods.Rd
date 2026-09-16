\name{compare_mr_methods}
\alias{compare_mr_methods}
\title{Compare ECC-MR with mainstream MR methods on one data set}
\usage{
compare_mr_methods(beta_X, beta_Y, se_Y, blocks = NULL, ld_mat = NULL,
                   r2_thresh = 0.3,
                   methods = c("ECC-MR", "IVW", "MR-Egger",
                               "Weighted median", "MR-Lasso"),
                   prune = TRUE, n_boot = 100L, median_boot = 200L,
                   lasso_boot = 200L, nlambda = 10L,
                   cv_folds = 5L, seed = NULL,
                   verbose = FALSE)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of exposure and outcome association
estimates.}
\item{se_Y}{Numeric vector of standard errors of \code{beta_Y}.}
\item{blocks}{Optional list of integer vectors of LD blocks; derived from
\code{ld_mat} when omitted.}
\item{ld_mat}{Optional p x p matrix of pairwise LD correlations.}
\item{r2_thresh}{Numeric; r^2 threshold used when deriving blocks from
\code{ld_mat}.}
\item{methods}{Character vector; any of \code{"ECC-MR"}, \code{"IVW"},
\code{"MR-Egger"}, \code{"Weighted median"}, \code{"MR-Lasso"}.}
\item{prune}{Logical; apply conventional methods to one LD-pruned
representative per block (default \code{TRUE}). If \code{FALSE}, they use
all instruments, violating their independence assumption.}
\item{n_boot}{Integer; bootstrap replicates for the ECC-MR standard error.}
\item{median_boot, lasso_boot}{Integer; bootstrap replicates for the
weighted-median and MR-Lasso standard errors.}
\item{nlambda, cv_folds}{Tuning controls passed to
\code{\link{eccmr}}.}
\item{seed}{Optional seed for reproducibility.}
\item{verbose}{Logical; print progress messages.}
}
\value{
A data frame with one row per method: \code{method}, \code{theta},
\code{se}, \code{ci_lo}, \code{ci_hi}, and the number of instruments used
(\code{n_instruments}).
}
\description{
Applies ECC-MR and a set of conventional summary-data MR methods (IVW,
MR-Egger, weighted median, MR-Lasso) to the same data set. By default the
conventional methods are applied to one LD-pruned representative instrument
per LD block (standard clumping practice), whereas ECC-MR uses all
instruments and models the LD structure directly -- this contrast is
precisely where ECC-MR gains its information advantage.
}
\examples{
sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                      seed = 42)
compare_mr_methods(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
                   n_boot = 50, seed = 1)
}
