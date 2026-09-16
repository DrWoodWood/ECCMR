\name{eccmr}
\alias{eccmr}
\title{ECC-MR: Error-correcting code inspired robust Mendelian randomization}
\usage{
eccmr(beta_X, beta_Y, se_Y, blocks = NULL, ld_mat = NULL,
      r2_thresh = 0.3, lambda = NULL, gamma = NULL,
      tuning = c("bic", "cv"),
      nlambda = 20L, lambda_min_ratio = 0.01,
      cv_folds = 10L, n_boot = 500L, snp_names = NULL,
      seed = NULL, max_iter = 200L, tol = 1e-08, verbose = FALSE)
}
\arguments{
\item{beta_X, beta_Y}{Numeric vectors of length p; exposure and outcome
association estimates for the genetic instruments.}
\item{se_Y}{Numeric vector of length p; standard errors of \code{beta_Y}.}
\item{blocks}{Optional list of integer vectors giving the LD block of each
instrument (e.g. external annotations). Required unless \code{ld_mat} is
given.}
\item{ld_mat}{Optional p x p matrix of pairwise LD correlations, used to
derive blocks with \code{\link{partition_ld_blocks}} and to threshold the
checks.}
\item{r2_thresh}{Numeric; r^2 threshold for edges/checks (default 0.3).}
\item{lambda, gamma}{Optional tuning parameters. \code{lambda} is the L1
decoding penalty (BIC/CV-selected when \code{NULL}). \code{gamma} is
deprecated since v0.3.0: the syndrome-consistency penalty is no longer part
of the core objective and defaults to 0; a positive value is still honoured
as an optional regulariser. When \code{NULL}
(the default), they are selected according to \code{tuning}.}
\item{tuning}{Character; tuning-parameter selection strategy. \code{"bic"}
(default) minimizes a BIC criterion along the \code{lambda}
path, which rewards actually fitting the pleiotropic effects; \code{"cv"}
uses K-fold cross-validation over LD blocks with the one-standard-error
rule.}
\item{nlambda}{Integer; size of the lambda grid (default 20).}
\item{lambda_min_ratio}{Numeric; smallest lambda as a fraction of the
largest (default 0.01).}
\value{
An object of class \code{eccmr}: a list with components
\describe{
  \item{theta}{Estimated causal effect from the relaxed (debiased) refit;
  \code{theta_l1} gives the L1-penalized estimate before debiasing.}
  \item{se, ci}{Bootstrap standard error and 95\% percentile confidence
  interval (when \code{n_boot > 0}).}
  \item{alpha}{Estimated pleiotropic effect of each instrument (relaxed
  refit; \code{alpha_l1} gives the penalized version used for support
  selection).}
  \item{pleiotropic}{Estimated nonzero pleiotropic effects (named when
  \code{snp_names} is supplied); \code{pleiotropic_index} gives the
  instrument indices.}
  \item{syndrome}{Result of \code{\link{syndrome_test}}: global pleiotropy
  test and standardized check statistics.}
  \item{lambda, gamma}{Selected tuning parameters.}
  \item{tune_table}{Data frame of BIC (or cross-validation) scores along
  the tuning path.}
  \item{boot_theta}{Bootstrap replicates of the causal estimate.}
  \item{n_instruments, n_checks}{Problem dimensions.}
  \item{fit}{The final relaxed \code{\link{ecc_fit}} result.}
  \item{call}{The matched call.}
}
}
\description{
Fits the ECC-MR model to GWAS summary statistics. Horizontal pleiotropy is
modelled as sparse corruption of SNP-level causal signals; a parity-check
matrix built from the LD structure annihilates the causal signal
(H beta_X = 0), so the syndrome H beta_Y isolates pleiotropic errors. The
causal effect and the pleiotropic effects are estimated jointly by convex
optimization; tuning parameters are selected by a BIC criterion along the
penalty path (or optionally by LD-block cross-validation), the selected
support is refitted without the L1 penalty to remove shrinkage bias, and
standard errors are obtained from an LD-block bootstrap. Unlike
MR-PRESSO-style approaches, pleiotropic instruments are corrected and
retained rather than discarded, and correlated instruments can be used
directly without LD pruning.
}
\examples{
sim <- simulate_eccmr(n_snps = 120, block_size = 10, scenario = "A",
                      seed = 42)
fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
             nlambda = 10, n_boot = 100, seed = 1)
print(fit)
summary(fit)
}
\references{
Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
Framework for Robust Mendelian Randomization. Manuscript.
}
