\name{run_boundary_scan}
\alias{run_boundary_scan}
\alias{print.eccmr_boundary}
\title{Scan the operating boundary of ECC-MR across pleiotropy regimes}
\usage{
run_boundary_scan(pis = c(0.1, 0.2, 0.3), directionals = c(0.05, 0),
                  tau_alphas = c(0.02, 0.04, 0.08), nsim = 10L,
                  n_snps = 300L, block_size = 10L, theta = 0.3,
                  se_y = 0.015,
                  methods = c("ECC-MR", "IVW", "MR-Egger",
                              "Weighted median", "MR-Lasso"),
                  prune = TRUE, nlambda = 12L,
                  cores = 1L, seed = 1L, verbose = TRUE)
}
\arguments{
\item{pis}{Numeric vector; fractions of pleiotropic instruments.}
\item{directionals}{Numeric vector; mean direct effects. Use 0 for
balanced pleiotropy and a positive value for directional pleiotropy.}
\item{tau_alphas}{Numeric vector; SD of the direct effects.}
\item{nsim}{Integer; Monte-Carlo replicates per grid cell (default 10; use
20+ for stable Monte-Carlo error).}
\item{n_snps, block_size, theta, se_y}{Passed to
\code{\link{simulate_eccmr}}.}
\item{methods}{Character vector; subset of \code{"ECC-MR"}, \code{"IVW"},
\code{"MR-Egger"}, \code{"Weighted median"}, \code{"MR-Lasso"}.}
\item{prune}{Logical; apply conventional methods to one LD-pruned
representative per block (default \code{TRUE}), while ECC-MR always uses
all instruments.}
\item{nlambda}{Tuning controls passed to \code{\link{eccmr}}.}
\item{cores}{Integer; parallel workers over replicates (Unix only).}
\item{seed}{Integer; base seed (each grid cell gets a derived seed).}
\item{verbose}{Logical; print per-cell progress.}
}
\value{
An object of class \code{eccmr_boundary}: a data frame with one row per
regime x method, columns \code{directional}, \code{pi}, \code{tau_alpha},
\code{method}, \code{bias}, \code{mc_sd}, \code{rmse}, \code{n_reps}. The
\code{print} method summarises each regime and marks the lowest-RMSE
method with \code{*}.
}
\description{
Runs a grid of simulation regimes varying (i) the fraction of pleiotropic
instruments, (ii) the mean direct effect (0 gives balanced pleiotropy) and
(iii) the per-variant pleiotropy scale, and reports bias / Monte-Carlo SD /
RMSE of ECC-MR and conventional summary-data MR methods in every regime.
This answers the question "when does modelling LD with parity checks
actually help": ECC-MR is expected to dominate whenever pleiotropic effects
are individually detectable (roughly |alpha| >~ 2-3 se_Y) or balanced, and
to fall back to IVW-like behaviour when pleiotropy is both frequent and
weaker than the noise floor. Only point estimates are computed (no
bootstrap), so the scan is fast.
}
\examples{
bd <- run_boundary_scan(pis = c(0.1, 0.3), directionals = c(0.05, 0),
                        tau_alphas = c(0.02, 0.08), nsim = 4)
print(bd)
}
