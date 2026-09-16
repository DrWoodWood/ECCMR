\name{eccmr-package}
\alias{eccmr-package}
\alias{_PACKAGE}
\title{eccmr: Error-Correcting Code Inspired Robust Mendelian Randomization}
\description{
ECC-MR reframes horizontal pleiotropy as sparse, correctable corruption of
SNP-level causal signals and leverages linkage disequilibrium (LD) among
instruments as structured redundancy. A parity-check matrix built from the
LD graph annihilates the causal signal, so the syndrome of the observed
outcome associations isolates pleiotropic errors; sparse syndrome decoding
within a convex optimization framework jointly estimates the causal effect
and the pleiotropic effects without discarding instruments.
}
\details{
The main entry point is \code{\link{eccmr}}. Lower-level building blocks are
\code{\link{partition_ld_blocks}}, \code{\link{build_parity_check}},
\code{\link{syndrome_test}} and \code{\link{ecc_fit}}.
\code{\link{simulate_eccmr}} generates synthetic summary statistics for
examples and testing; \code{\link{compare_mr_methods}} and
\code{\link{run_simulation}} benchmark ECC-MR against mainstream MR methods
(IVW, MR-Egger, weighted median, MR-Lasso) on simulated data.
}
\references{
Jiang J, Hu D, Zhang Q, Lin Z. ECC-MR: An Error-Correcting Code Inspired
Framework for Robust Mendelian Randomization. Manuscript.
}
\keyword{package}
