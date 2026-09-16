\name{run_simulation}
\alias{run_simulation}
\title{Simulation benchmark of ECC-MR against mainstream MR methods}
\usage{
run_simulation(nsim = 20L, scenario = "A",
               methods = c("ECC-MR", "IVW", "MR-Egger",
                           "Weighted median", "MR-Lasso"),
               prune = TRUE, n_snps = 300L, block_size = 10L,
               theta = 0.3, ..., n_boot = 100L, nlambda = 10L,
               cv_folds = 5L, cores = 1L,
               seed = 1L, verbose = TRUE)
}
\arguments{
\item{nsim}{Integer; number of simulation replicates (default 20; use 50+
for stable Monte-Carlo error).}
\item{scenario}{Character; pleiotropy scenario passed to
\code{\link{simulate_eccmr}} (\code{"A"}, \code{"B"} or \code{"C"}).}
\item{methods}{Character vector of methods, see
\code{\link{compare_mr_methods}}.}
\item{prune}{Logical; LD-prune instruments for conventional methods
(default \code{TRUE}).}
\item{n_snps, block_size, theta, ...}{Simulation parameters passed to
\code{\link{simulate_eccmr}} (e.g. \code{pi_pleio}, \code{directional},
\code{delta}).}
\item{n_boot, nlambda, cv_folds}{ECC-MR controls for each
replicate; defaults are deliberately light to keep the benchmark fast.}
\item{cores}{Integer; number of parallel workers (Unix only, via
\code{parallel::mclapply}); \code{1} (default) runs serially.}
\item{seed}{Integer base seed; replicate \code{k} uses \code{seed + k - 1}.}
\item{verbose}{Logical; print progress messages.}
}
\value{
An object of class \code{eccmr_sim}: a list with
\describe{
  \item{summary}{Data frame with one row per method: mean estimate, bias,
  empirical SD, RMSE, mean reported SE, 95\% CI coverage and rejection
  rate.}
  \item{raw}{Data frame with one row per method per replicate.}
  \item{theta, scenario, nsim}{Benchmark meta-information.}
}
The \code{print} method displays a formatted comparison table (lowest RMSE
marked with \code{*}); the \code{plot} method draws bias and RMSE panels.
}
\description{
Repeatedly simulates GWAS summary statistics with
\code{\link{simulate_eccmr}}, applies ECC-MR and conventional MR methods to
each replicate via \code{\link{compare_mr_methods}}, and aggregates bias,
Monte-Carlo variability, RMSE, confidence-interval coverage and rejection
rates. This reproduces the simulation study of the ECC-MR manuscript and
demonstrates where ECC-MR outperforms the alternatives.
}
\examples{
\dontrun{
bench <- run_simulation(nsim = 50, scenario = "A", directional = 0.05,
                        seed = 1)
print(bench)
plot(bench)
}
}
