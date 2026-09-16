\name{simulate_eccmr}
\alias{simulate_eccmr}
\title{Simulate GWAS summary statistics with LD-block structure}
\usage{
simulate_eccmr(n_snps = 300L, block_size = 10L, theta = 0.3,
               se_y = 0.015, scenario = c("A", "B", "C"),
               pi_pleio = 0.3, pi_block = 0.2, directional = 0.03,
               tau_alpha = 0.04, delta = 0.02, seed = NULL)
}
\arguments{
\item{n_snps}{Integer; number of instruments.}
\item{block_size}{Integer; instruments per LD block.}
\item{theta}{Numeric; true causal effect.}
\item{se_y}{Numeric; standard error of the outcome associations.}
\item{scenario}{Character; one of \code{"A"} (independent directional
pleiotropy; IVW biased), \code{"B"} (LD-block clustered pleiotropy with
random sign; balanced) or \code{"C"} (LD-block clustered, directional,
individually weak pleiotropy).}
\item{pi_pleio, pi_block, directional, tau_alpha, delta}{Scenario
parameters.}
\item{seed}{Optional random seed.}
}
\value{
A list with \code{beta_X}, \code{beta_Y}, \code{se_Y}, \code{alpha} (true
pleiotropic effects), \code{blocks} (LD-block membership) and the true
\code{theta}.
}
\description{
Generates synthetic summary statistics mimicking the simulation scenarios
of the ECC-MR manuscript. Instruments are organised in equally sized LD
blocks.
}
\examples{
sim <- simulate_eccmr(n_snps = 100, block_size = 10, scenario = "A",
                      seed = 1)
str(sim)
}
