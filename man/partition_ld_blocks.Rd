\name{partition_ld_blocks}
\alias{partition_ld_blocks}
\title{Partition instruments into LD blocks}
\usage{
partition_ld_blocks(ld_mat, r2_thresh = 0.3)
}
\arguments{
\item{ld_mat}{A p x p matrix (or sparse \code{Matrix}) of pairwise LD
correlations r between instruments, estimated from a reference panel of
matched ancestry.}
\item{r2_thresh}{Numeric; edges connect instruments with r^2 >
\code{r2_thresh}.}
}
\value{
A list of integer vectors; each vector contains the instrument indices of
one LD block (singletons allowed).
}
\description{
Builds approximately independent LD blocks from a pairwise LD (correlation)
matrix by thresholding squared correlations and taking connected components
of the resulting graph.
}
\examples{
ld <- diag(6); ld[1, 2] <- ld[2, 1] <- 0.9; ld[3, 4] <- ld[4, 3] <- 0.8
partition_ld_blocks(ld, r2_thresh = 0.3)
}
