% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/QPtransform.R
\name{QPtransform}
\alias{QPtransform}
\alias{QPtransform,SummarizedExperiment-method}
\alias{QPtransform,matrix-method}
\title{Quantile matching and power transformation}
\usage{
QPtransform(x, ...)

\S4method{QPtransform}{SummarizedExperiment}(x, assay_from = "counts", assay_to = "processed", ...)

\S4method{QPtransform}{matrix}(x)
}
\arguments{
\item{x}{the matrix of counts (n times p) or a SummarizedExperiment
containing such matrix (transposed).}

\item{...}{Additional arguments (currently not used).}

\item{assay_from}{The assay with the input count data.}

\item{assay_to}{The assay in which to store the processed data.}
}
\value{
if x is a matrix, a matrix with the processed data; if x
  is a SummarizedExperiment, a SummarizedExperiment object with the processed
  matrix as an additional assay.
}
\description{
This function implements the preprocessing strategy discussed in Nguyen et
al. (2020). We recommend this transformation when applying the PCzinb
algorithm to real dataset.
}
\details{
Briefly, the transformation consists of two steps: (i) matching of the 95
percentile across cells to account for sequencing depth; (ii) adjusting the
data to be closer to a zinb distribution by using a power transformation
\eqn{X^\alpha}, where \eqn{\alpha \in [0,1]} is chosen to minimize the
Kolmogorov-Smirnov statistic.
}
\examples{
library(SummarizedExperiment)
se <- SummarizedExperiment(assays=list(counts=matrix(rpois(50, 5), ncol=10)))
suppressWarnings(se <- QPtransform(se))
se
mat <- matrix(rpois(50, 5), nrow=10)
suppressWarnings(QPtransform(mat))
}
\references{
Nguyen, T. K. H., Berge, K. V. D., Chiogna, M., & Risso, D.
  (2020). Structure learning for zero-inflated counts, with an application to
  single-cell RNA sequencing data. arXiv:2011.12044.
}
