% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCzinb0noT.R
\name{zinb0.noT}
\alias{zinb0.noT}
\title{This function estimates adjacency matrix of a zinb0 models given a matrix of counts, using optim.
zinb0 models are zinb models with the structure is infered from mu (only)}
\usage{
zinb0.noT(X, maxcard, alpha, extend)
}
\arguments{
\item{X}{the matrix of counts.}

\item{maxcard}{the uper bound for cardinalities of conditional sets K}

\item{alpha}{the sisnificant level tests}

\item{extend}{TRUE if we consider the union of tests}
}
\value{
the adj matrix.
}
\description{
This function estimates adjacency matrix of a zinb0 models given a matrix of counts, using optim.
zinb0 models are zinb models with the structure is infered from mu (only)
}
