% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCzinb0noT.R
\name{zinb0.noT}
\alias{zinb0.noT}
\title{Structure learning with zero-inflated negative binomial model (mean only)}
\usage{
zinb0.noT(X, maxcard, alpha, extend)
}
\arguments{
\item{X}{the matrix of counts (n times p).}

\item{maxcard}{the uper bound of the cardinality of the conditional sets K}

\item{alpha}{the significant level of the tests}

\item{extend}{if TRUE it considers the union of the tests, otherwise it
considers the intersection.}
}
\value{
the estimated adjacency matrix of the graph.
}
\description{
This function estimates the adjacency matrix of a ZINB model given a matrix
of counts, using the optim function.
}
\details{
This approach assumes that the structure of the graph only depends on the
mean parameter, treating zero inflation as a technical noise effect. We call
this model `zinb0`.
}
