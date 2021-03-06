% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCnbinom.R
\name{nb.wald}
\alias{nb.wald}
\title{This function estimates adjacency matrix of a NB models given a matrix of counts, using glm.}
\usage{
nb.wald(X, maxcard, alpha, extend)
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
This function estimates adjacency matrix of a NB models given a matrix of counts, using glm.
}
