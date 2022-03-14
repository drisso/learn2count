#' Structure learning for count data
#'
#' This method implements a set of local algorithms to estimate the structure of
#' a graph, starting from a matrix (or a SummarizedExperiment) of counts,
#' typically representing the gene expression of p genes in n cells.
#'
#' This approach is based on the PC algorithm and assumes that the
#' node-conditional distribution of the data is zero-inflated negative binomial
#' (`zinb`), or its special cases negative binomial (`nb`) and Poisson (`poi`).
#'
#' There are two version of the algorithm for `zinb`: one that assumes that the
#' structure of the graph depends only on the mean parameter (`zinb0`) and one
#' that assumes that the structure of the graph depends on both the mean
#' parameter and the zero inflation parameter (`zinb1`).
#'
#' Obviously, `zinb1` is the most general case and should be the more accurate
#' in most cases. However, the `poi` and `nb` algorithms offer some
#' computational advantages.
#'
#' The default value of the significance level alpha is based on a rule that
#' ensure control of the type I error of all tests. See Nguyen et al. (2020) for
#' details.
#'
#' @references Nguyen, T. K. H., Berge, K. V. D., Chiogna, M., & Risso, D.
#'   (2020). Structure learning for zero-inflated counts, with an application to
#'   single-cell RNA sequencing data. arXiv:2011.12044.
#'
#' @import methods
#' @rdname PCzinb
setGeneric(
    name = "PCzinb",
    def = function(x, ...) standardGeneric("PCzinb")
)

#' Quantile matching and power transformation
#'
#' This function implements the preprocessing strategy discussed in Nguyen et
#' al. (2020). We recommend this transformation when applying the PCzinb
#' algorithm to real dataset.
#'
#' Briefly, the transformation consists of two steps: (i) matching of the 95
#' percentile across cells to account for sequencing depth; (ii) adjusting the
#' data to be closer to a zinb distribution by using a power transformation
#' \eqn{X^\alpha}, where \eqn{\alpha \in [0,1]} is chosen to minimize the
#' Kolmogorov-Smirnov statistic.
#'
#' @param x the matrix of counts (n times p) or a SummarizedExperiment
#'   containing such matrix (transposed).
#' @param ... Additional arguments (currently not used).
#' @references Nguyen, T. K. H., Berge, K. V. D., Chiogna, M., & Risso, D.
#'   (2020). Structure learning for zero-inflated counts, with an application to
#'   single-cell RNA sequencing data. arXiv:2011.12044.
#'
#' @rdname QPtransform
setGeneric(
    name = "QPtransform",
    def = function(x, ...) standardGeneric("QPtransform")
)
