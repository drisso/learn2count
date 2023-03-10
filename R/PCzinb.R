#' @rdname PCzinb
#' @importFrom SummarizedExperiment assay assayNames `assay<-`
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata `metadata<-`
#' @param ... Arguments to pass to the matrix method.
#' @param whichAssay The assay to use as input for the matrix method.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(matrix(rpois(50, 5), ncol=10))
#' PCzinb(se, method="poi")
setMethod(
    f = "PCzinb",
    signature = signature(x = "SummarizedExperiment"),
    definition = function(x, whichAssay = "processed", ...){
        if(whichAssay == "processed" & !whichAssay %in% assayNames(x)) {
            warning("We recommend to use QPtransform() before learning the graph.")
            whichAssay <- 1
        }
        adj <- PCzinb(t(assay(x, whichAssay)), ...)
        metadata(x)$adj_mat <- as(adj, "dgCMatrix")
        rownames(metadata(x)$adj_mat) <- rownames(x)
        colnames(metadata(x)$adj_mat) <- rownames(x)

        return(x)
})

#' @param x the matrix of counts (n times p) or a SummarizedExperiment containing such matrix (transposed).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @param method the algorithm used to estimate the graph: `poi`, `nb`, `zinb0`,
#'   or `zinb1`. See details below.
#' @return if x is a matrix, the estimated adjacency matrix of the graph; if x
#'   is a SummarizedExperiment, a SummarizedExperiment object with the adjacency
#'   matrix as metadata.
#' @rdname PCzinb
#' @importFrom stats pnorm
#' @export
#' @examples
#' mat <- matrix(rpois(50, 5), nrow=10)
#' PCzinb(mat, method="poi")
setMethod(
    f = "PCzinb",
    signature = signature(x ="matrix"),
    definition = function(x,
                          method=c("poi", "nb", "zinb0", "zinb1"),
                          alpha=2*pnorm(nrow(x)^.2,lower.tail=FALSE),
                          maxcard=2,
                          extend=TRUE) {

        method <- match.arg(method)

        switch(method,
               poi = pois.wald(x, maxcard, alpha, extend),
               nb = nb.wald(x, maxcard, alpha, extend),
               zinb0 = zinb0.noT(x, maxcard, alpha, extend),
               zinb1 = zinb1.noT(x, maxcard, alpha, extend))

})

#' @rdname PCzinb
#' @importFrom SingleCellExperiment `rowPair<-`
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @param ... Arguments to pass to the matrix method.
#' @param whichAssay The assay to use as input for the matrix method.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' se <- SingleCellExperiment(matrix(rpois(50, 5), ncol=10))
#' se <- PCzinb(se, method="poi")
#' rowPair(se)
setMethod(
    f = "PCzinb",
    signature = signature(x = "SingleCellExperiment"),
    definition = function(x, whichAssay = "processed", ...){
        if(whichAssay == "processed" & !whichAssay %in% assayNames(x)) {
            warning("We recommend to use QPtransform() before learning the graph.")
            whichAssay <- 1
        }
        adj <- PCzinb(t(assay(x, whichAssay)), ...)
        rowPair(x) <- as(adj, "dgCMatrix")

        return(x)
})
