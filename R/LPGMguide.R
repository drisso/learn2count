#' Structure learning with Poisson models using the Or-LPGM algorithm
#'
#' This function estimates the adjacency matrix of a Local Poisson Graphical
#' Model (LPGM) given a matrix of counts and topological ordering, using
#' generalized linear models. It is a natural extension of the LPGM of Allen and
#' Liu (2013). See Nguyen et al. (2022) for details.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests.
#' @param order the topological ordering of variables (names of nodes).
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients glm
#' @references
#' Allen and Liu (2013). A local Poisson graphical model for inferring
#' networksfrom sequencing data. IEEE Transactions on NanoBioscience, 12(3),
#' 189–198.
#' @references
#' Nguyen, Chiogna, Risso, Banzato (2022). Guided structure learning of
#' DAGs for count data. arXiv:2206.09754.
LPGM.ord <- function(X,order,alpha){

    p <- ncol(X)

    Tfit <- foreach(i = 1:p, .combine="cbind", .packages = c("stats"),
                    .export = c("neighborord.LP")) %dopar% {
                        fit  <- neighborord.LP(i, X=X, order, alpha,p)
                        return(fit)
                    }

    Beta   <- matrix(unlist(Tfit), p, p)

    return(Beta)
}

######### This function estimate parent set for each nodes using glm
neighborord.LP <- function(i, X, order, alpha, p){

    Beta <- matrix(0,p,p)

    if (length(colnames(X))>0) {
        rownames(Beta) <- colnames(X)
    } else {
        rownames(Beta) <- as.character(seq(1, p, 1))
    }

    temp <- which(order == colnames(X)[i])

    if(temp == 1) {
        Beta[,i] <- 0
    } else {
        Y <- X[,order[seq_len(temp-1)]]
        fit <- glm(X[,i] ~ Y, family = "poisson")
        ind.edge <- which(coefficients(summary(fit))[-1,4] < alpha)
        Beta[order[ind.edge], i] <- 1
    }
    return(Beta[,i])
}
