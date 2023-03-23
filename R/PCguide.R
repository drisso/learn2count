#' Structure learning with Poisson models using the Or-PPGM algorithm
#'
#' This function estimates the adjacency matrix of a Poisson model given a
#' matrix of counts and topological ordering, using the Or-PPGM algorithm
#' proposed in Nguyen et al. (2022).
#'
#' @references
#' Nguyen, Chiogna, Risso, Banzato (2022). Guided structure learning of
#' DAGs for count data. arXiv:2206.09754.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests.
#' @param maxcard the uper bound of the cardinality of the conditional sets K.
#' @param order the topological ordering of variables (names of nodes).
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients glm
#' @importFrom utils combn
pois.ord <- function(X, maxcard, alpha, order){

    p <- ncol(X)
    adj.ord <- matrix(0, p, p)
    adj.ord[upper.tri(adj.ord, diag = FALSE)] <- 1
    rownames(adj.ord) <- colnames(adj.ord) <- as.character(order)
    adj <- adj.ord[colnames(X), colnames(X)]

    ncard <- 0
    while (ncard <= maxcard) {
        V <-  foreach(i = 1:p, .combine = "cbind") %dopar% {

            neighbor <- which (adj[, i] == 1)

            if (length(neighbor) > ncard){
                for (j in seq_len(neighbor)){
                    neighbor.tem <- setdiff(neighbor, neighbor[j])

                    if (length(neighbor.tem)>= ncard){
                        condset.temp <- combn(neighbor.tem, ncard, FUN = list)
                        indcond <- FALSE
                        k <- 1
                        while (indcond == FALSE & k <= length(condset.temp)){
                            fit <- glm(X[,i] ~ X[,c(neighbor[j],
                                                    condset.temp[[k]])],
                                       family="poisson")
                            if (coefficients(summary(fit))[2,4] > alpha){
                                adj[neighbor[j], i] <- 0
                                indcond <- TRUE
                            }
                            k = k+1
                        }
                    }
                }
            }
            return(adj[,i])
        }

        adj <- V
        ncard <- ncard + 1
    }

    return(adj)
}
