#' Structure learning with Poisson models
#'
#' This function estimates the adjacency matrix of a Poisson model given a
#' matrix of counts, using the glm function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients
pois.wald <- function(X,maxcard,alpha,extend){
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0

  while (ncard <= maxcard) {
    V <-  foreach(i = 1:p, .combine = "cbind") %dopar%{
      neighbor <- which (adj[, i] == 1)
      if (length(neighbor) >= ncard){
        condset <- combn(neighbor, ncard, FUN = list)
        for (j in 1: length(neighbor)){
          condset.temp <- condset
          indcond <- FALSE
          k <- 1
          while (!indcond & k <= length(condset.temp)){
            if (!(neighbor[j] %in% condset.temp[[k]])){
              fit <- glm(X[,i] ~ scale(X[,c(neighbor[j], condset.temp[[k]])]),
                         family="poisson")
              if (coefficients(summary(fit))[2,4] > alpha){
                adj[neighbor[j], i] <- 0
                indcond <- TRUE
              }
            }
            k <- k+1
          }
        }
      }
      return(adj[,i])
    }

    if (extend == TRUE){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else

      adj <- V * t(V)

    ncard <- ncard + 1

  }
  return(adj)
}

