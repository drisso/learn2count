#' Structure learning with Poisson models using neighbourhood selection X_i|X_{pre(i)}
#'
#' This function estimates the adjacency matrix of a Poisson model given a
#' matrix of counts and topological ordering, using the glm function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param order the topological ordering of variables (names of nodes)
#' @return the estimated adjacency matrix of the graph.
#' @export neighborord.LP function
#' @importFrom stats coefficients

LPGM.ord <- function(X,order,alpha){
  p <- ncol(X)
  T <- foreach(i = 1:p, .combine="cbind",.packages = c("stats"),.export = c("neighborord.LP")) %dopar%{
    fit  <- neighborord.LP(i,X=X,order,alpha,p)
    return(fit)
  }
  Beta   <- matrix(unlist(T),p,p)
  return(Beta)
}

######### This function estimate parent set for each nodes using glm function
neighborord.LP <- function(i,X,order,alpha,p){
  Beta      <- matrix(0,p,p)
  if (length(colnames(X))>0){
    rownames(Beta) <- colnames(X)
  }else{
    rownames(Beta) <- as.character(seq(1,p,1))
  }
  temp <- which(order==colnames(X)[i])
  if(temp==1) {
    Beta[,i] <- 0
  }else{
    Y <- X[,order[1:(temp-1)]]
    fit <- glm(X[,i]~Y,family = "poisson")
    ind.edge <- which(coefficients(summary(fit))[-1,4] < alpha)
    Beta[order[ind.edge], i] <- 1
  }
  return(Beta[,i])
}
