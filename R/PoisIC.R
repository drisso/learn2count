
#' Structure learning with Poisson models using K2 algorithm
#'
#' This function finds the best fitting structure of a Poisson model given a
#' matrix of counts and topological ordering, using a given criterion ("AIC", "BIC").
#'
#' @param X the matrix of counts (n times p).
#' @param maxcard the uper bound of the cardinality of the parent sets 
#' @param order the topological ordering of variables (names of nodes)
#' @param criterion the score function that measure the fitting of structures, could be "AIC" or "BIC"
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients
Poisk2 <- function(X, order, criterion="BIC", maxcard) {
  if (length(colnames(X)) >0){ 
    nodes <- colnames(X)
  }else{
    nodes <- seq(1,dim(X)[2],1)
  }
  p <- length(nodes)
  pa_list <- list()
  
  ###### first auxiliary function to calculate the score
  f <- function(i,pa){
    Y <- X[,pa]
    q <- nrow(as.matrix(unique(Y)))
    if (q==0){
      fit <- glm(X[,i]~1,family="poisson")
    }else{
      fit <- glm(X[,i]~Y,family="poisson")
    }
    if (criterion =="BIC"){
      res <- -BIC(fit)
    }
    if(criterion =="AIC"){
      res <- -AIC(fit)
    }
    return (res)
  }
  
  
  ###### second auxiliary function to find new candidates for parent sets
  findmax <- function(i,parents) {
    gmax <- f(i, parents)
    z <- integer()
    if (pos==1) {return(z)} else{
      candidates <- setdiff(order[1:(pos-1)], parents)
      for (j in 1:length(candidates)) {
        pa <- c(parents,candidates[j])
        gnew <- f(i,pa)
        if (gnew > gmax) {gmax <-gnew
        z <- candidates[j]
        }}
      return(z)}}
  
  
  ### estimate the adjacency matrix
  Adj <- matrix(0,nrow = p, ncol = p)
  colnames( Adj) <- rownames( Adj)<- nodes
  Adj.est <- foreach(i = 1:p, .combine = "cbind") %dopar%{
    pa_list[[i]] <- integer()
    pos <- which( order==nodes[i])
    gold <- f(i, c())
    OK <- TRUE
    counter <- 0
    
    while ((OK) & (length(pa_list[[i]]) < min(maxcard,pos-1)) ){
      counter <- counter+1
      z <- findmax(i, pa_list[[i]])
      if (length(z) == 0) {
        OK <- FALSE
      }
      else {
        pa_list[[i]] <- c(pa_list[[i]],z)
      }
    }
    
    Adj[ pa_list[[i]],i] <- 1
    return(Adj[,i])
  }
  
  colnames(Adj.est) <- rownames(Adj.est) <- nodes
  list(adjm = Adj.est, graphN = as(Adj.est, "graphNEL"))
}




