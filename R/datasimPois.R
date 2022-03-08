#' Generate Poisson data
#'
#' Generate data matrix (nxp) under Poisson assumption based on Karlis (or LPGM. 2013) with adjancency matrix B,
#' @param n number of samples
#' @param p number of variables (nodes)
#' @param lambda mean of Poisson
#' @param lambda.c mean of noise
#' @param B adjacency matrix
#' @importFrom stats rpois
pois.simdata <- function(n, p,B,lambda,lambda.c){

  # create "adjacency" matrix A from the adjacency matrix B
  if(nrow(B) != ncol(B)){
    print("not a symmetric matrix")
    return()
  }
  A <- diag(1,nrow=nrow(B),ncol=ncol(B))
  for(i in 1:(nrow(B)-1)){
    for( j in (i+1):ncol(B)){
       if(B[i,j] == 1){
        tmp <- rep(0,nrow(B))
        tmp[c(i,j)] <- 1
        A <- cbind(A,tmp)
      }

    }
  }
  ### vector mean
    sigma <- lambda*B
    ltri.sigma <- sigma[lower.tri(sigma)]
    nonzero.sigma <- ltri.sigma[which(ltri.sigma !=0 )]
    Y.lambda <- c(rep(lambda,nrow(sigma)), nonzero.sigma)
  ### data matrix X
    Y <- matrix(unlist( do.call(rbind, lapply(Y.lambda,function(i) { rpois(n,i)}))),length(Y.lambda),n)
    X <- A%*%Y
    # add the labmda.c to all the nodes.
    X <- X + matrix(unlist(do.call(rbind, lapply(rep(lambda.c,p),function(i) { rpois(n,i)}))) ,p,n)

  return(t(X))
}