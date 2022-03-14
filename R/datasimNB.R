#' Generate negatibe binomial (NB) data
#'
#' Generatedata matrix (nxp) under NB assumption based on Karlis (or LPGM. 2013) with adjancency matrix B.
#' @param n number of samples
#' @param p number of variables (nodes)
#' @param mu mean of negative binomial
#' @param mu.nois mean of noise
#' @param theta dispersion parameter of negative binomial, where \eqn{var=\mu+\mu^2/theta}
#' @param B adjacency matrix
nbinom.Simdata <- function(n, p,B,mu,mu.nois,theta){
  set.seed(123)
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

 ## vector mean
    sigma <- mu*B
    ltri.sigma <-sigma[lower.tri(sigma)]
    nonzero.sigma <- ltri.sigma[which(ltri.sigma !=0 )]
    Y.mu <- c(rep(mu,nrow(sigma)), nonzero.sigma)
    ## data matrix X
    Y <-  matrix(unlist(do.call(rbind, lapply(Y.mu,function(i) {rnbinom(n,mu=i,theta)}))),length(Y.mu),n)
    X <- A%*%Y
    # add the labmda.c to all the nodes.
    X <- X + matrix(unlist(do.call(rbind, lapply(rep(mu.nois,p),function(i) {rnbinom(n,mu=i,theta)}))),p,n)

    return(t(X))
}
