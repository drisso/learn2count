#' Generate zero-inflated negative binomial data
#'
#' Simulate zinb data based on Karlis (or LPGM)
#' @param n number of samples
#' @param p number of variables (nodes)
#' @param mu mean of negative binomial part
#' @param theta dispersion parameter of negative binomial part, where \eqn{var=\mu+\mu^2/theta}
#' @param pi parameter relate to zero part and negative part
#' @param B adjacency matrix
#' @param mu.nois mean of noise
#' @importFrom countreg rzinbinom
zinb.simdata <- function(n,p,B, mu,mu.nois,theta,pi){
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
   ltri.sigma <- sigma[lower.tri(sigma)]
   nonzero.sigma <- ltri.sigma[which(ltri.sigma !=0 )]
   Y.mu <- c(rep(mu,nrow(sigma)), nonzero.sigma)
##data matrix
   Y <-  matrix(unlist(do.call(rbind, lapply(Y.mu,function(i) {rzinbinom(n,mu=i,theta, pi=pi)}))),length(Y.mu),n)
   X <- A%*%Y
   # add noise to all the nodes.
   X <- X + matrix(unlist(do.call(rbind, lapply(rep(mu.nois,p),function(i) { rzinbinom(n,mu=i,theta=1, pi=pi)}))),p,n)

   return(t(X))
}



