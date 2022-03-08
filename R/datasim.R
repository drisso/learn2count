#' Generate Poisson, negative binomial (NB), and zero-inflated NB (ZINB) data
#'
#' This function generates a data matrix (\eqn{n \times p}) under Poisson, NB,
#' or ZINB assumptions, based on Karlis (or LPGM. 2013) with adjancency matrix B.
#' @param n number of samples
#' @param p number of variables (nodes)
#' @param B adjacency matrix
#' @param family the distribution from which to simulate (Poisson, NB, ZINB)
#' @param mu mean of the non-zero component
#' @param mu_noise mean of noise
#' @param theta dispersion parameter of the non-zero component,
#'  where \eqn{var=\mu+\mu^2/theta}
#' @param pi probability of zero inflation
#' @importFrom countreg rzinbinom
#' @importFrom stats rpois rnbinom
#' @importFrom parallel mclapply
#' @examples
#' p <- 4
#' B <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0), nrow=p)
#' n <- 10
#' simdata(n, p, B, mu = 5, mu_noise=1)
#' simdata(n, p, B, mu = 5, mu_noise=1, family="NB", theta=1)
#' simdata(n, p, B, mu = 5, mu_noise=1, family="ZINB", theta=1, pi=0.2)
#' @export
simdata <- function(n, p, B, family = c("Poisson", "NB", "ZINB"), mu, mu_noise,
                    theta=NA, pi=NA){

    family <- match.arg(family)

    if(family == "Poisson" & !is.na(theta)) {
        warning("The value of theta is ignored when using a Poisson.")
    } else if(family != "ZINB" & !is.na(pi)) {
        warning("The value of pi is ignored when not using a ZINB.")
    }

    if(family != "Poisson" & is.na(theta)) {
        stop("Theta must be specified when using NB or ZINB.")
    }

    if(family == "ZINB" & is.na(pi)) {
        stop("Pi must be specificied when using ZINB.")
    }

    # create "adjacency" matrix A from the adjacency matrix B
    if(!isSymmetric.matrix(B)){
        stop("B is not a symmetric matrix")
    }

    if(nrow(B) != p){
        stop("The adjacency matrix should have dimension equal to p.")
    }

    A <- diag(nrow(B))
    for(i in seq_len(nrow(B)-1)){
        for( j in seq(i+1, ncol(B))){
            if(B[i,j] == 1){
                tmp <- rep(0,nrow(B))
                tmp[c(i,j)] <- 1
                A <- cbind(A,tmp)
            }
        }
    }

    ## vector mean
    sigma <- mu*B
    ltri_sigma <- sigma[lower.tri(sigma)]
    nonzero_sigma <- ltri_sigma[which(ltri_sigma !=0 )]
    Y_mu <- c(rep(mu, nrow(sigma)), nonzero_sigma)

    ##data matrix
    Y <-  matrix(unlist(do.call(rbind, lapply(Y_mu, generate_data, theta=theta, pi=pi, family=family, n=n))),length(Y_mu),n)
    X <- A%*%Y
    # add noise to all the nodes.
    theta_noise <- ifelse(is.na(theta), NA, 1)
    X <- X + matrix(unlist(do.call(rbind, lapply(rep(mu_noise,p), generate_data, theta=theta_noise, pi=pi, family=family, n=n))),p,n)

    return(t(X))

}

generate_data <- function(mu, theta=NA, pi=NA, family, n) {
    switch(family,
           Poisson = rpois(n, lambda=mu),
           NB = rnbinom(n,mu=mu,size=theta),
           ZINB = rzinbinom(n, mu=mu, theta=theta, pi=pi))
}
