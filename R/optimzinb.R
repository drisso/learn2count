# Find a single dispersion parameter for a count by 1-dimensional optimization of the likelihood
# Given a vector of count, this function computes a single dispersion parameter
# (log(theta)) the counts under a zero-inflated negative binomial
#' (ZINB) model. The ZINB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param Y the vector of counts
#' @param mu the vector mean  of the negative binomial
#' @param theta is the vector dispersion parameter of the negative binomial.
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
#' @param logitPi the vector of logit of the probabilities of the zero component
#' @importFrom zinbwave zinb.loglik.dispersion zinb.loglik zinb.loglik.dispersion.gradient
#' @examples
#' n <- 10
#' mu <- seq(10,50,length.out=n)
#' logitPi <- rnorm(1)
#' zeta <- rnorm(1)
#' Y <- rnbinom(n=n, size=exp(zeta), mu=mu)
#' zinbOptimizeDispersion ( mu, logitPi,Y,n)
zinbOptimizeDispersion <- function( mu, logitPi,Y,n) {

  g <- optimize(f=zinb.loglik.dispersion, Y=Y, mu=mu,
                logitPi=logitPi, maximum=TRUE,interval=c(-100,100))

  zeta.op <- g$maximum

  zeta.ot <- try(optim( par=zeta.op, fn=zinb.loglik.dispersion ,
                        gr=zinb.loglik.dispersion.gradient,mu=mu,
                        logitPi=logitPi,Y=Y,control=list(fnscale=-1,trace=0),
                        method="BFGS")$par,silent = TRUE)
  if (class(zeta.ot) != "try-error"){
    zeta <- zeta.ot
  }else{
    zeta <- zeta.op
  }


  zeta <- rep((zeta),n)
  zeta
}

# Copied on 5/14/2019 from log1pexp of the copula package (v. 0.999.19) by
# Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan,
# Johanna G. Neslehova
# Copied here to avoid dependence on gsl which causes troubles.
log1pexp <- function (x, c0 = -37, c1 = 18, c2 = 33.3)
{
  if (has.na <- any(ina <- is.na(x))) {
    y <- x
    x <- x[ok <- !ina]
  }
  r <- exp(x)
  if (any(i <- c0 < x & (i1 <- x <= c1)))
    r[i] <- log1p(r[i])
  if (any(i <- !i1 & (i2 <- x <= c2)))
    r[i] <- x[i] + 1/r[i]
  if (any(i3 <- !i2))
    r[i3] <- x[i3]
  if (has.na) {
    y[ok] <- r
    y
  }
  else r
}




#' Parse ZINB regression model
#'
#' Given the parameters of a ZINB regression model, this function parses the
#' model and computes the vector of log(mu), logit(pi), and the dimensions of
#' the different components of the vector of parameters. See
#' \code{\link{zinb.loglik.regression}} for details of the ZINB regression model
#' and its parameters.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi) concatenated
#' @param A.mu matrix of the model (1,X_K\backslash\{s\},T, default=empty)
#' @param A.pi matrix of the model (1,X_K\backslash\{s\},T, default=empty)
#' @return A list with slots \code{logMu}, \code{logitPi}, \code{dim.alpha} (a
#'   vector of length 2 with the dimension of each of the vectors \code{a.mu},
#'   \code{a.pi}  in \code{alpha}), and \code{start.alpha} (a vector
#'   of length 2 with the starting indices of the 2 vectors in \code{alpha})
#' @seealso \code{\link{zinb.loglik.regression}}
zinb.regression.parseModel <- function(alpha, A.mu,A.pi) {

  n <- nrow(A.mu)
  logMu <-0
  logitPi <- 0
  dim.alpha <- rep(0,2)
  start.alpha <- rep(NA,2)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
    start.alpha[1] <- i+1
    i <- i+j
  }

  j <- ncol(A.pi)
  if (j>0) {
    logitPi <- logitPi - A.pi %*% alpha[(i+1):(i+j)]
    dim.alpha[2] <- j
    start.alpha[2] <- i+1
    i <- i+j
  }


  return(list(logMu=logMu, logitPi=logitPi, dim.alpha=dim.alpha,
              start.alpha=start.alpha))
}

#####################
#########optimal functions

###zinb1: there are structures both in $mu$ and in $pi$
#' @export
#' @importFrom zinbwave zinb.loglik.regression zinb.loglik.regression.gradient
optim_funnoT <- function(beta_mu, gamma_pi, Y, X_mu, zeta, n) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(beta_mu, gamma_pi),
           Y=Y, A.mu=cbind(rep(1,n),X_mu),
           A.pi=cbind(rep(1,n),X_mu),
           C.theta=matrix(zeta, nrow = n, ncol = 1),
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}
###zinb0: there are only structures in $mu$
#' @export
optim_fun0noT <- function(beta_mu, gamma_pi, Y, X_mu, zeta, n) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(beta_mu, gamma_pi),
           Y=Y, A.mu=cbind(rep(1,n),X_mu),
           A.pi= matrix(rep(1,n),n,1),
           C.theta=matrix(zeta, nrow = n, ncol = 1),
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}

