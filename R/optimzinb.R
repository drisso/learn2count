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



#' Log-likelihood of the zero-inflated negative binomial model
#'
#' Given a vector of counts, this function computes the sum of the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (ZINB) model. The ZINB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param Y the vector of counts
#' @param mu the vector of mean parameters of the negative binomial
#' @param theta the vector of dispersion parameters of the negative binomial, or
#'   a single scalar is also possible if the dispersion parameter is constant.
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
#' @param logitPi the vector of logit of the probabilities of the zero component
#' @export
#' @return the log-likelihood of the model.
#' @importFrom copula log1pexp
#' @importFrom stats dnbinom optim optimize rbinom rnbinom runif var
#' @examples
#' n <- 10
#' mu <- seq(10,50,length.out=n)
#' logitPi <- rnorm(10)
#' zeta <- rnorm(10)
#' Y <- rnbinom(n=n, size=exp(zeta), mu=mu)
#' zinb.loglik(Y, mu, exp(zeta), logitPi)
#' zinb.loglik(Y, mu, 1, logitPi)
zinb.loglik <- function(Y, mu, theta, logitPi) {
  
  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  
  # contribution of zero inflation
  lognorm <- - copula::log1pexp(logitPi)
  
  # log-likelihood
  sum(logPnb[Y>0]) + sum(logPnb[Y==0] + copula::log1pexp(logitPi[Y==0] -
                                                           logPnb[Y==0])) + sum(lognorm)
}




#' Log-likelihood of the zero-inflated negative binomial model, for a fixed
#' dispersion parameter
#'
#' Given a unique dispersion parameter and a set of counts, together with a
#' corresponding set of mean parameters and probabilities of zero components,
#' this function computes the sum of the log-probabilities of the counts under
#' the ZINB model. The dispersion parameter is provided to the function through
#' zeta = log(theta), where theta is sometimes called the inverse dispersion
#' parameter. The probabilities of the zero components are provided through
#' their logit, in order to better numerical stability.
#'
#' @param zeta a vector, the log of the inverse dispersion parameters of the
#'   negative binomial model
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @param logitPi a vector of logit of the probabilities of the zero component
#' @export
#' @seealso \code{\link{zinb.loglik}}.
#' @return the log-likelihood of the model.
#' @examples
#' mu <- seq(10,50,length.out=10)
#' logitPi <- rnorm(10)
#' zeta <- rnorm(10)
#' Y <- rnbinom(n=10, size=exp(zeta), mu=mu)
#' zinb.loglik.dispersion(zeta, Y, mu, logitPi)
zinb.loglik.dispersion <- function(zeta, Y, mu, logitPi) {
  zinb.loglik(Y, mu, exp(zeta), logitPi)
}

#' Derivative of the log-likelihood of the zero-inflated negative binomial model
#' with respect to the log of the inverse dispersion parameter
#'
#' @param zeta the log of the inverse dispersion parameters of the negative
#'   binomial
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @param logitPi a vector of the logit of the probability of the zero component
#' @export
#' @return the gradient of the inverse dispersion parameters.
#' @seealso \code{\link{zinb.loglik}}, \code{\link{zinb.loglik.dispersion}}.
zinb.loglik.dispersion.gradient <- function(zeta, Y, mu, logitPi) {
  theta <- exp(zeta)
  
  # Check zeros in the count vector
  Y0 <- Y <= 0
  Y1 <- Y > 0
  has0 <- !is.na(match(TRUE,Y0))
  has1 <- !is.na(match(TRUE,Y1))
  
  grad <- 0
  if (has1) {
    grad <- grad + sum( theta * (digamma(Y[Y1] + theta) - digamma(theta) +
                                   zeta - log(mu[Y1] + theta) + 1 -
                                   (Y[Y1] + theta)/(mu[Y1] + theta) ) )
  }
  
  if (has0) {
    logPnb <- suppressWarnings(dnbinom(0, size = theta, mu = mu[Y0],
                                       log = TRUE))
    grad <- grad + sum( theta * (zeta - log(mu[Y0] + theta) + 1 -
                                   theta/(mu[Y0] + theta)) / (1+exp(logitPi[Y0] - logPnb)))
    # *exp(- copula::log1pexp( -logPnb + logitPi[Y0])))
    
  }
  
  grad
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




#'log-likelihood of the ZINB regression model
#'
#' This function computes the log-likelihood of a ZINB regression
#' model given a vector of counts.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi,) concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param A.pi matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (log(\theta), default=zero)
#' @details The regression model is parametrized as follows: \deqn{log(\mu) =
#'   A_\mu * a_\mu} \deqn{logit(1-\Pi) = A_\pi * a_\pi} \deqn{log(\theta) = C_\theta} 
#'   where \eqn{\mu, \Pi, \theta} are
#'   respectively the vector of mean parameters of the NB distribution, the
#'   vector of probabilities of the zero component, and the vector of inverse
#'   dispersion parameters.  The
#'   log-likelihood of a vector of parameters \eqn{\alpha = (a_\mu; a_\pi)}
#' @return the  log-likelihood.
#' @export
zinb.loglik.regression <- function(alpha, Y,
                                   A.mu = matrix(nrow=length(Y), ncol=0),
                                   A.pi = matrix(nrow=length(Y), ncol=0),
                                   C.theta = matrix(0, nrow=length(Y), ncol=1)) {
  
  # Parse the model
  r <- zinb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu,
                                  A.pi = A.pi)
  
  # Call the log likelihood function
  z <- zinb.loglik(Y, exp(r$logMu), exp(C.theta), r$logitPi)
  #return z
  z
}

#' Gradient of the  log-likelihood of the ZINB regression model
#'
#' This function computes the gradient of the log-likelihood of a ZINB
#' regression model given a vector of counts.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi) concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param A.pi matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (see Details, default=zero)
#' @details The regression model is described in
#'   \code{\link{zinb.loglik.regression}}.
#' @seealso \code{\link{zinb.loglik.regression}}
#' @return The gradient of the penalized log-likelihood.
#' @export
zinb.loglik.regression.gradient <- function(alpha, Y,
                                            A.mu = matrix(nrow=length(Y), ncol=0),
                                            A.pi = matrix(nrow=length(Y), ncol=0),
                                            C.theta = matrix(0, nrow=length(Y), ncol=1)) {
  
  # Parse the model
  r <- zinb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu,
                                  A.pi = A.pi)
  
  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)
  
  # Check zeros in the count matrix
  Y0 <- Y <= 0
  Y1 <- Y > 0
  has0 <- !is.na(match(TRUE,Y0))
  has1 <- !is.na(match(TRUE,Y1))
  
  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0 
  need.wres.pi <- r$dim.alpha[2] >0 
  
  # Compute some useful quantities
  muz <- 1/(1+exp(-r$logitPi))
  clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)
  
  lognorm <- -r$logitPi - copula::log1pexp(-r$logitPi)
  
  dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)
  
  # Compute the partial derivatives we need
  ## w.r.t. mu
  if (need.wres.mu) {
    wres_mu <- numeric(length = n)
    if (has1) {
      wres_mu[Y1] <- Y[Y1] - mu[Y1] *
        (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1])
    }
    if (has0) {
      
      wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 +
                            C.theta[Y0] - log(mu[Y0] + theta[Y0]) +
                            log(mu[Y0]))
    }
  }
  
  ## w.r.t. pi
  if (need.wres.pi) {
    wres_pi <- numeric(length = n)
    if (has1) {
     
      wres_pi[Y1] <-  muz[Y1]
      
    }
    if (has0) {
     
      wres_pi[Y0] <- -(1 - exp(clogdens0)) * muz[Y0] * (1-muz[Y0]) / dens0
    }
  }
  
  # Make gradient
  grad <- numeric(0)
  
  ## w.r.t. a_mu
  if (r$dim.alpha[1] >0) {
    istart <- r$start.alpha[1]
    iend <- r$start.alpha[1]+r$dim.alpha[1]-1
    grad <- c(grad , colSums(wres_mu * A.mu))
  }
  
  ## w.r.t. a_pi
  if (r$dim.alpha[2] >0) {
    istart <- r$start.alpha[2]
    iend <- r$start.alpha[2]+r$dim.alpha[2]-1
    grad <- c(grad , colSums(wres_pi * A.pi) )
  }
  
  
  
  grad
}

#####################
#########optimal functions
 ###zinb1: there are structures both in $mu$ and in $pi$
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

