# Find a single dispersion parameter for a count by 1-dimensional optimization of the likelihood
# Given a vector of count, this function computes a single dispersion parameter
# (log(theta)) the counts under a  negative binomial
#' (NB) model. The NB distribution is parametrized by two
#' parameters: the mean value and the dispersion of the negative binomial distribution
#'
#' @param Y the vector of counts
#' @param mu the vector mean  of the negative binomial
#' @param n the length of the vector to return
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
nb.OptimizeDispersion <- function( mu,Y,n) {

  g <- optimize(f=nb.loglik.dispersion, Y=Y, mu=mu,
                maximum=TRUE,interval=c(-100,100))

  zeta.op <- g$maximum

  zeta.ot <- try(optim( par=zeta.op, fn=nb.loglik.dispersion ,
                        gr=nb.loglik.dispersion.gradient,mu=mu,
                        Y=Y,control=list(fnscale=-1,trace=0),
                        method="BFGS")$par,silent = TRUE)
  if (!is(zeta.ot, "try-error")){
    zeta <- zeta.ot
  }else{
    zeta <- zeta.op
  }


  zeta <- rep((zeta),n)
  zeta
}



#' Log-likelihood of the  negative binomial model
#' Given a vector of counts, this function computes the sum of the
#' log-probabilities of the counts under a  negative binomial
#' (NB) model. The NB distribution is parametrized by two
#' parameters: the mean value and the dispersion of the negative binomial distribution
#' @param Y the vector of counts
#' @param mu the vector of mean parameters of the negative binomial
#' @param theta the vector of dispersion parameters of the negative binomial, or
#'   a single scalar is also possible if the dispersion parameter is constant.
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
#'
#' @return the log-likelihood of the model.
#' @importFrom stats dnbinom optim optimize rbinom rnbinom runif var
nb.loglik <- function(Y, mu, theta) {

  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

  sum(logPnb)

}


#' Log-likelihood of  negative binomial model, for a fixed
#' dispersion parameter
#'
#' Given a unique dispersion parameter and a set of counts, together with a
#' corresponding set of mean parameters,
#' this function computes the sum of the log-probabilities of the counts under
#' the NB model. The dispersion parameter is provided to the function through
#' zeta = log(theta), where theta is sometimes called the inverse dispersion
#' parameter.
#'
#' @param zeta a vector, the log of the inverse dispersion parameters of the
#'   negative binomial model
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @return the log-likelihood of the model.
nb.loglik.dispersion <- function(zeta, Y, mu){

  nb.loglik(Y, mu, exp(zeta))

}

#' Parse ZINB regression model
#'
#' Given the parameters of a NB regression model, this function parses the
#' model and computes the vector of log(mu), and the dimensions of
#' the different components of the vector of parameters. See
#' \code{\link{nb.loglik.regression}} for details of the NB regression model
#' and its parameters.
#'
#' @param alpha the vectors of parameters c(a.mu) concatenated
#' @param A.mu matrix of the model (default=empty)
#' @return A list with slot \code{logMu},
#' @seealso \code{\link{nb.loglik.regression}}
nb.regression.parseModel <- function(alpha, A.mu) {

  n <- nrow(A.mu)
  logMu <-0
  dim.alpha <- rep(0,1)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
  }


  return(list(logMu=logMu, dim.alpha=dim.alpha))

}



#'log-likelihood of the NB regression model
#'
#' This function computes the log-likelihood of a NB regression
#' model given a vector of counts.
#'
#' @param alpha the vectors of parameters a.mu concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (\eqn{log(\theta)}, default=zero)
#' @details The regression model is parametrized as follows: \deqn{log(\mu) =
#'   A_\mu * a_\mu}  \deqn{log(\theta) = C_\theta}
#'   where \eqn{\mu, \theta} are
#'   respectively the vector of mean parameters of the NB distribution,
#'    and the vector of inverse   dispersion parameters.  The
#'   log-likelihood of a vector of parameters \eqn{\alpha = a_\mu}
#' @return the log-likelihood.
nb.loglik.regression <- function(alpha, Y,
                                   A.mu = matrix(nrow=length(Y), ncol=0),
                                   C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu)

  # Call the log likelihood function
  z <- nb.loglik(Y, exp(r$logMu), exp(C.theta))
  #return z
  z
}

#' Gradient of the  log-likelihood of the NB regression model
#'
#' This function computes the gradient of the log-likelihood of a NB
#' regression model given a vector of counts.
#'
#' @param alpha the vectors of parameters a.mu concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (see Details, default=zero)
#' @details The regression model is described in
#'   \code{\link{nb.loglik.regression}}.
#' @seealso \code{\link{nb.loglik.regression}}
#' @return The gradient of the log-likelihood.
nb.loglik.regression.gradient <- function(alpha, Y,
                                            A.mu = matrix(nrow=length(Y), ncol=0),
                                            C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu)

  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)

  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0

  # Compute the partial derivatives we need
  ## w.r.t. mu
  if (need.wres.mu) {
    wres_mu <- numeric(length = n)
    wres_mu <- Y - mu *
      (Y + theta)/(mu + theta)
    wres_mu <- as.vector(wres_mu)
  }


  # Make gradient
  grad <- numeric(0)

  ## w.r.t. a_mu
  if (r$dim.alpha[1] >0) {
    grad <- c(grad , colSums(wres_mu * A.mu) )
  }



  grad
}

nb.loglik.dispersion.gradient <- function(zeta, Y, mu) {
  theta <- exp(zeta)

  grad <- 0
  grad <- grad + sum( theta * (digamma(Y + theta) - digamma(theta) +
                                 zeta - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta) ) )

  grad
}


nb.optim_funnoT <- function(beta_mu, Y, X_mu, zeta, n) {
  optim( fn=nb.loglik.regression,
         gr=nb.loglik.regression.gradient,
         par=beta_mu,
         Y=Y, A.mu=cbind(rep(1,n),X_mu),
         C.theta=matrix(zeta, nrow = n, ncol = 1),
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}


