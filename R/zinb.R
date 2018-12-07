#' Log-likelihood of the zero-inflated negative binomial model
#'
#' Given a vector of counts, this function computes the sum of the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (ZINB) model. For each count, the ZINB distribution is parametrized by three
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
#' @importFrom stats dnbinom
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
    sum(logPnb[Y>0]) + sum(logPnb[Y==0] + copula::log1pexp(logitPi[Y==0] -  logPnb[Y==0])) + sum(lognorm)
}
