#' Structure learning with zero-inflated negative binomial model (mean only)
#'
#' This function estimates the adjacency matrix of a ZINB model given a matrix
#' of counts, using the optim function.
#'
#' This approach assumes that the structure of the graph only depends on the
#' mean parameter, treating zero inflation as a technical noise effect. We call
#' this model `zinb0`.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @export
zinb0.noT <- function(X,maxcard,alpha, extend){
  p <- ncol(X)
  n <- nrow(X)
  ######estimate dispersion parameter
  iter.theta <- 2
  stop.epsilon <- .0001
  zeta <- try(foreach(i = 1:p, .combine = "cbind",.export = c("zinb.regression.parseModel",
                                                          "zinbOptimizeDispersion",
                                                          "zinb.loglik.dispersion",
                                                          "zinb.loglik",
                                                          "zinb.loglik.dispersion.gradient",
                                                          "optim_fun0noT",
                                                          "zinb.loglik.regression.gradient",
                                                          "zinb.loglik.regression")) %dopar%{
          iter <- 1
          local.lik <- rep(NA,iter.theta)
          zeta.i <- rep(mean(X[,i])^2/(var(X[,i])-mean(X[,i])),n)
          #2. Estimate parameters of ZINB model with zeta.i given by the first step
          fitadd <- try(fitadd <- optim_fun0noT (beta_mu= rep(1,p), gamma_pi=1, Y=X[,i],
                                                X_mu=X[,-i], zeta.i, n),silent = TRUE)
          if(is(fitadd, "try-error")){
             fit <- glm(X[,i]~X[,-i],family = "poisson")
             fitadd <- optim_fun0noT (beta_mu= fit$coefficients, gamma_pi=1, Y=X[,i],
                                    X_mu=X[,-i], zeta.i, n)
            }
          ####### Calculate loglikelihood at the first iteration with
          #    alpha=fitadd, and C.theta = zeta.i obtained from the above step

            local.lik[1] <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                   A.mu = cbind(rep(1,n),X_mu=X[, -i]),
                                                    A.pi= matrix(rep(1,n),n,1),
                                                    C.theta = zeta.i)
            for (iter in 2:iter.theta) {
                #1. Estimate zeta with initial value alpha=fitadd given by previuous iteration
                r <- zinb.regression.parseModel (alpha=fitadd, A.mu=cbind(rep(1,n),X[,-i]),
                                                 A.pi= matrix(rep(1,n),n,1))
                zeta.temp <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
               #2. Estimate parameters of ZINB model with zeta given by the first step
                fitadd.temp <- optim_fun0noT (beta_mu=fitadd[1:p],
                                            gamma_pi=fitadd[(p+1)],Y=X[,i],
                                            X_mu=X[,-i], zeta.temp, n)
                local.lik[iter] <- zinb.loglik.regression (alpha=fitadd.temp, Y=X[,i],
                                                           A.mu = cbind(rep(1,n),X_mu=X[,-i]),
                                                           A.pi= matrix(rep(1,n),n,1),
                                                           C.theta = zeta.temp)

                if (local.lik[iter] > local.lik[iter-1]){
                   fitadd <- fitadd.temp
                   zeta.i <- zeta.temp
                  }else break
                if(abs((local.lik[iter]-local.lik[iter-1]) /local.lik[iter-1])< stop.epsilon)
                   break

                iter <- iter+1
                 }

             return(zeta.i)
        }, silent = TRUE)

  if (!is.matrix(zeta)){
    zeta <-  foreach(i = 1:p, .combine = "cbind",.export = c("zinb.regression.parseModel",
                                                             "zinbOptimizeDispersion",
                                                             "zinb.loglik.dispersion",
                                                             "zinb.loglik",
                                                             "zinb.loglik.dispersion.gradient",
                                                             "optim_fun0noT",
                                                             "zinb.loglik.regression.gradient",
                                                             "zinb.loglik.regression")) %dopar%{
      r <- zinb.regression.parseModel (alpha=rep(1,2*p), A.mu=cbind(rep(1,n),scale(X[,-i])),
                                       A.pi= matrix(rep(1,n),n,1))
      zeta.i <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
    }
  }
  ############### Estimate adjacency matrix

  adj <- matrix(1,p,p)
  diag(adj) <- 0

  ncard <- 0
  while (ncard <= maxcard) {
    V <-  foreach(i = 1:p, .combine = "cbind",.export = c("zinb.regression.parseModel",
                                                          "zinbOptimizeDispersion",
                                                          "zinb.loglik.dispersion",
                                                          "zinb.loglik",
                                                          "zinb.loglik.dispersion.gradient",
                                                          "optim_fun0noT",
                                                          "zinb.loglik.regression.gradient",
                                                          "zinb.loglik.regression")) %dopar%{

    neighbor <- which (adj[, i] == 1)
    if (length(neighbor) >= ncard){
      condset <- combn(neighbor, ncard, FUN = list)
      for (j in 1: length(neighbor)){
        condset.temp <- condset
        indcond <- FALSE
        k <- 1
        while (!indcond & k <= length(condset.temp)){
           if (!(neighbor[j] %in% condset.temp[[k]])){

             # initial value
             beta_mu <- c(glm(X[,i]~scale(X[,c(neighbor[j], condset.temp[[k]])])
                              ,family = "poisson")$coefficients)
             gamma_pi <- 0.5

             # fit model with new edges
               fitadd <- optim_fun0noT (beta_mu=beta_mu, gamma_pi, Y=X[,i],
                                    X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])]), zeta[,i], n)
                # calculate loglikelihood of new model
                zinb.loglik.add <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                          A.mu = cbind(rep(1,n),X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])])),
                                                          A.pi= matrix(rep(1,n),n,1),
                                                          C.theta = zeta[,i])
                # fit model without adding new edges
                if (length(condset.temp[[k]])>0){
                    fitnoadd <- optim_fun0noT (beta_mu=beta_mu[-2], gamma_pi, Y=X[,i],
                                           X_mu=scale(X[, condset.temp[[k]]]), zeta[,i], n)
                    # calculate loglikelihood of model without adding new edges
                    zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitnoadd, Y=X[,i],
                                                                A.mu = cbind(rep(1,n),X_mu=scale(X[, condset.temp[[k]]])),
                                                                A.pi= matrix(rep(1,n),n,1),
                                                                C.theta = zeta[,i])
                     }else{
                          fitnoadd <- optim_fun0noT (beta_mu=beta_mu[c(1,2)], gamma_pi, Y=X[,i],
                                                 X_mu=rep(0,n), zeta[,i], n)
                          # calculate loglikelihood of model without adding new edges
                          zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                                       A.mu = cbind(rep(1,n),rep(0,n)),
                                                                       A.pi= matrix(rep(1,n),n,1),
                                                                       C.theta = zeta[,i])
                        }

                      goodfit.Deviance <- 2*abs(zinb.loglik.noadd -zinb.loglik.add )

                      if (1-pchisq(goodfit.Deviance,1) > alpha){
                          adj[neighbor[j], i] <- 0
                          indcond <- TRUE
                        }
                      }
                      k <- k+1
                    }
                  }
                }
              return(adj[,i])
              }
    if (extend){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else
      adj <- V * t(V)

    ncard <- ncard + 1
  }
  return(adj)
}
