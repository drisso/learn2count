#' Structure learning with negative binomial model using optim
#'
#' This function estimates the adjacency matrix of a NB model given a matrix of
#' counts, using the optim function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @import foreach
#' @importFrom stats pchisq
nbscale.noT <- function(X,maxcard,alpha, extend){
  p <- ncol(X)
  n <- nrow(X)
  ######estimate dispersion parameter
  iter.theta <- 2
  stop.epsilon <- .0001
  zeta <-  try(foreach(i = 1:p, .combine = "cbind",
                       .export = c("nb.regression.parseModel",
                                   "nb.OptimizeDispersion",
                                   "nb.loglik.dispersion",
                                   "nb.loglik",
                                   "nb.loglik.dispersion.gradient",
                                   "nb.optim_funnoT",
                                   "nb.loglik.regression.gradient",
                                   "nb.loglik.regression")) %dopar%{
                                     iter <- 1
                                     local.lik <- rep(NA,iter.theta)
                                     zeta.i <- rep(mean(X[,i])^2/(var(X[,i])-mean(X[,i])),n)
                                     ##2. Estimate parameters of nb model with zeta.i given by the first step
                                     fitadd <- try(fitadd <- nb.optim_funnoT (beta_mu= rep(1,p), Y=X[,i], X_mu=X[,-i], zeta.i, n),silent = TRUE)
                                     if(is(fitadd, "try-error")){
                                       fit <- glm(X[,i]~X[,-i],family = "poisson")
                                       fitadd <- nb.optim_funnoT (beta_mu= fit$coefficients, Y=X[,i], X_mu=X[,-i], zeta.i, n)
                                     }
                                     ####### calculate loglikelihood at the first iteration


                                     local.lik[1] <- nb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                                           A.mu = cbind(rep(1,n),X_mu=X[, -i]),
                                                                           C.theta = zeta.i)
                                     for (iter in 2:iter.theta) {
                                       #1. Estimate zeta with initial value alpha=fitadd given by previuous iteration
                                       r <- nb.regression.parseModel (alpha=fitadd, A.mu=cbind(rep(1,n),X[,-i]))
                                       zeta.temp <- nb.OptimizeDispersion ( mu=r$logMu,Y=X[,i],n)
                                       #2. Estimate parameters of nb model with zeta given by the first step
                                       fitadd.temp <- nb.optim_funnoT (beta_mu=fitadd[1:p], Y=X[,i],X_mu=X[,-i], zeta.temp, n)
                                       local.lik[iter] <- nb.loglik.regression (alpha=fitadd.temp, Y=X[,i],
                                                                                A.mu = cbind(rep(1,n),X_mu=X[,-i]),
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
  if (is.matrix(zeta)==FALSE){
  zeta <-  foreach(i = 1:p, .combine = "cbind") %dopar%{
    data <- data.frame(cbind(X[,i],X[,-i]))
    colnames(data) <- paste("V", 1:p, sep="")
    fmla <- as.formula(paste("V1 ~ ", paste(colnames(data)[-1], collapse= "+")))
    zeta.i <- rep(glm.nb(fmla, data = data,link = "log")$theta,n)
  }
}
##### estimate adj matrix
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  while (ncard <= maxcard) {
    V <-  foreach(i = 1:p, .combine = "cbind",.export = c("nb.regression.parseModel",
                                                          "nb.OptimizeDispersion",
                                                          "nb.loglik.dispersion",
                                                          "nb.loglik",
                                                          "nb.loglik.dispersion.gradient",
                                                          "nb.optim_funnoT",
                                                          "nb.loglik.regression.gradient",
                                                          "nb.loglik.regression")) %dopar%{
    neighbor <- which (adj[, i] == 1)
    if (length(neighbor) >= ncard){
      condset <- combn(neighbor, ncard, FUN = list)
      for (j in 1: length(neighbor)){
        condset.temp <- condset
        indcond <- FALSE
        k <- 1
        while (indcond == FALSE & k <= length(condset.temp)){
           if (neighbor[j] %in% condset.temp[[k]] == FALSE){

             # initial value
             beta_mu <- c(glm(X[,i]~scale(X[,c(neighbor[j], condset.temp[[k]])])
                              ,family = "poisson")$coefficients)

             # fit model with new edges
               fitadd <- nb.optim_funnoT (beta_mu=beta_mu, Y=X[,i],
                                    X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])]), zeta[,i], n)
              # calculate loglikelihood of model with new edges
                nb.loglik.add <- nb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                          A.mu = cbind(rep(1,n),X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])])),
                                                          C.theta = zeta[,i])
                # fit model without adding new edges
                if (length(condset.temp[[k]])>0){
                    fitnoadd <- nb.optim_funnoT (beta_mu=beta_mu[-2], Y=X[,i],
                                           X_mu=scale(X[, condset.temp[[k]]]), zeta[,i], n)
                    # calculate loglikelihood of model without adding new edges
                    nb.loglik.noadd <- nb.loglik.regression (alpha=fitnoadd, Y=X[,i],
                                                                A.mu = cbind(rep(1,n),X_mu=scale(X[, condset.temp[[k]]])),
                                                                C.theta = zeta[,i])
                     }else{
                          fitnoadd <- nb.optim_funnoT (beta_mu=beta_mu[c(1,2)], Y=X[,i],
                                                 X_mu=rep(0,n), zeta[,i], n)
                          # calculate loglikelihood of model without adding new edges
                          nb.loglik.noadd <- nb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                                       A.mu = cbind(rep(1,n),rep(0,n)),
                                                                       C.theta = zeta[,i])
                        }

                      goodfit.Deviance <- 2*(nb.loglik.add -nb.loglik.noadd )

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
    if (extend == TRUE){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else
      adj <- V * t(V)

    ncard <- ncard + 1
  }
  return(adj)
}
