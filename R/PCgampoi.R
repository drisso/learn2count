#' Structure learning with negative binomial model using the glmGamPoi package
#'
#' This function estimates the adjacency matrix of a NB model given a matrix of
#' counts, using the glm_gp function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @param overdispersion if TRUE a negative binomial is fit, if FALSE a Poisson.
#' @return the estimated adjacency matrix of the graph.
#' @importFrom glmGamPoi glm_gp
#' @importFrom MASS glm.nb negative.binomial
#' @importFrom utils combn
#' @importFrom stats as.formula glm
gp.lr <- function(X,maxcard,alpha, extend, overdispersion){
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  while (ncard <= maxcard) {
    adj.est <-  foreach(i = 1:p, .combine = "cbind",.packages = c("glmGamPoi")) %dopar%{
      neighbor <- which (adj[, i] == 1)
      if (length(neighbor) >= ncard){
        condset <- combn(neighbor, ncard, FUN = list)
        for (j in 1: length(neighbor)){
          condset.temp <- condset
          indcond <- FALSE
          k <- 1
          while (!indcond & k <= length(condset.temp)){
            if (!(neighbor[j] %in% condset.temp[[k]])){
              # fit model with new edges c(neighbor[j]
              X_new <- scale(as.matrix(cbind(X[,c(neighbor[j], condset.temp[[k]])]),
                                 nrow=n, ncol=ncard+1))
              data <- data.frame(X_new)
              colnames(data) <- paste("V", 1:(ncard+1), sep="")
              fmla <- as.formula(paste(" ~ ", paste(colnames(data), collapse= "+")))
              fitadd <- glm_gp(data = t(X[,i,drop=FALSE]),
                               design = fmla, col_data = data,
                               overdispersion = overdispersion)

              ########## QLR tests
              fit0 <- glm_gp(data = t(X[,i,drop=FALSE]),
                             design = ~1, overdispersion = overdispersion)
              dev <- fit0$deviances-fitadd$deviances
              pval <- 1-pchisq(dev, df=1)
              if (pval>alpha){
                adj[neighbor[j], i] <- 0
                indcond <- TRUE
              }
            }
            k <- k+1
          }
        }
      }
      #adj[,i]
    #})
      return(adj[,i])
    }

    if (extend == TRUE){
      adj <- adj.est + t(adj.est)
      adj[which(adj != 0)] <-1
    }else
    adj <- adj.est * t(adj.est)
    ncard <- ncard + 1
  }
  return(adj)
}
