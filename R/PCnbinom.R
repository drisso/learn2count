
#' This function estimates adjacency matrix of a NB models given a matrix of counts, using glm.
#'
#' @param alpha the sisnificant level tests
#' @param X the matrix of counts.
#' @param maxcard the uper bound for cardinalities of conditional sets K
#' extend=TRUE if we consider the union of tests

#' @return the adj matrix.

nb.wald <- function(X,maxcard,alpha, extend){
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  while (ncard <= maxcard) {
    adj.est <-  foreach(i = 1:p, .combine = "cbind",.packages = c("MASS")) %dopar%{
      neighbor <- which (adj[, i] == 1)
      if (length(neighbor) >= ncard){
        condset <- combn(neighbor, ncard, FUN = list)
        for (j in 1: length(neighbor)){
          condset.temp <- condset
          indcond <- FALSE
          k <- 1
          while (indcond == FALSE & k <= length(condset.temp)){
            if (neighbor[j] %in% condset.temp[[k]] == FALSE){
              # fit model with new edges c(neighbor[j]
              X_new <- scale(as.matrix(cbind(X[,c(neighbor[j], condset.temp[[k]])]),
                                 nrow=n, ncol=ncard+1))
              data <- data.frame(cbind(X[,i],X_new))
              colnames(data) <- paste("V", 1:(ncard+2), sep="")
              fmla <- as.formula(paste("V1 ~ ", paste(colnames(data)[-1], collapse= "+")))    
              fitadd <- try(glm.nb(fmla, data = data,link = "log"),silent = TRUE)
              if (class(fitadd) =="try-error"){
                fitadd <- glm(X[,i]~scale(X_new),family = negative.binomial(theta = 1))
              }
              
              ########## wald type tests
            
              if (summary(fitadd)$coefficients[2,4] > alpha){
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
