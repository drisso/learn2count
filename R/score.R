#' Prediction scores for estimated graph
#'
#' This function returns a set of metrics useful to evaluate the goodness of the
#' estimation of the graph on simulated data, or whenever a true graph is
#' available.
#'
#' @return A vector with the number of true positives (TP), the number of false
#'   positives (FP), the number of false negatives (FN), the positive predictive
#'   value (PPV), the sensitivity (Se), and the F1 score (F1).
#'
#' @param trueG the adjacency matrix of the true graph.
#' @param estimatedG the adjacency matrix of the estimated graph.
#' @param type the type of the true graph, could be "UG" undirected graph, or "DAG" directed acyclic graph.
#'
#' @export
#' @examples
#' set.seed(123)
#' adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow=3)
#' mat <- simdata(n=100, p=3, B=adj, mu=5, mu_noise=1)
#' res <- PCzinb(mat, method="poi", alpha=0.05)
#' prediction_scores(adj, res, type= "UG")
prediction_scores <- function(trueG ,estimatedG, type){
  if (type == "UG"){
    TP <-  sum(trueG *estimatedG)/2
    FP <-  sum((estimatedG-trueG)==1)/2
    FN <- sum((estimatedG-trueG)==-1)/2
    PPV  <- TP/(TP+FP)
    Se <- TP/(TP+FN)
    F1 <- 2*(PPV*Se)/(PPV+Se)
    return (c(TP=TP,FP=FP,FN=FN,PPV=PPV,Se=Se,F1=F1))
  }else{
    TP <-  sum(trueG *estimatedG)
    FP <-  sum((estimatedG-trueG)==1)
    FN <- sum((estimatedG-trueG)==-1)
    PPV  <- TP/(TP+FP)
    Se <- TP/(TP+FN)
    F1 <- 2*(PPV*Se)/(PPV+Se)
    return (c(TP=TP,FP=FP,FN=FN,PPV=PPV,Se=Se,F1=F1))
  }
}

