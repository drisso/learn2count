#' @rdname QPtransform
#' @param assay_from The assay with the input count data.
#' @param assay_to The assay in which to store the processed data.
#' @export
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays=list(counts=matrix(rpois(50, 5), ncol=10)))
#' suppressWarnings(se <- QPtransform(se))
#' se
setMethod(
    f = "QPtransform",
    signature = signature(x = "SummarizedExperiment"),
    definition = function(x, assay_from = "counts", assay_to = "processed", ...){
        assay(x, assay_to) <- t(QPtransform(t(assay(x, assay_from))))
        return(x)
})

#' @return if x is a matrix, a matrix with the processed data; if x
#'   is a SummarizedExperiment, a SummarizedExperiment object with the processed
#'   matrix as an additional assay.
#' @rdname QPtransform
#' @importFrom stats pnorm quantile
#' @export
#' @examples
#' mat <- matrix(rpois(50, 5), nrow=10)
#' suppressWarnings(QPtransform(mat))
setMethod(
    f = "QPtransform",
    signature = signature(x ="matrix"),
    definition = function(x) {

    fX <- match_quant(x)
    power_transf(fX)
})


match_quant <- function(x, quant=.95) {

    p <- ncol(x)
    dis <- apply(x, 1, quantile, quant)

    toskip <- which(dis==0)

    if(length(toskip)>0){
        qnum <- mean(dis[-toskip])
        todoX <- x[-toskip,]
        fX <- todoX / (dis[-toskip] %o% rep(1, p) / qnum)
    } else{
        qnum <- mean(dis)
        todoX <- x
        fX <- todoX / (dis %o% rep(1, p) / qnum)
    }

    return(fX)
}

power_transf <- function(fX) {

    alpha.ks <- seq(0.01, .5, length=100)

    KSstat_est <- foreach(i = 1:ncol(fX),
                          .combine = "cbind",
                          .packages=c("iZID")) %dopar%{

        ks.stat <- rep(NA, length(alpha.ks))

        for (j in 1:length(alpha.ks)) {

            x <- fX[,i]^alpha.ks[j]
            p.nb <- 1-mean(x)/var(x)
            r <- mean(x)*(1-p.nb)/p.nb
            ks.stat[j] <- KS.statistic(x,r,p.nb)
        }
        ks.stat
                          }

    sum.ks <- rowSums(KSstat_est)
    ind.alpha <- which(sum.ks==min(sum.ks))

    #### transform X to X^\alpha with alpha obtained from KS test
    datatest <- floor(fX^alpha.ks[ind.alpha])
    gene.skip <- which(colMeans(datatest)==0)
    if (length(gene.skip)>0) {
        datatest <- datatest[,-gene.skip]
    }
    return(datatest)
}

#' @importFrom iZID nb.zihmle
#' @importFrom stats pnbinom
KS.statistic <- function(x,r,p){
    mle_ori <- suppressWarnings(nb.zihmle(x, r, p, type = "zi", lowerbound=1e-04,
                         upperbound=10^6))
    probs_ori <- mle_ori[3] + (1 - mle_ori[3]) * stats::pnbinom(0:max(x),
                                                                size = ceiling(mle_ori[1]),
                                                                prob = mle_ori[2])
    step_ori = stats::stepfun(0:max(x), c(0, probs_ori))
    z = stats::knots(step_ori)
    dev = c(0, (stats::ecdf(x))(z) - step_ori(z))
    max(abs(dev))
}