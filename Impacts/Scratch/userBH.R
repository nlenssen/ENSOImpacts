userBH <- function (p, n = length(p)){
    p <- as.numeric(p)

    if (all(nna <- !is.na(p)))
        nna <- TRUE
    nm <- names(p)
    p0 <- setNames(p, nm)

    p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)

    if (n <= 1)
        return(p0)
   
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    p0[nna] <- pmin(1, cummin(n/i * p[o]))[ro]
    
    return(p0)
}


fdrFunction <- function(tempField, FDR = 0.25, plotName = NULL){
    goodTests <- which(!is.na(c(tempField)))

    dataVec <- c(1-tempField)[goodTests]

    iRank <- rank(dataVec)
    iOrder <- order(dataVec)

    mTotalTests <- length(goodTests)


    # (i/m)Q, where i is the rank, m is the total number 
    # of tests, and Q is the false discovery rate you choose.
    sigLevel <-  iRank * FDR /mTotalTests

    sigTests <- which(dataVec < sigLevel)

    # The largest P value that has P<(i/m)Q is significant, and all of the P values smaller than
    # it are also significant, even the ones that aren't less than their
    # Benjamini-Hochberg critical value.

    maxTest <- max(dataVec[sigTests])

    finalSig <- which(dataVec < maxTest)

    sigInds <- goodTests[finalSig]

    finalVec <- rep(NA,length(tempField))
    finalVec[sigInds] <- 1
    sigMat <- matrix(finalVec, length(lon), length(lat))

    if(!is.null(plotName)){
        pdf(plotName,6,6)
        plot(dataVec, sigLevel,xlim=c(0,1), ylim=c(0,1), xlab= 'P-value', ylab= 'Significance Cutoff',
            main=paste('FDR =',FDR,'Cutoff =', round(maxTest,4)))
        points(dataVec[finalSig], sigLevel[finalSig], col = 'blue')
        points(dataVec[sigTests], sigLevel[sigTests], col='red')
        abline(0,1)
        dev.off()
    }
    return(sigMat)
}