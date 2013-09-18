#!/usr/bin/Rscript

#max likehood

options <- commandArgs(trailingOnly = TRUE)

params.filename <- options[1]
params.H <- as.double(options[2])
params.maxtime <- as.double(options[3])


pos <- read.table(col.names=c("x", "y", "z"), sep="\t", file=params.filename)

normal.lik1<-function(theta, y) {

    mu     <- theta[1]
    sigma2 <- theta[2]
    n      <- length(y)
    logl   <- -.5*n*log(2*pi) -.5*n*log(sigma2) - (1/(2*sigma2))*sum((y-mu)**2)

    return <- -logl
}

results.x <- optim(c(0,1), normal.lik1, y = pos$x)$par
results.y <- optim(c(0,1), normal.lik1, y = pos$y)$par
results.z <- optim(c(0,1), normal.lik1, y = pos$z)$par

results.dx <- results.x[2]/(2.*params.maxtime)
results.dy <- results.y[2]/(2.*params.maxtime)
results.dz <- results.z[2]/(2.*params.maxtime)


cat(c(params.H, params.maxtime, results.dx, results.dy, results.dz))
cat("\n")
