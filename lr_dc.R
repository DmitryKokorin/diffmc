#!/usr/bin/Rscript

#linear regression

options <- commandArgs(trailingOnly = TRUE)

params.filename <- options[1]
params.H <- as.double(options[2])

#threshold <- 1000
threshold <- 100


data <- read.table(col.names=c("t",
                               "x",  "y",  "z",
                               "x2", "y2", "z2",
                               "x3", "y3", "z3",
                               "x4", "y4", "z4",
                               "x5", "y5", "z5",
                               "x6", "y6", "z6",
                               "n", "photons"),
                   sep="\t", file=params.filename)



lmfunc <- function(momentum, time) {

    d <- data.frame(momentum=momentum, time=time)
    lmr <- lm(0.5*momentum ~ time, data=d, subset=(threshold+1):(length(d$time) - 1))
#    nonlinear <- head(d, threshold)
#    preds <- as.data.frame(predict(lmr, newdata=nonlinear, interval="prediction"))
#    tmp1 <- cbind(nonlinear, 2.0*preds)
#    tmp <- subset(tmp1, (momentum > fit + (upr - lwr)*3.0 | momentum < fit - (upr - lwr)*3.0))

#    return <- tail(tmp, 1)$time
    #print (summary(lmr))
    return <- (coef(lmr)[[2]])
}



lmr.t <- lmfunc(data$x2, data$t)
#lmr.n <- lmfunc(data$x2, data$n)


cat(c(params.H, lmr.t, "\n"))
#print (lmr.t)
