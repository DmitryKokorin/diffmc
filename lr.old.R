#!/usr/bin/Rscript

#linear regression

options <- commandArgs(trailingOnly = TRUE)

params.filename <- options[1]
params.H <- as.double(options[2])

threshold <- 2000

data <- read.table(col.names=c("t",
                               "x",  "y",  "z",
                               "x2", "y2", "z2",
                               "x3", "y3", "z3",
                               "x4", "y4", "z4",
                               "x5", "y5", "z5",
                               "x6", "y6", "z6",
                               "n", "photons"),
                   sep="\t", file=params.filename)


lmfunc <- function(x2, t) {

    d <- data.frame(x2=x2, t=t)
    lmr <- lm(0.5*x2 ~ t, data=d, subset=(threshold+1):length(data$t))
    nonlinear <- head(d, threshold)
    preds <- as.data.frame(predict(lmr, newdata=nonlinear, interval="prediction"))
    tmp <- subset(cbind(nonlinear, 2.0*preds), (x2 > upr))

    return <- tail(tmp, 1)$t
}


lmr.t <- lmfunc(data$x2, data$t)
lmr.n <- lmfunc(data$x2, data$n)


cat(c(params.H, lmr.t, lmr.n, "\n"))
