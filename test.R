#!/usr/bin/Rscript

#linear regression

options <- commandArgs(trailingOnly = TRUE)

params.filename <- options[1]

threshold <- as.integer(options[2])
print(threshold)

data <- read.table(col.names=c("x2","n"),
                   sep=" ", file=params.filename)


#data <- read.table(col.names=c("t",
#                               "x",  "y",  "z",
#                               "x2", "y2", "z2",
#                               "x3", "y3", "z3",
#                               "x4", "y4", "z4",
#                               "x5", "y5", "z5",
#                               "x6", "y6", "z6",
#                               "n", "photons"),
#                   sep="\t", file=params.filename)


#print(data)

lmfunc <- function(momentum, time) {

    d <- data.frame(momentum=momentum, time=time)
#    print(d)
#    lmr <- lm(0.5*x2 ~ t, data=d, subset=(threshold+1):length(data$t))
#    print(length(d$time))
    lmr <- lm(0.5*momentum ~ time, data=d, subset=(threshold+1):(length(d$time) - 1))
    print(summary(lmr))
    nonlinear <- head(d, threshold)
#    print(nonlinear)
    preds <- as.data.frame(predict(lmr, newdata=nonlinear, interval="prediction"))
    tmp1 <- cbind(nonlinear, 2.0*preds)
#    print(tmp1)
    tmp <- subset(tmp1, (momentum > upr | momentum < lwr))
#    print(tmp)

    return <- tail(tmp, 1)$time
}


lmr.n <- lmfunc(data$x2, data$n)

print(lmr.n)
