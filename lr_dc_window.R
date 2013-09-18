#!/usr/bin/Rscript

#linear regression

options <- commandArgs(trailingOnly = TRUE)

params.filename <- options[1]
params.width <- as.integer(options[2])
params.start <- as.integer(options[3])


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
#    print (params.start)
    lmr <- lm(0.5*momentum ~ time, data=d, subset=(params.start):(params.start + params.width))
#    print (summary(lmr))
    return <- (coef(lmr)[[2]])
}



lmr.x2.t <- lmfunc(data$x2, data$t)
lmr.y2.t <- lmfunc(data$y2, data$t)


cat(c(lmr.x2.t, lmr.y2.t))
