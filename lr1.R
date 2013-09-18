#!/usr/bin/Rscript

#linear regression
#library(ggplot2)

#options <- commandArgs(trailingOnly = TRUE)

#params.filename <- options[1]
#params.H <- as.double(options[2])
params.filename="workstation-2012.10.28/output-workstation-starkexp-1.0T-512x32768-01-e-pi_2-20121028/out.txt"

threshold <- 5000

data <- read.table(col.names=c("t",
                               "x",  "y",  "z",
                               "x2", "y2", "z2",
                               "x3", "y3", "z3",
                               "x4", "y4", "z4",
                               "x5", "y5", "z5",
                               "x6", "y6", "z6",
                               "n", "photons"),
                   sep="\t", file=params.filename)
X11()
#plot(data$x2, data$t)

lmfunc <- function(x2, t) {

    d <- data.frame(x2=x2, t=t)
    lmr <- lm(0.5*x2 ~ t, data=d, subset=(threshold+1):length(data$t))
    nonlinear <- head(d, threshold)
#    print(lmr)
#    autoplot(lmr)
#    print(summary(lmr))
#    print(anova(lmr))
    tmpdata <- nonlinear

    plot(tmpdata$t, tmpdata$x2, lwd=2, type="l")
    grid()
    preds <- as.data.frame(predict(lmr, newdata=tmpdata, interval="prediction"))
#    print(preds)
    lines(tmpdata$t, 2.0*preds$lwr, col="red", lwd=2)
    lines(tmpdata$t, 2.0*preds$upr, col="blue", lwd=2)

    tmp <- subset(cbind(nonlinear, 2.0*preds), (x2 > upr))

    return <- tail(tmp, 1)$t
}


#lmr.t <- lmfunc(data$x2, data$t)
lmr.n <- lmfunc(data$x2, data$n)
#plot(lmr.n)

#message("Press Return To Continue")
#invisible(readLines("stdin", n=1))


#cat(c(params.H, lmr.t, lmr.n, "\n"))
