#!/usr/bin/Rscript

options <- commandArgs(trailingOnly = TRUE)

params.inputFilename <- options[1]
params.outputFilename <- options[2]

params.width   <- as.integer(options[3])

data <- read.table(col.names=c("t",
                               "x",  "y",  "z",
                               "x2", "y2", "z2",
                               "x3", "y3", "z3",
                               "x4", "y4", "z4",
                               "x5", "y5", "z5",
                               "x6", "y6", "z6",
                               "n", "photons"),
                   sep="\t", file=params.inputFilename)

rmse <- function(lmr)
{
    residuals <- resid(lmr)
    return <- sqrt(sum(residuals^2)/(length(residuals) - 2))
}

sigma <- function(momentum, time, idx)
{
    d <- data.frame(momentum=momentum, time=time)
    lmr <- lm(0.5*momentum ~ time, data=d, subset=idx:(idx + params.width))

    return <- rmse(lmr)
}


range <- 1:(length(data$n) - params.width)

sigmas <- sapply(range, function(x) sigma(data$x2, data$n, x))
write.table(cbind(data$n[range], sigmas), file=params.outputFilename, row.names=FALSE, col.names=FALSE)

