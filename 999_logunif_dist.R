x <- exp(runif(100000, log(20000), log(1000000)))
hist(x, breaks=100, xlim=c(20000, 1000000))
