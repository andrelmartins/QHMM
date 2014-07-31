library(rqhmm)

x = rgamma(10000, shape=5, scale=1)
sum(dgamma(x, shape=1, scale=2, log=TRUE))

emission.test.qhmm("gamma", c(1, 2), x)
