
library(rqhmm)

expDif <- function(ln_x1, ln_x2) {
    return(ln_x1 + log(1 - exp(ln_x2 - ln_x1)));
}

ddgamma <- function(x, shape = 1, scale = 2, shift = 1, offset = 0) {
    x = x + offset
    xlow = x + (shift - 1)
    xhigh = x + shift
    
    txlow = xlow <= 0
    vhigh = pgamma(xhigh, shape = shape, scale = scale, lower.tail = TRUE, log.p = TRUE)
    vals = expDif(vhigh, pgamma(xlow, shape = shape, scale = scale, lower.tail = TRUE, log.p = TRUE))
    vals[txlow] = vhigh
    
    return(vals)
}

x = rgamma(1000, shape=5, scale=1)
xx = round(x)

sum(ddgamma(xx, shape=5, scale=1))

emission.test.qhmm("dgamma", c(1, 2), xx)

#
#
library(rhmm)

x = rgamma(1000, shape=5, scale=1)
xx = round(x)

res = rhmm.etest("gamma", xx, params = c(1, 5, 1))
sum(res$final.probs)

res = rhmm.etest("gamma2", xx, params = c(1, 5, 1))
sum(res$final.probs)

# equivalent to gamma2
emission.test.qhmm("dgamma", c(5, 1), xx)
# equivalent to gamma
emission.test.qhmm("dgamma", c(5, 1), xx, options = list(shift = 0.5))

