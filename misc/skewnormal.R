#
# Test emission function skewnorm
#
library(rqhmm)
library(sn)

#
# NOTE: Tests reveal that the initial value of the skewness should 
#       be pushed in the correct direction (negative if dist skew is negative
#       and vice-versa, otherwise optimization tends to converge to positive
#       skew).
#

#
# load library fGarch for alternative code to fit snorm


# test log-prob
run.logprob.test <- function(N, location, scale, skewness) {
  params = c(location, scale, skewness)
  attr(params, "fixed") <- c(TRUE, TRUE, TRUE)

  data = rsn(N, location, scale, skewness)

  res = emission.test.qhmm("skew_normal", params, data)

  ref = dsn(data, location, scale, skewness)

  c(sum(log(ref)), res$result$loglik[1])
}

run.logprob.test(100, 1.5, 2, 0.5)
run.logprob.test(100, 1.5, 2, -0.5)


#
# indep opt test
#
run.opt <- function(N, location, scale, skewness, params, alt.method=F) {
  loglik <- function(x, y) {
    f = -sum(log(2) + dnorm(y, x[1], x[2], log=TRUE) + pnorm(y*x[3], x[3]*x[1], x[2], log.p=TRUE))
    if (!is.finite(f)) {
      print(x)
    }

    f
  }

  data = rsn(N, location, scale, skewness)

  if (alt.method)
    optim(params, loglik, y=data)
  else
    optim(params, loglik, y=data, method="L-BFGS-B", lower=c(-Inf, 0.01, 0))
}
