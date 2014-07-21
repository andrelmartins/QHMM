#
# quick test of normal distribution emissions
#

library(rqhmm)

#
# 1. simple test
#
emission.test.qhmm("normal", c(0, 1), rnorm(1000, mean=2, sd=sqrt(3)))

#
# 2. two state test
#

normal.test <- function(means, vars, values, states) {
  N = length(means)
  stopifnot(length(means) == length(vars))

  #
  tm = matrix(0, nrow=N, ncol=N)
  for (i in 1:(N - 1)) {
    tm[i, i] = 1
    tm[i, i + 1] = 2
  }
  tm[N, N] = 1
  tm[N, 1] = 2

  # 
  hmm <- new.qhmm(list(c(1, 1), NULL), # 2 emissions to set state path
                  tm, 
                  rep("discrete", N),
                  rep(list(c("normal", "discrete")), N))

  #
  # initial state
  init = c(1, rep(0, N-1))
  set.initial.probs.qhmm(hmm, init)

  # set transition parameters
  for (i in 1:(N - 1))
      set.transition.params.qhmm(hmm, i, c(0.5, 0.5))
  set.transition.params.qhmm(hmm, N, c(0.5, 0.5))

  # set emission parameters
  for (i in 1:N) {
    # 1st emission track (actual normal emissions)
    set.emission.params.qhmm(hmm, i, c(means[i], vars[i]), slot = 1)

    # 2nd emission track forces state path
    pars = rep(0, N)
    pars[i] = 1
    set.emission.params.qhmm(hmm, i, pars, fixed = rep(T, N), slot = 2)
  }

  dset = rbind(values, states)
  result = em.qhmm(hmm, list(dset))
    
  # get means
  mus = sapply(1:N, function(i) get.emission.params.qhmm(hmm, i)[1])
  vars = sapply(1:N, function(i) get.emission.params.qhmm(hmm, i)[2])
  
  return(list(result = result, mu = mus, var = vars))
}

#
#

normal.test(c(2, 3), c(1, 1), c(rnorm(100, mean = 2), rnorm(100, mean=3)), c(rep(1, 100), rep(2, 100)))
