#
#

library(rqhmm)

f <- function(x) {
  1 / (1 + exp(5 - 10 * x)) # midpoint at 0.5
}

# alternate
states = rep(c(rep(1, 10), rep(2, 10)), 10)

obs = vector(mode="numeric", length=length(states))

# state 1
obs[states == 1] = sample(c(0, 1), sum(states == 1), prob=c(0.2, 0.8), replace=T)

# state 2
obs[states == 2] = sample(c(0, 1), sum(states == 2), prob=c(0.6, 0.4), replace=T)

# covariate
x = vector(mode="numeric", length=length(states))
x[states == 2] = runif(sum(states == 1), 0, 0.5)
x[states == 1] = runif(sum(states == 2), 0.5, 1)

# HMM
hmm = new.qhmm(list(1, 1),
               rbind(c(1,2), c(1,2)),
               c("logistic", "discrete"),
               list("discrete", "discrete"))

set.initial.probs.qhmm(hmm, c(1, 0))
set.emission.params.qhmm(hmm, 1, c(0.2, 0.8), fixed = c(T, T))
set.emission.params.qhmm(hmm, 2, c(0.6, 0.4), fixed = c(T, T))
set.emission.option.qhmm(hmm, 1, "offset", 0)
set.emission.option.qhmm(hmm, 2, "offset", 0)
set.transition.params.qhmm(hmm, 1, c(1, -1))
set.transition.params.qhmm(hmm, 2, c(0.5, 0.5))

# fit
em.qhmm(hmm, list(obs), covar.lst = list(x))

plot.fitted <- function(hmm) {
  pars = get.transition.params.qhmm(hmm, 1)
  a = pars[1]
  b = pars[2]
  plot(1 / (1 + exp(a + b * x)))
}
