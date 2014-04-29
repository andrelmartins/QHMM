library(rqhmm)

hmm = new.qhmm(list(1, NULL),
  rbind(c(1, 2, 0), c(0, 1, 2), c(2, 0, 1)),
  rep("autocorr", 3),
  rep(list("discrete"), 3))
  
set.initial.probs.qhmm(hmm, c(1, 0, 0))
set.transition.params.qhmm(hmm, 1:3, 0.5)
set.emission.params.qhmm(hmm, 1:3, c(0.5, 0.50))

seq = as.numeric(sample(2, 10, replace=TRUE))

tmp = NULL
for (i in 1:1000) tmp = rbind(tmp, stochastic.backtrace.qhmm(hmm, seq))

sapply(1:3, function(state) colSums(tmp == state))
