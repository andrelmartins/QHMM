#
# Dishonest casino example from Durbin, et. al.
#
library(rqhmm)

# our encoding: die dots 1,2,3,4,5,6
#               states: Fair die 1, Loaded die 2 

data = readLines("casino.txt", 12)
rolls = do.call(c, lapply(strsplit(data, "")[seq(1,10,2)], as.numeric))
true.path = (do.call(c, strsplit(data, "")[seq(2,10,2)]) == "L") + 1

# 
pretty.print <- function(rolls, path, true.path) {
  labels = c("F", "L")
  for (i in 1:5) {
    start = (i - 1)*60 + 1
    end = start + 60 - 1
    cat("Rolls   ", rolls[start:end], "\n", sep='')
    cat("Die     ", labels[true.path[start:end]], "\n", sep='')
    cat("Viterbi ", labels[path[start:end]], "\n", sep='')
    cat("\n")
  }
}

# HMM

hmm = new.qhmm(list(1, NULL),
  rbind(c(1,2),c(1,2)),
  c("discrete", "discrete"),
  list("discrete", "discrete"))

set.initial.probs.qhmm(hmm, c(1, 0))
set.transition.params.qhmm(hmm, 1, c(0.95, 0.05))
set.transition.params.qhmm(hmm, 2, c(0.1, 0.9))
set.emission.params.qhmm(hmm, 1, rep(1/6, 6))
set.emission.params.qhmm(hmm, 2, c(rep(1/10, 5), 1/2))

# run viterbi

path = viterbi.qhmm(hmm, rolls)

pretty.print(rolls, path, true.path) # outputs figure 3.5 from Durbin's book

# posterior decoding
fw = forward.qhmm(hmm, rolls)
bk = backward.qhmm(hmm, rolls)

logPx = attributes(fw)$loglik

posterior = exp(fw + bk - logPx)

plot(1:300, posterior[1,], type='l', xlab="roll", ylab="P(fair)") # figure 3.6 (except the gray overlays)
