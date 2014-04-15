library(rqhmm)

# should fail
hmm = new.qhmm(list(1, NULL),
    rbind(c(1,2), c(1,2)),
    rep("autocorr", 2),
    rep(list("discrete"), 2))

# should work
hmm = new.qhmm(list(1, NULL),
    rbind(c(1,2), c(2,1)),
    rep("autocorr", 2),
    rep(list("discrete"), 2))
