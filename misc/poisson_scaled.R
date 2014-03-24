#
# test poisson scaled emissions

library(rqhmm)


scaled.poisson.test <- function(lambda, scales, values, states) {
    N = length(scales)

    #
    tm = matrix(0, nrow=N, ncol=N)
    for (i in 1:(N - 1)) {
        tm[i, i] = 1
        tm[i, i + 1] = 2
    }
    tm[N, N] = 1
    tm[N, 1] = 2
   
    #
    egrp = new.emission.groups(N, 2)
    egrp = add.emission.groups(egrp, group = list(1, 1:N))
    
    hmm <- new.qhmm(list(c(1, 1), NULL), # 2 emissions to set state path
       tm, 
       rep("discrete", N),
       rep(list(c("poisson_scaled", "discrete")), N),
       emission.groups = egrp)
    
    # initial state
    init = c(1, rep(0, N-1))
    set.initial.probs.qhmm(hmm, init)
    
    # set transition parameters
    for (i in 1:(N - 1))
        set.transition.params.qhmm(hmm, i, c(0.5, 0.5))
    set.transition.params.qhmm(hmm, N, c(0.5, 0.5))
    
    # set emission parameters
    for (i in 1:N) {
        # 1st emission track (actual poisson emissions)
        set.emission.params.qhmm(hmm, i, lambda, slot = 1)
        set.emission.option.qhmm(hmm, i, "scale", scales[i])

        # 2nd emission track forces state path
        pars = rep(0, N)
        pars[i] = 1
        set.emission.params.qhmm(hmm, i, pars, fixed = rep(T, N), slot = 2)
    }

    dset = rbind(values, states)
    result = em.qhmm(hmm, list(dset))
    
    # get lambda
    lambda = get.emission.params.qhmm(hmm, 1)
    
    return(list(result = result, lambda = lambda))
}


#
# 1. all sharing the same scale
scaled.poisson.test(0.1, c(1,1), rpois(200, lambda=1), c(rep(1, 100), rep(2, 100)))

#
# 2. different scales per state
scaled.poisson.test(0.1, c(1,2), c(rpois(100, lambda=1), rpois(100, lambda=2)), c(rep(1, 100), rep(2, 100)))
