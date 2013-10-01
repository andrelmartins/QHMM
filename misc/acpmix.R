library(rqhmm)

estimate.alpha <- function(seq, k) {
  N = length(seq)
  N1 = sum(seq[1:(N-1)] == k)
  N2 = sum(seq[1:(N-1)] == k & seq[2:N] == k)

  N2 / N1
}

# stress test loop

R = 1

for (r in 1:R) {


#
# sample a state path
#
p1 = 0.2
p2 = 0.6

K = 1000

lens1 = rgeom(K, 1-p1)
lens2 = rgeom(K, 1-p2)

states.aux = sapply(1:K, function(i) c(rep(1, lens1[i]), rep(2, lens2[i])))
states = do.call("c", states.aux)

#
# sample observations
#

p11 = 0.1
p12 = 0.9
p21 = 0.9
p22 = 0.1

em1 = sample(c(1,2), length(states), prob = c(p11, p12), replace=T)
em2 = sample(c(1,2), length(states), prob = c(p21, p22), replace=T)
obs = em1
obs[states == 2] = em2[states == 2]

# unif prior
prior1 = rep(0.5, length(states))
prior2 = rep(0.5, length(states))
prior = rbind(prior1, prior2)

#
# define HMM
#
gamma = 1 # only auo-correlation
max.iters = 10
tol = 1e-4

#
hmm <- new.qhmm(list(1, 2),
                rbind(c(1, 2), c(1, 2)),
                c("acpmix", "acpmix"),
                list("discrete", "discrete"))

set.initial.probs.qhmm(hmm, c(0.5, 0.5))
set.emission.params.qhmm(hmm, 1, c(p11, p12), fixed = c(T, T))
set.emission.params.qhmm(hmm, 2, c(p21, p22), fixed = c(T, T))

set.transition.params.qhmm(hmm, 1, c(0.5, gamma))
set.transition.params.qhmm(hmm, 2, c(0.5, gamma))
set.transition.option.qhmm(hmm, 1, "maxIters", max.iters)
set.transition.option.qhmm(hmm, 2, "maxIters", max.iters)
set.transition.option.qhmm(hmm, 1, "tolerance", tol)
set.transition.option.qhmm(hmm, 2, "tolerance", tol)

tmp = em.qhmm(hmm, list(obs), covar.lst = list(log(prior)))

alpha1 = get.transition.params.qhmm(hmm, 1)[1]
alpha2 = get.transition.params.qhmm(hmm, 2)[1]

cat("p1 =", p1, "p1* =", round(estimate.alpha(states, 1), 2), "alpha1 =", alpha1, "\n")
cat("p2 =", p2, "p2* =", round(estimate.alpha(states, 2), 2), "alpha2 =", alpha2, "\n")

}

