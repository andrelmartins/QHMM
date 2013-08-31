#
# RQHMM: Quick Hidden Markov Model Package
#

new.qhmm <- function(data.shape, valid.transitions, transition.functions, emission.functions, transition.groups = NULL, emission.groups = NULL, support.missing = FALSE) {
  #
  # HMM structure validation
  
  # number of states
  stopifnot(is.matrix(valid.transitions))
  stopifnot(length(dim(valid.transitions)) == 2)
  stopifnot(dim(valid.transitions)[1] == dim(valid.transitions)[2])
  n.states = dim(valid.transitions)[1]
  
  stopifnot(n.states == length(transition.functions))
  stopifnot(n.states == length(emission.functions))
  
  # data shape
  stopifnot(is.list(data.shape))
  stopifnot(length(data.shape) == 2)
  stopifnot(is.vector(data.shape[[1]]))
  stopifnot(is.null(data.shape[[2]]) || is.vector(data.shape[[2]]))
  
  emission.slot.dims = as.integer(data.shape[[1]])
  covar.slot.dims = as.integer(data.shape[[2]])
  
  # it's a bit weird, but I'm allowing more than one distribution to share an emission
  stopifnot(all(emission.slot.dims >= 0))
  stopifnot(all(covar.slot.dims > 0))
  stopifnot(all(is.finite(emission.slot.dims)))
  stopifnot(all(is.finite(covar.slot.dims)))

  # emission functions match shape
  n.slots = length(emission.slot.dims)
  stopifnot(all(sapply(emission.functions, length) == n.slots))
  
  # coerce valid transition matrix to integer
  storage.mode(valid.transitions) <- "integer"
  
  # valid transition setup
  for (i in 1:n.states) {
    low.i = min(valid.transitions[i,])
    max.i = max(valid.transitions[i,])
    
    stopifnot(low.i >= 0)
    stopifnot(max.i <= n.states)
    if (max.i > 0) {
      for (j in 1:max.i)
        stopifnot(sum(valid.transitions[i,] == j) == 1)
    } else if (max.i == 0)
      warning("state ", i, " has no valid outgoing transitions")
  }
  
  # transition.groups
  if (!is.null(transition.groups)) {
    stopifnot(is.list(transition.groups))
    stopifnot(all(sapply(transition.groups, is.vector)))
    
    # no duplicates
    flat = do.call("c", transition.groups)
    if (length(unique(flat)) != length(flat))
      stop("a state can only appear once and in a single transition group")
  
    # all valid state numbers
    valid.states = is.finite(flat) & flat > 0 & flat <= n.states
    if (!all(valid.states))
      stop("invalid state numbers in transition groups: ", do.call("paste", as.list(flat[!valid.states])))
    
    # transform group state IDs from R 1-based to C 0-based
    transition.groups = lapply(transition.groups, function(grp) as.integer(grp - 1))
  }
  
  # emission groups
  if (!is.null(emission.groups)) {
    stopifnot(is.list(emission.groups))
    stopifnot(all(sapply(transition.groups, is.vector)))

    emission.groups = lapply(emission.groups, as.integer)

    # prefixed by slot number
    valid.slot = sapply(emission.groups, function(grp) length(grp) > 1 && grp[1] > 0 && grp[1] <= n.slots)
    if (!all(valid.slot))
      stop("invalid emission groups: groups must be prefixed by valid slot numbers")

    # per slot number, a state can only appear in a single group
    for (slotId in 1:n.slots) {
      groups = emission.groups[sapply(emission.groups, function(grp) grp[1] == slotId)]
      group.states = lapply(groups, function(grp) grp[2:length(grp)])
      flat = do.call("c", group.states)
      
      # no duplicates
      if (length(unique(flat)) != length(flat))
        stop("invalid slot groups for slot ", slotId, ": a state appears in more than one group")
      
      # all valid state numbers
      valid.states = is.finite(flat) & flat > 0 & flat <= n.states
      if (!all(valid.states))
        stop("in slot ", slotId, ", invalid state numbers in emission groups: ", do.call("paste", as.list(flat[!valid.states])))
    }
    
    # transform group slot and state IDs from R 1-based to C 0-based
    emission.groups = lapply(emission.groups, function(grp) as.integer(grp - 1))
  }
  
  # check if function names are valid
  for (fname in transition.functions)
    if (!.Call(rqhmm_transition_exists, fname))
      stop("unknown transition function: ", fname);
  
  for (state.emissions in emission.functions) {
    for (fname in state.emissions)
      if (!.Call(rqhmm_emission_exists, fname))
        stop("unknown emission function: ", fname)
  }
  
  # check if emission functions accept given dimensions
  
  # create instance
  res = .Call(rqhmm_create_hmm,
              list(emission.slot.dims, covar.slot.dims),
              t(valid.transitions), transition.functions, emission.functions,
              emission.groups, transition.groups, as.logical(support.missing))
  class(res) <- "rqhmm"
  return(res)
}

get.transition.params.qhmm <- function(hmm, state) {
  state = as.integer(state)
  stopifnot(length(state) == 1)
  stopifnot(state >= 1)
  .Call(rqhmm_get_transition_params, hmm, state);
}

# actually, you can pass more than one state to set multiple states
# to the same params
set.transition.params.qhmm <- function(hmm, state, params, fixed = attr(params, "fixed")) {
  params = as.numeric(params)
  state = as.integer(state)
  stopifnot(all(state > 0))
  stopifnot(length(params) > 0)
  stopifnot(is.null(fixed) || length(params) == length(fixed))
  invisible(.Call(rqhmm_set_transition_params, hmm, state, params, fixed))
}

get.emission.params.qhmm <- function(hmm, state, slot = 1) {
  state = as.integer(state)
  slot = as.integer(slot)
  stopifnot(length(state) == 1 && state >= 1)
  stopifnot(length(slot) == 1 && slot >= 1)
  .Call(rqhmm_get_emission_params, hmm, state, slot)
}

set.emission.params.qhmm <- function(hmm, state, params, slot = 1, fixed = attr(params, "fixed")) {
  params = as.numeric(params)
  state = as.integer(state)
  slot = as.integer(slot)
  stopifnot(all(state > 0))
  stopifnot(all(slot > 0))
  stopifnot(length(params) > 0)
  stopifnot(is.null(fixed) || length(params) == length(fixed))
  invisible(.Call(rqhmm_set_emission_params, hmm, state, slot, params, fixed))
}

get.initial.probs.qhmm <- function(hmm) {
  .Call(rqhmm_get_initial_probs, hmm)
}

set.initial.probs.qhmm <- function(hmm, probs) {
  probs = as.double(probs)
  stopifnot(sum(probs) == 1)
  stopifnot(all(probs >= 0 & probs <= 1))
  
  invisible(.Call(rqhmm_set_initial_probs, hmm, probs))
}

null.or.integer <- function(vec) {
  if (is.null(vec))
    vec
  else
    as.integer(vec)
}

forward.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_forward, hmm, emissions, covars, null.or.integer(missing));
}

backward.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_backward, hmm, emissions, covars, null.or.integer(missing));
}

viterbi.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_viterbi, hmm, emissions, covars, null.or.integer(missing));
}

posterior.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_posterior, hmm, emissions, covars, null.or.integer(missing));
}

em.qhmm <- function(hmm, emission.lst, covar.lst = NULL, missing.lst = NULL, tolerance = 1e-5) {
  stopifnot(is.list(emission.lst) && (is.null(covar.lst) || is.list(covar.lst))
            && (is.null(missing.lst) || is.list(missing.lst)))
  if (!is.null(covar.lst))
    stopifnot(length(emission.lst) == length(covar.lst))

  if (!is.null(missing.lst)) {
    stopifnot(length(emission.lst) == length(missing.lst))
    missing.lst = lapply(missing.lst, null.or.integer)
  }
  
  # do the actual call
  .Call(rqhmm_em, hmm, emission.lst, covar.lst, missing.lst, tolerance)
}

emission.test.qhmm <- function(emission.name, emission.params, values, covars = NULL) {
  values.shape = NULL
  if (is.vector(values)) {
    values.shape = 1
  } else {
    values.shape = dim(values)[1]
  }

  covar.shape = NULL
  if (!is.null(covars)) {
    if (is.vector(covars))
      covar.shape = 1
    else
      covar.shape = dim(covars)[1]
  }
  hmm <- new.qhmm(list(values.shape, covar.shape), as.matrix(1), "discrete", list(emission.name))
  set.transition.params.qhmm(hmm, 1, 1)
  set.emission.params.qhmm(hmm, 1, emission.params)
  
  loglik = 0
  if (is.null(covars))
    loglik = em.qhmm(hmm, list(values))
  else
    loglik = em.qhmm(hmm, list(values), list(covars))
  
  return(list(loglik = loglik, params = get.emission.params.qhmm(hmm, 1)))
}
