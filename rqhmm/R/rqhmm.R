#
# RQHMM: Quick Hidden Markov Model Package
#

new.qhmm <- function(data.shape, valid.transitions, transition.functions, emission.functions, transition.groups = NULL, emission.groups = NULL, support.missing = FALSE, enable.debug = FALSE) {
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
              emission.groups, transition.groups,
              as.logical(support.missing), as.logical(enable.debug))
  class(res) <- "qhmm"
  names(res) <- c("n.states", "n.emission.slots", "valid.transitions")
  res$valid.transitions = t(res$valid.transitions) # undo transpose
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
  if (!is.null(fixed))
    fixed = as.logical(fixed) # weid R behaviour. Without this (or some function that touches the param (except the stopifnot below) the variable is not initialized
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
  if (!is.null(fixed))
    fixed = as.logical(fixed) # weid R behaviour. Without this (or some function that touches the param (except the stopifnot below) the variable is not initialized
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

set.transition.covars.qhmm <- function(hmm, state, covarIdxs) {
  covarIdxs = as.integer(covarIdxs)
  state = as.integer(state)
  stopifnot(length(state) > 0)
  stopifnot(length(covarIdxs) > 0)
  stopifnot(all(state > 0))
  stopifnot(all(covarIdxs > 0))
  invisible(.Call(rqhmm_set_transition_covars, hmm, state, covarIdxs))
}

set.emission.covars.qhmm <- function(hmm, state, covarIdxs, slot = 1) {
  covarIdxs = as.integer(covarIdxs)
  state = as.integer(state)
  slot = as.integer(slot)
  stopifnot(length(state) > 0)
  stopifnot(length(slot) > 0)
  stopifnot(length(covarIdxs) > 0)
  stopifnot(all(state > 0))
  stopifnot(all(slot > 0))
  stopifnot(all(covarIdxs > 0))
  invisible(.Call(rqhmm_set_emission_covars, hmm, state, slot, covarIdxs))
}

# Maybe extend this to apply function to more than one state at a time?
#
set.transition.option.qhmm <- function(hmm, state, optName, value) {
  state = as.integer(state)
  optName = as.character(optName)
  value = as.numeric(value)

  stopifnot(length(state) == 1)
  stopifnot(all(state > 0))
  stopifnot(length(optName) > 0)
  if (length(optName) != length(value)) {
    if (length(valuee) == 1)
      value = rep(value, length(optName))
    if (length(optName) == 1)
      optName = rep(optName, length(value))
    
    stop("mismatching lengths between optName and value")
  }

  .Call(rqhmm_set_transition_option, hmm, state, optName, value)
}

get.transition.option.qhmm <- function(hmm, state, optName) {
  state = as.integer(state)
  optName = as.character(optName)

  stopifnot(length(state) == 1)
  stopifnot(all(state > 0))
  stopifnot(length(optName) > 0)

  .Call(rqhmm_get_transition_option, hmm, state, optName)
}

set.emission.option.qhmm <- function(hmm, state, optName, value, slot = 1) {
  state = as.integer(state)
  slot = as.integer(slot)
  optName = as.character(optName)
  value = as.numeric(value)
  
  stopifnot(length(state) == 1)
  stopifnot(length(slot) == 1)
  stopifnot(all(state > 0))
  stopifnot(all(slot > 0))
  stopifnot(length(optName) > 0)
  if (length(optName) != length(value)) {
    if (length(valuee) == 1)
    value = rep(value, length(optName))
    if (length(optName) == 1)
    optName = rep(optName, length(value))
    
    stop("mismatching lengths between optName and value")
  }
  
  .Call(rqhmm_set_emission_option, hmm, state, slot, optName, value)
}

get.emission.option.qhmm <- function(hmm, state, optName, slot = 1) {
  state = as.integer(state)
  slot = as.integer(slot)
  optName = as.character(optName)
  
  stopifnot(length(state) == 1)
  stopifnot(length(slot) == 1)
  stopifnot(all(state > 0))
  stopifnot(all(slot > 0))
  stopifnot(length(optName) > 0)
  
  .Call(rqhmm_get_emission_option, hmm, state, slot, optName)
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

emission.test.qhmm <- function(emission.name, emission.params, values, covars = NULL, options = NULL) {
  values = as.numeric(values) # for now all HMMs take numeric vectors as values
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

  if (!is.null(emission.params))
    set.emission.params.qhmm(hmm, 1, emission.params)

  if (!is.null(options)) {
    stopifnot(is.list(options))
    optNames = names(options)
    stopifnot(!is.null(optNames))

    for (optName in optNames)
      set.emission.option.qhmm(hmm, 1, optName, options[[optName]])
  }
  
  loglik = 0
  if (is.null(covars))
    loglik = em.qhmm(hmm, list(values))
  else
    loglik = em.qhmm(hmm, list(values), list(covars))
  
  return(list(loglik = loglik, params = get.emission.params.qhmm(hmm, 1)))
}

transition.test.qhmm <- function(transition.name, transition.params, n.states, state.sequence, covars = NULL, options = NULL) {
  stopifnot(n.states > 0)
  stopifnot(length(state.sequence) > 1)

  covar.shape = NULL
  if (!is.null(covars)) {
    if (is.vector(covars))
      covar.shape = 1
    else
      covar.shape = dim(covars)[1]
  }

  tm = NULL
  for (i in 1:n.states)
    tm = rbind(tm, 1:n.states)
  
  hmm <- new.qhmm(list(1, covar.shape), tm, rep(transition.name, n.states), as.list(rep("discrete", n.states)))

  # set emission params
  for (i in 1:n.states) {
    pars = rep(0, n.states)
    pars[i] = 1
    set.emission.params.qhmm(hmm, i, pars, fixed = rep(T, n.states))
  }

  # set initial state
  init = rep(0, n.states)
  init[state.sequence[1]] = 1
  set.initial.probs.qhmm(hmm, init)

  # set transition params & options
  for (i in 1:n.states)
    set.transition.params.qhmm(hmm, i, transition.params)
  
  if (!is.null(options)) {
    stopifnot(is.list(options))
    optNames = names(options)
    stopifnot(!is.null(optNames))

    for (optName in optNames)
      set.transition.option.qhmm(hmm, 1:n.states, optName, options[[optName]])
  }
  
  # run EM
  loglik = 0
  if (is.null(covars))
    loglik = em.qhmm(hmm, list(state.sequence))
  else
    loglik = em.qhmm(hmm, list(state.sequence), list(covars))

  # obtain transition matrix / parameter table
  tm = NULL
  for (i in 1:n.states)
    tm = rbind(tm, get.transition.params.qhmm(hmm, i))

  return(list(loglik = loglik, params = tm))
}

path.blocks.qhmm <- function(path, states, in.sequence = FALSE) {
  stopifnot(length(path) > 0)
  stopifnot(length(states) > 0)
  .Call(rqhmm_path_blocks, path, states, in.sequence)
}

path.blocks2.qhmm <- function(path, start.states, middle.states, end.states) {
  stopifnot(length(path) > 0)
  stopifnot(length(start.states) > 0)
  stopifnot(length(middle.states) > 0)
  stopifnot(length(end.states) > 0)
  .Call(rqhmm_path_blocks_ext, path, start.states, middle.states, end.states)
}

# collect/restore
#

collect.params.qhmm <- function(hmm) {
  transitions = lapply(1:(hmm$n.states), function(state)
    get.transition.params.qhmm(hmm, state))
  emissions = lapply(1:(hmm$n.states), function(state)
    lapply(1:(hmm$n.emission.slots), function(slot)
           get.emission.params.qhmm(hmm, state, slot = slot)))

  return(list(transitions = transitions, emissions = emissions))
}

restore.params.qhmm <- function(hmm, saved) {
  for (state in 1:(hmm$n.states)) {
    set.transition.params.qhmm(hmm, state, saved$transitions[[state]])

    state.params = saved$emissions[[state]]
    
    for (slot in 1:(hmm$n.emission.slots))
      set.emission.params.qhmm(hmm, state, state.params[[slot]], slot = slot)
  }
}

print.qhmm <- function(object, ...) {
  summary.qhmm(object, ...)
}

summary.qhmm <- function(object, digits = 3, nsmall = 0L, ...) {
  param.cat <- function(param.vec) {
    is.fixed = attr(param.vec, "fixed")
    for (i in 1:length(param.vec)) {
      valStr = format(param.vec[i], digits = digits, nsmall = nsmall)
      if (is.fixed[i])
        cat(" [", valStr, "]", sep='')
      else
        cat(" ", valStr, sep='')
    }
  }
  
  values = collect.params.qhmm(object)

  for (state in 1:(object$n.states)) {
    targets = which(object$valid.transitions[state,] > 0)
    cat("( state", state, ") ->", targets, "\n")
    cat("  transitions:")
    param.cat(values$transitions[[state]])
    cat("\n")
    cat("  emissions:\n")
    for (slot in 1:(object$n.emission.slots)) {
      cat("    [", slot, "]:", sep='')
      param.cat(values$emissions[[state]][[slot]])
      cat("\n")
    }
  }
}
