#
# RQHMM: Quick Hidden Markov Model Package
#

is.num.vec <- function(obj) is.vector(obj, mode = "numeric")

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
    stopifnot(is.matrix(emission.groups))
    stopifnot(all(emission.groups >= 0))
    stopifnot(dim(emission.groups) == c(n.states, n.slots))
    
    # coerce valid matrix to integer
    storage.mode(emission.groups) <- "integer"
    
    if (max(emission.groups) == 0) {
      # no groups
      emission.groups = NULL
    } else {
      # validate group properties
      for (i in 1:max(emission.groups)) {
        # 1. groups must have more than one element
        n.elems = sum(emission.groups == i)
        if (n.elems < 2)
          stop("invalid emission groups: group ", i, " only has ", n.elems, " members. Must have at least two.")
        
        # 2. all slots in a given group must have the same dimension
        slots = which(apply(emission.groups, 2, function(col) any(col == i)))
        if (length(unique(emission.slot.dims[slots])) != 1)
          stop("invalid emission groups: group ", i, " spans slots of differing dimension")
        
        # 3. all state/slot pairs must share the same emission distribution
        fnames = vector(mode="list", length = max(emission.groups))
        for (i in 1:n.states) {
          for (j in 1:n.slots) {
            grp = emission.groups[i, j]
            if (grp > 0)
              fnames[[grp]] = c(fnames[[grp]], emission.functions[[i]][j])
          }
        }
        bad = which(sapply(fnames, function(v) length(unique(v)) != 1))
        if (length(bad) > 0)
          stop("invalid emission groups: groups ", paste(bad, "", sep = " "), " each have inconsistent emission function names")
      }
      
      # convert to list form
      n.grps = max(emission.groups)
      state.vecs = vector(mode = "list", length = n.grps)
      slot.vecs = vector(mode = "list", length = n.grps)
      for (i in 1:n.states) {
        for (j in 1:n.slots) {
          grp = emission.groups[i, j]
          if (grp > 0) {
            state.vecs[[grp]] = c(state.vecs[[grp]], i)
            slot.vecs[[grp]] = c(slot.vecs[[grp]], j)
          }
        }
      }
      emission.groups = lapply(1:n.grps, function(grp) {
        # transform group slot and state IDs from R 1-based to C 0-based
        c.slots = as.integer(slot.vecs[[grp]] - 1)
        c.states = as.integer(state.vecs[[grp]] - 1)
        list(c.slots, c.states)
      })
    }
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
  else {
    if (is.matrix(vec)) {
      storage.mode(vec) <- "integer"
      return(vec)
    } else
      as.integer(vec)
  }
}

forward.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_forward, hmm, emissions, covars, null.or.integer(missing))
}

backward.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_backward, hmm, emissions, covars, null.or.integer(missing))
}

viterbi.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL) {
  .Call(rqhmm_viterbi, hmm, emissions, covars, null.or.integer(missing))
}

posterior.qhmm <- function(hmm, emissions, covars = NULL, missing = NULL, n_threads = 1) {
  .Call(rqhmm_posterior, hmm, emissions, covars, null.or.integer(missing), as.integer(n_threads))
}

em.qhmm <- function(hmm, emission.lst, covar.lst = NULL, missing.lst = NULL, tolerance = 1e-5, n_threads = 1) {
  stopifnot(is.list(emission.lst) && (is.null(covar.lst) || is.list(covar.lst))
            && (is.null(missing.lst) || is.list(missing.lst)))
  if (!is.null(covar.lst))
    stopifnot(length(emission.lst) == length(covar.lst))

  if (!is.null(missing.lst)) {
    stopifnot(length(emission.lst) == length(missing.lst))
    missing.lst = lapply(missing.lst, null.or.integer)
  }
  
  # do the actual call
  .Call(rqhmm_em, hmm, emission.lst, covar.lst, missing.lst, tolerance, as.integer(n_threads))
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
  
  result = NULL
  if (is.null(covars))
    result = em.qhmm(hmm, list(values))
  else
    result = em.qhmm(hmm, list(values), list(covars))
  
  return(list(result = result, params = get.emission.params.qhmm(hmm, 1)))
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
  result = NULL
  if (is.null(covars))
    result = em.qhmm(hmm, list(state.sequence))
  else
    result = em.qhmm(hmm, list(state.sequence), list(covars))

  # obtain transition matrix / parameter table
  tm = NULL
  for (i in 1:n.states)
    tm = rbind(tm, get.transition.params.qhmm(hmm, i))

  return(list(result = result, params = tm))
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

  return(list(transitions = transitions,
              emissions = emissions,
              initial.probs = get.initial.probs.qhmm(hmm)))
}

restore.params.qhmm <- function(hmm, saved) {
  for (state in 1:(hmm$n.states)) {
    set.transition.params.qhmm(hmm, state, saved$transitions[[state]])

    state.params = saved$emissions[[state]]
    
    for (slot in 1:(hmm$n.emission.slots))
      set.emission.params.qhmm(hmm, state, state.params[[slot]], slot = slot)
  }

  # for purposes of backward compatibility
  # some, in use, code did not same this value
  if ("initial.probs" %in% names(saved))
    set.initial.probs.qhmm(hmm, saved$initial.probs)
}

print.qhmm <- function(object, ...) {
  summary.qhmm(object, ...)
}

summary.qhmm <- function(object, digits = 3, nsmall = 0L, ...) {
  param.cat <- function(param.vec) {
    if (length(param.vec) == 0) {
      cat(" N/A")
      return()
    }
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
  init.probs = values$initial.probs

  for (state in 1:(object$n.states)) {
    targets = which(object$valid.transitions[state,] > 0)
    cat("( state", state, ") ->", targets, "\n")
    cat("  initial prob:", values$initial.probs[state], "\n")
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

##
## CGD: Build a matrix specifying which emissions are shared.
##      States (rows) X Slices (cols).  
##
## Matrix will have the following properties:
##  * 0s specify no grouping.
##  * Groups are specified as sequential integers (1:NGroups). 
##  * Emission functions MUST be set the same within groups.
##  * Slot dimensions MUST be the same.
new.emission.groups <- function(nstates, nslots) {
  return(matrix(0, nrow=nstates, ncol=nslots))
}


## States and slots are vectors of integers, representing the combinations of states and slots shared.
## Note that this function is NOT exported and availiable ONLY to functions inside of rqhmm.
add.emission.groups.slots.states <- function(emission_sharing_matrix, states, slots) {
  stopifnot(length(states) == length(slots))  ## States and slots indices must be shared.
  stopifnot(length(states) > 1) ## At least two states sharing emissions parameters.
  nextGroup <- max(emission_sharing_matrix) + 1
  
  for(i in 1:length(states)) {
    if (emission_sharing_matrix[states[i], slots[i]] > 0)
      stop("invalid operation: (state:", states[i], " slot:", slots[i], ") already belongs to a group")
    emission_sharing_matrix[states[i], slots[i]] <- nextGroup
  }
  return(emission_sharing_matrix)
}

## groups: 
# group is defined in one of two ways: 
# EITHER: an integer vector (slot number followed by two or more state numbers) 
# OR:     a list with two integer vectors (slots) and (states)
add.emission.groups <- function(emission_sharing_matrix, group=NULL, states=NULL, slots=NULL) {
  if(!is.null(states) & !is.null(slots) & is.null(group)) {
    return( add.emission.groups.slots.states(emission_sharing_matrix, states, slots) )
  }
  else if(is.null(states) & is.null(slots) & !is.null(group)) {
    if(!is.list(group)) {
      stopifnot(length(group) > 2)
      n_states <- length(group)-1
      return( add.emission.groups.slots.states(emission_sharing_matrix, states = group[2:length(group)], slots = rep(group[1], n_states)) )
    }
    else {
      stopifnot(length(group) == 2)
      slots <- group[[1]]
      states<- group[[2]]
      n_slots <- length(slots)
      n_states <- length(states)
      slots <- c(sapply(slots, function(x) {rep(x, n_states)}))
      states <- rep(states, n_slots)
      return( add.emission.groups.slots.states(emission_sharing_matrix, states= states, slots= slots)  )
    }
  }
  else {
    stop("ERROR: Either (group) OR (states and slots) must be specified.  Not both.")
  }
 
}
