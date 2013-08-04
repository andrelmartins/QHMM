#
# RQHMM: Quick Hidden Markov Model Package
#

new.qhmm <- function(data.shape, valid.transitions, transition.functions, emission.functions, transition.groups = NULL, emission.groups = NULL) {
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
    for (slot in 1:n.slots) {
      groups = emission.groups[sapply(emission.groups, function(grp) grp[1] == slot)]
      group.states = lapply(groups, function(grp) grp[2:length(grp)])
      flat = do.call("c", group.states)
      
      # no duplicates
      if (length(unique(flat)) != length(flat))
        stop("invalid slot groups for slot ", slot, ": a state appears in more than one group")
      
      # all valid state numbers
      valid.states = is.finite(flat) & flat > 0 & flat <= n.states
      if (!all(valid.states))
        stop("in slot ", slot, ", invalid state numbers in emission groups: ", do.call("paste", as.list(flat[!valid.states])))
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
  
  # create instance
  
}
