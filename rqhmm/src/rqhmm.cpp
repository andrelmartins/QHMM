#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "func_entry.hpp"
#include <transitions/discrete.hpp>
#include <transitions/autocorr.hpp>
#include <emissions/poisson.hpp>
#include <emissions/discrete.hpp>
#include <emissions/geometric.hpp>
#include <emissions/direct.hpp>
#include <emissions/fixed.hpp>
#include <emissions/discrete_gamma.hpp>
#include <hmm.hpp>
#include <utils.hpp>
#include <vector>

static std::vector<FuncEntry*> __emissions;
static std::vector<FuncEntry*> __transitions;

class RQHMMData {
public:
  RQHMMData(int n_states, SEXP data_shape) {
    /* initial state probabilities */
    double lp = -log(n_states);
    
    init_log_probs = new double[n_states];
    for (int i = 0; i < n_states; ++i)
      init_log_probs[i] = lp;
    
    /* store copy of data shape for Iterator instantiation */
    SEXP emission_shape = VECTOR_ELT(data_shape, 0);
    SEXP covar_shape = VECTOR_ELT(data_shape, 1);

    emission_slots = Rf_length(emission_shape);
    emission_size = 0;
    if (emission_slots == 0)
      e_slot_dim = NULL;
    else {
      e_slot_dim = new int[emission_slots];
      for (int i = 0; i < emission_slots; ++i) {
        e_slot_dim[i] = INTEGER(emission_shape)[i];
        emission_size += e_slot_dim[i];
      }
    }
    
    covar_slots = Rf_length(covar_shape);
    covar_size = 0;
    if (covar_slots == 0)
      c_slot_dim = NULL;
    else {
      c_slot_dim = new int[covar_slots];
      for (int i = 0; i < covar_slots; ++i) {
        c_slot_dim[i] = INTEGER(covar_shape)[i];
        covar_size += c_slot_dim[i];
      }
    }
    
    this->n_states = n_states;
    this->supports_missing = false;
  }
  
  ~RQHMMData() {
    delete hmm;
    delete[] init_log_probs;
    if (e_slot_dim)
      delete[] e_slot_dim;
    if (c_slot_dim)
      delete[] c_slot_dim;
  }
  
  void get_dims(SEXP data, int & n, int & m) {
    if (data == R_NilValue) {
      n = 0;
      m = 0;
    }
      
    if (isMatrix(data)) {
      n = nrows(data);
      m = ncols(data);
    } else {
      n = 1;
      m = Rf_length(data);
    }
  }
  
  Iter * create_iterator(SEXP emissions, SEXP covars, SEXP missing) {
    /* check missing data support */
    if (!supports_missing && missing != R_NilValue)
      error("HMM instance does not support missing data!");

    /* validate emissions & covars */
    double * eptr;
    double * cptr = NULL;
    int * mptr = NULL;
    int L, N;
    get_dims(emissions, N, L);

    if (N != emission_size)
      error("emissions don't match data shape: n.rows = %d, required = %d", N, emission_size);

    eptr = REAL(emissions);
    
    if (covar_size > 0) {
      int Lc, Nc;
      get_dims(covars, Nc, Lc);
      
      if (Lc != L)
        error("covars length must match emissions length: %d != %d\n", Lc, L);
      if (Nc != covar_size)
        error("covars don't match data shape: n.rows = %d, required = %d", Nc, covar_size);
      
      cptr = REAL(covars);
    }

    if (missing != R_NilValue)
      mptr = INTEGER(missing);
    
    return new Iter(L, emission_slots, e_slot_dim, eptr,
                    covar_slots, c_slot_dim, cptr, mptr);
  }
  
  void fill_iterator_list(std::vector<Iter*> & iterators, SEXP emission_list, SEXP covar_list, SEXP missing_list) {
    int len = Rf_length(emission_list);

    for (int i = 0; i < len; ++i) {
      Iter * iter_i;
      SEXP covar_i = R_NilValue;
      SEXP missing_i = R_NilValue;

      if (covar_list != R_NilValue)
        covar_i = VECTOR_ELT(covar_list, i);
      if (missing_list != R_NilValue)
        missing_i = VECTOR_ELT(missing_list, i);

      iter_i = create_iterator(VECTOR_ELT(emission_list, i), covar_i, missing_i);
	
      iterators.push_back(iter_i);
    }

  }
  
  int * valid_covar_slots_copy(int * slots, int length) {
    if (covar_slots == 0)
      error("this hmm does not have any covariates!");
    
    int * shifted_copy = new int[length];
    
    for (int i = 0; i < length; ++i) {
      shifted_copy[i] = slots[i] - 1;
      if (shifted_copy[i] >= covar_slots) {
        delete shifted_copy;
        error("invalid covar slot[%d]: %d > %d", i + 1, slots[i], covar_slots);
      }
    }
    
    return shifted_copy;
  }

  void set_transition_covars(int state, int * idxs, int length) {
    int * shifted_copy = valid_covar_slots_copy(idxs, length);
    
    bool result = hmm->set_transition_covars(state, shifted_copy, length);
    delete shifted_copy;
    if (!result)
      error("covar slots not valid for state: %d", state + 1);
  }
  
  void set_emission_covars(int state, int slot, int * idxs, int length) {
    int * shifted_copy = valid_covar_slots_copy(idxs, length);
    
    bool result = hmm->set_emission_covars(state, slot, shifted_copy, length);
    delete shifted_copy;
    if (!result)
      error("covar slots not valid for state: %d [slot: %d]", state + 1, slot + 1);
  }
  
  int emission_size;
  int emission_slots;
  int * e_slot_dim;
  int covar_size;
  int covar_slots;
  int * c_slot_dim;
  HMM * hmm;
  double * init_log_probs;
  int n_states;
  bool supports_missing;
};

FuncEntry * get_entry(std::vector<FuncEntry*> & table, const char * name) {
  std::vector<FuncEntry*>::iterator it;
  for (it = table.begin(); it != table.end(); ++it)
    if (!strcmp((*it)->name, name)) {
      return *it;
    }
  return NULL;
}

int * create_target_vector(int * transition_row, int n_states, int & out_n_targets) {
  /* first pass: get number of targets */
  out_n_targets = 0;
  for (int i = 0; i < n_states; ++i)
    if (transition_row[i] > out_n_targets)
      out_n_targets = transition_row[i];
  
  /* second pass: create target array */
  int * targets = new int[out_n_targets];
  for (int i = 0; i < n_states; ++i)
    if (transition_row[i] > 0)
      targets[transition_row[i] - 1] = i;
  
  return targets;
}

template<typename TType>
void fill_transitions(TType * ttable, int n_states, SEXP valid_transitions, SEXP transitions) {
  for (int i = 0; i < n_states; ++i) {
    const char * tfunc_name = CHAR(STRING_ELT(transitions, i));
    FuncEntry * tfunc = get_entry(__transitions, tfunc_name);
    int n_targets = 0;
    int * targets;
    int * row = INTEGER(valid_transitions) + i * n_states;
    
    targets = create_target_vector(row, n_states, n_targets);
    
    ttable->insert(tfunc->create_transition_instance(n_states, i, n_targets, targets));
    
    delete[] targets;
  }
}

void process_transition_groups(TransitionTable * ttable, SEXP groups) {
  if (groups != R_NilValue) {
    int n_groups = Rf_length(groups);
      
    for (int i = 0; i < n_groups; ++i) {
      SEXP group_i = VECTOR_ELT(groups, i);
      int * igrp = INTEGER(group_i);
      int grp_len = Rf_length(group_i);
      
      ttable->makeGroup(igrp, grp_len);
    }
  }
  ttable->commitGroups(); // turn remaining singletons into unitary groups
}

template<typename EType>
RQHMMData * _create_hmm_transitions(SEXP data_shape, EType * emissions, SEXP valid_transitions, SEXP transitions, SEXP transition_groups, bool with_missing) {
  /* make choice about transition table */
  int n_states = Rf_length(transitions);
  bool needs_covars = false;
  RQHMMData * data = new RQHMMData(n_states, data_shape);
  data->supports_missing = with_missing;
  
  /* first pass: check if any transition function needs covars */
  for (int i = 0; i < n_states; ++i) {
    const char * tfunc_name = CHAR(STRING_ELT(transitions, i));
    FuncEntry * tfunc = get_entry(__transitions, tfunc_name);
    
    if (tfunc->needs_covars) {
      needs_covars = true;
      break;
    }
  }

  /* second pass: create instances */
  if (needs_covars) {
    NonHomogeneousTransitions * ttable = new NonHomogeneousTransitions(n_states);
    
    fill_transitions(ttable, n_states, valid_transitions, transitions);
    process_transition_groups(ttable, transition_groups); 

    data->hmm = HMM::create(ttable, emissions, data->init_log_probs);
  } else {
    HomogeneousTransitions * ttable = new HomogeneousTransitions(n_states);
    
    fill_transitions(ttable, n_states, valid_transitions, transitions);
    process_transition_groups(ttable, transition_groups);

    data->hmm = HMM::create(ttable, emissions, data->init_log_probs);
  }
  
  return data;
}

RQHMMData * _create_hmm(SEXP data_shape, SEXP valid_transitions, SEXP transitions, SEXP emissions, SEXP transition_groups, SEXP emission_groups, SEXP support_missing, SEXP enable_debug) {
  /* make choice about emission table */
  SEXP emission_shape = VECTOR_ELT(data_shape, 0);
  int n_emissions = Rf_length(emission_shape);
  int n_states = Rf_length(transitions);
  bool with_missing = support_missing != R_NilValue && LOGICAL(support_missing)[0] == TRUE;
  bool with_debug = enable_debug != R_NilValue && LOGICAL(enable_debug)[0] == TRUE;
  
  if (n_emissions == 1) {
    Emissions * etable = new Emissions(n_states);
    int dim = INTEGER(emission_shape)[0];
    
    for (int i = 0; i < n_states; ++i) {
      const char * efunc_name = CHAR(STRING_ELT(VECTOR_ELT(emissions, i), 0));
      FuncEntry * efunc = get_entry(__emissions, efunc_name);
      etable->insert(efunc->create_emission_instance(i, 0, dim, with_missing, with_debug));
    }
    
    // process emission groups
    if (emission_groups != R_NilValue) {
      int n_groups = Rf_length(emission_groups);
      
      for (int i = 0; i < n_groups; ++i) {
        SEXP group_i = VECTOR_ELT(emission_groups, i);
        int * igrp = INTEGER(group_i);
        int grp_len = Rf_length(group_i) - 1; // exclude slot indicator
        
        etable->makeGroup(igrp + 1, grp_len); // exclude slot indicator
      }
    }
    etable->commitGroups(); // turn remaining singletons into unitary groups
    
    return _create_hmm_transitions(data_shape, etable, valid_transitions, transitions, transition_groups, with_missing);
  } else {
    MultiEmissions * etable = new MultiEmissions(n_states, n_emissions);
    
    for (int i = 0; i < n_states; ++i) {
      std::vector<EmissionFunction *> funcs_i;
      SEXP emissions_i = VECTOR_ELT(emissions, i);
      
      for (int j = 0; j < n_emissions; ++j) {
        int dim = INTEGER(emission_shape)[j];
        const char * efunc_name = CHAR(STRING_ELT(emissions_i, j));
        FuncEntry * efunc = get_entry(__emissions, efunc_name);
        funcs_i.push_back(efunc->create_emission_instance(i, j, dim, with_missing));
      }
      
      etable->insert(funcs_i);
    }

    // process emission groups
    if (emission_groups != R_NilValue) {
      int n_groups = Rf_length(emission_groups);
      
      for (int i = 0; i < n_groups; ++i) {
        SEXP group_i = VECTOR_ELT(emission_groups, i);
        int * igrp = INTEGER(group_i);
        int grp_slot = *igrp;
        int grp_len = Rf_length(group_i) - 1; // exclude slot indicator
        
        etable->makeGroup(igrp + 1, grp_len, grp_slot); // exclude slot indicator
      }
    }
    etable->commitGroups(); // turn remaining singletons into unitary groups
    
    return _create_hmm_transitions(data_shape, etable, valid_transitions, transitions, transition_groups, with_missing);
  }
}

/* utility class */

class RQAux {
private:
  SEXP ptr;
  int n, m;
  int i, j;
  int * states;
  int * slots;
  bool has_slot;
  
  void update() {
    stateID = states[i] - 1;
    
    if (stateID < 0 || stateID >= data->n_states)
      error("invalid state number: %d", states[i]);
    
    if (has_slot) {
      slotID = slots[j] - 1;
      
      if (slotID < 0 || slotID >= data->emission_slots)
        error("invalid slot number: %d", slots[j]);
    }
  }
  
public:
  ~RQAux() {
    UNPROTECT(1);
  }
  
  void init(SEXP obj, SEXP states, SEXP slots) {
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    hmm = data->hmm;
    
    i = 0;
    j = 0;
    
    n = Rf_length(states);
    this->states = INTEGER(states);
    if (slots == R_NilValue) {
      m = 0;
      has_slot = false;
      this->slots = NULL;
    } else {
      m = Rf_length(slots);
      has_slot = true;
      this->slots = INTEGER(slots);
    }
    
    update();
  }
  
  bool is_finished() {
    return (i == n);
  }
  
  void next() {
    if (!is_finished()) {
      if (has_slot) {
        ++j;
        if (j == m) {
          j = 0;
          ++i;
          
          if (!is_finished())
            update();
        } else
	  update();
      } else {
        ++i;
        if (!is_finished())
          update();
      }
    }
  }
  
  // properties
  int stateID;
  int slotID;
  HMM * hmm;
  RQHMMData * data;
};


static SEXP convert_block_vector(std::vector<block_t> * blocks) {
  SEXP result;
  std::vector<block_t>::iterator it;
  int * ptr;
  
  if (blocks->size() == 0)
    return R_NilValue;
  
  PROTECT(result = allocMatrix(INTSXP, 2, blocks->size()));
  ptr = INTEGER(result);
  
  for (it = blocks->begin(); it != blocks->end(); ++it) {
    *ptr = (*it).start + 1; // convert to one-based
    ++ptr;
    *ptr = (*it).end + 1; // convert to one-based
    ++ptr;
  }
  
  UNPROTECT(1);
  return result;
}


extern "C" {
#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
# define attr_default __attribute__ ((visibility ("default")))
#else
# define attr_hidden
# define attr_default
#endif
#define RREGDEF(name)  R_RegisterCCallable("rqhmm", #name, (DL_FUNC) name)

  void attr_hidden register_emission(FuncEntry * function) {
    // TODO: add check to avoid duplicates
    __emissions.push_back(function);
  }
  
  void attr_hidden register_transition(FuncEntry * function) {
    // TODO: add check to avoid duplicates
    __transitions.push_back(function);
  }
  
  void attr_hidden unregister_all(const char * package) {
    std::vector<FuncEntry*>::iterator it;
    
    // drop emissions from package
    it = __emissions.begin();
    while (it != __emissions.end()) {
      FuncEntry * entry = *it;
      if (!strcmp(entry->package, package))
        it = __emissions.erase(it);
      else
        ++it;
    }
    
    // drop transitions from package
    it = __transitions.begin();
    while (it != __transitions.end()) {
      FuncEntry * entry = *it;
      if (!strcmp(entry->package, package))
        it = __transitions.erase(it);
      else
        ++it;
    }
  }
  
  void attr_hidden empty_tables() {
    std::vector<FuncEntry*>::iterator it;
    
    for (it = __emissions.begin(); it != __emissions.end(); it = __emissions.erase(it));
    for (it = __transitions.begin(); it != __transitions.end(); it = __transitions.erase(it));
  }
  
  SEXP attr_hidden fentry_exists(SEXP name, std::vector<FuncEntry*> & table) {
    SEXP result;
    const char * name_ptr;
    
    PROTECT(name = AS_CHARACTER(name));
    PROTECT(result = NEW_LOGICAL(1));
    LOGICAL(result)[0] = FALSE;
    
    name_ptr = CHAR(STRING_ELT(name, 0));
    
    std::vector<FuncEntry*>::iterator it;
    for (it = table.begin(); it != table.end(); ++it)
      if (!strcmp((*it)->name, name_ptr)) {
        LOGICAL(result)[0] = TRUE;
        break;
      }
    
    UNPROTECT(2);
    
    return result;

  }
  
  SEXP rqhmm_transition_exists(SEXP name) {
    return fentry_exists(name, __transitions);
  }
  
  SEXP rqhmm_emission_exists(SEXP name) {
    return fentry_exists(name, __emissions);
  }
  
  void rqhmm_finalizer(SEXP ptr) {
    RQHMMData * data;
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    if (!data) return;
    delete data;
    R_ClearExternalPtr(ptr);
  }
  
  // TODO: add support for missing data
  SEXP rqhmm_create_hmm(SEXP data_shape, SEXP valid_transitions, SEXP transitions, SEXP emissions, SEXP emission_groups, SEXP transition_groups, SEXP support_missing, SEXP enable_debug) {
    RQHMMData * data = _create_hmm(data_shape, valid_transitions, transitions, emissions, transition_groups, emission_groups, support_missing, enable_debug);
    SEXP ans;
    SEXP ptr;
    SEXP n_states;
    SEXP n_slots;
    
    PROTECT(ans = NEW_LIST(3));
    ptr = R_MakeExternalPtr(data, install("RQHMM_struct"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, rqhmm_finalizer, (Rboolean) TRUE);
    setAttrib(ans, install("handle_ptr"), ptr);

    PROTECT(n_states = NEW_INTEGER(1));
    INTEGER(n_states)[0] = Rf_length(transitions);

    PROTECT(n_slots = NEW_INTEGER(1));
    INTEGER(n_slots)[0] = data->emission_slots;

    SET_VECTOR_ELT(ans, 0, n_states);
    SET_VECTOR_ELT(ans, 1, n_slots);
    SET_VECTOR_ELT(ans, 2, Rf_duplicate(valid_transitions));
    
    UNPROTECT(4);
    
    return ans;
  }
  
  SEXP rqhmm_forward(SEXP rqhmm, SEXP emissions, SEXP covars, SEXP missing) {
    SEXP result;
    RQHMMData * data;
    Iter * iter;
    SEXP ptr;
    SEXP loglik;
    
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    
    /* create data structures */
    iter = data->create_iterator(emissions, covars, missing);
    PROTECT(result = allocMatrix(REALSXP, data->n_states, iter->length()));
    
    /* invoke forward */
    double log_lik = data->hmm->forward((*iter), REAL(result));

    /* clean up */
    delete iter;
    
    /* prepare result */
    PROTECT(loglik = NEW_NUMERIC(1));
    REAL(loglik)[0] = log_lik;
    setAttrib(result, install("loglik"), loglik);
    
    UNPROTECT(3);
    
    return result;
  }
  
  SEXP rqhmm_backward(SEXP rqhmm, SEXP emissions, SEXP covars, SEXP missing) {
    SEXP result;
    RQHMMData * data;
    Iter * iter;
    SEXP ptr;
    SEXP loglik;
    
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    
    /* create data structures */
    iter = data->create_iterator(emissions, covars, missing);
    PROTECT(result = allocMatrix(REALSXP, data->n_states, iter->length()));
    
    /* invoke backward */
    double log_lik = data->hmm->backward((*iter), REAL(result));
    
    /* clean up */
    delete iter;
    
    /* prepare result */
    PROTECT(loglik = NEW_NUMERIC(1));
    REAL(loglik)[0] = log_lik;
    setAttrib(result, install("loglik"), loglik);
    
    UNPROTECT(3);
    
    return result;
  }
  
  SEXP rqhmm_viterbi(SEXP rqhmm, SEXP emissions, SEXP covars, SEXP missing) {
    SEXP result;
    RQHMMData * data;
    Iter * iter;
    SEXP ptr;
    
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    
    /* create data structures */
    iter = data->create_iterator(emissions, covars, missing);
    PROTECT(result = NEW_INTEGER(iter->length()));
    
    /* invoke viterbi */
    data->hmm->viterbi((*iter), INTEGER(result));
    
    /* 0-based -> 1-based state numbers */
    int * rptr = INTEGER(result);
    for (int i = 0; i < iter->length(); ++i)
      ++rptr[i];

    /* clean up */
    delete iter;
    
    UNPROTECT(2);
    
    return result;
  }
  
  SEXP rqhmm_get_transition_params(SEXP rqhmm, SEXP state) {
    RQAux aux;
    SEXP result = R_NilValue;

    aux.init(rqhmm, state, R_NilValue);

    Params * params = aux.hmm->get_transition_params(aux.stateID);

    /* convert parameters into result vector */
    if (params != NULL) {
      SEXP fixed;
      PROTECT(result = NEW_NUMERIC(params->length()));
      PROTECT(fixed = NEW_LOGICAL(params->length()));

      for (int i = 0; i < params->length(); ++i) {
        REAL(result)[i] = (*params)[i];
        LOGICAL(fixed)[i] = (params->isFixed(i) ? TRUE : FALSE);
      }

      /* store fixed status as attribute */
      setAttrib(result, install("fixed"), fixed);

      delete params;
      UNPROTECT(2);
    }

    return result;
  }

  void apply_fixed(Params & params, SEXP fixed) {
    if (fixed != R_NilValue) {
      PROTECT(fixed = AS_LOGICAL(fixed));
      int * ptr = LOGICAL(fixed);
      int n = length(fixed);

      for (int i = 0; i < n; ++i)
        params.setFixed(i, ptr[i] == TRUE);

      UNPROTECT(1);
    }
  }

  SEXP rqhmm_set_transition_params(SEXP rqhmm, SEXP state, SEXP params, SEXP fixed) {
    RQAux aux;
    
    aux.init(rqhmm, state, R_NilValue);
    
    Params params_obj = Params(Rf_length(params), REAL(params));
    apply_fixed(params_obj, fixed);
    
    for (; !aux.is_finished(); aux.next()) {
      /* check if valid */
      bool is_valid = aux.hmm->valid_transition_params(aux.stateID, params_obj);
      if (!is_valid)
        error("param vector not valid for state: %d", aux.stateID + 1);
      else
        aux.hmm->set_transition_params(aux.stateID, params_obj);
    }

    return R_NilValue;
  }

  SEXP rqhmm_get_emission_params(SEXP rqhmm, SEXP state, SEXP slot) {
    RQAux aux;
    SEXP result = R_NilValue;
    
    aux.init(rqhmm, state, slot);

    Params * params = aux.hmm->get_emission_params(aux.stateID, aux.slotID);

    /* convert parameters into result vector */
    if (params != NULL) {
      SEXP fixed;
      PROTECT(result = NEW_NUMERIC(params->length()));
      PROTECT(fixed = NEW_LOGICAL(params->length()));

      for (int i = 0; i < params->length(); ++i) {
        REAL(result)[i] = (*params)[i];
        LOGICAL(fixed)[i] = (params->isFixed(i) ? TRUE : FALSE);
      }

      /* store fixed status as attribute */
      setAttrib(result, install("fixed"), fixed);

      delete params;
      UNPROTECT(2);
    }

    return result;
  }
  
  SEXP rqhmm_set_emission_params(SEXP rqhmm, SEXP state, SEXP slot, SEXP params, SEXP fixed) {
    RQAux aux;
    
    aux.init(rqhmm, state, slot);

    Params params_obj = Params(Rf_length(params), REAL(params));
    apply_fixed(params_obj, fixed);

    for (; !aux.is_finished(); aux.next()) {
      /* check if valid */
      bool is_valid = aux.hmm->valid_emission_params(aux.stateID, aux.slotID, params_obj);
      if (!is_valid)
        error("param vector not valid for state: %d slot: %d", aux.stateID + 1, aux.slotID + 1);
      else
        aux.hmm->set_emission_params(aux.stateID, aux.slotID, params_obj);
    }

    return R_NilValue;
  }
  
  SEXP rqhmm_set_transition_covars(SEXP rqhmm, SEXP state, SEXP covarIdxs) {
    RQAux aux;
    
    aux.init(rqhmm, state, R_NilValue);

    for (; !aux.is_finished(); aux.next())
      aux.data->set_transition_covars(aux.stateID, INTEGER(covarIdxs), Rf_length(covarIdxs));

    return R_NilValue;
  }
  
  SEXP rqhmm_set_emission_covars(SEXP rqhmm, SEXP state, SEXP slot, SEXP covarIdxs) {
    RQAux aux;
    
    aux.init(rqhmm, state, slot);

    for (; !aux.is_finished(); aux.next())
      aux.data->set_emission_covars(aux.stateID, aux.slotID, INTEGER(covarIdxs), Rf_length(covarIdxs));
    
    return R_NilValue;
  }
  
  SEXP rqhmm_get_transition_option(SEXP rqhmm, SEXP state, SEXP name) {
    RQAux aux;
    SEXP result;
    
    aux.init(rqhmm, state, R_NilValue);
    
    PROTECT(result = NEW_NUMERIC(Rf_length(name)));
    
    for (int i = 0; i < Rf_length(name); ++i) {
      double value;
      bool res = aux.hmm->get_transition_option(aux.stateID, CHAR(STRING_ELT(name, i)), &value);
      if (res)
        REAL(result)[i] = value;
      else
        REAL(result)[i] = NA_REAL;
    }
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP rqhmm_set_transition_option(SEXP rqhmm, SEXP state, SEXP name, SEXP value) {
    RQAux aux;
    SEXP result;
    
    aux.init(rqhmm, state, R_NilValue);
    
    PROTECT(result = NEW_LOGICAL(Rf_length(name)));
    
    for (int i = 0; i < Rf_length(name); ++i) {
      bool res = aux.hmm->set_transition_option(aux.stateID, CHAR(STRING_ELT(name, i)), REAL(value)[i]);
      if (res)
        LOGICAL(result)[i] = TRUE;
      else
        LOGICAL(result)[i] = FALSE;
    }
    
    UNPROTECT(1);
    
    return result;
  }

  SEXP rqhmm_get_emission_option(SEXP rqhmm, SEXP state, SEXP slot, SEXP name) {
    RQAux aux;
    SEXP result;
    
    aux.init(rqhmm, state, slot);
    
    PROTECT(result = NEW_NUMERIC(Rf_length(name)));
    
    for (int i = 0; i < Rf_length(name); ++i) {
      double value;
      bool res = aux.hmm->get_emission_option(aux.stateID, aux.slotID, CHAR(STRING_ELT(name, i)), &value);
      if (res)
        REAL(result)[i] = value;
      else
        REAL(result)[i] = NA_REAL;
    }
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP rqhmm_set_emission_option(SEXP rqhmm, SEXP state, SEXP slot, SEXP name, SEXP value) {
    RQAux aux;
    SEXP result;
    
    aux.init(rqhmm, state, slot);
    
    PROTECT(result = NEW_LOGICAL(Rf_length(name)));
    
    for (int i = 0; i < Rf_length(name); ++i) {
      bool res = aux.hmm->set_emission_option(aux.stateID, aux.slotID, CHAR(STRING_ELT(name, i)), REAL(value)[i]);
      if (res)
        LOGICAL(result)[i] = TRUE;
      else
        LOGICAL(result)[i] = FALSE;
    }
    
    UNPROTECT(1);
    
    return result;
  }
  
  SEXP rqhmm_posterior(SEXP rqhmm, SEXP emissions, SEXP covars, SEXP missing) {
    SEXP result;
    RQHMMData * data;
    Iter * iter;
    SEXP ptr;
    SEXP loglik;
    double * fw, * bk;
    
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);
    
    /* create data structures */
    iter = data->create_iterator(emissions, covars, missing);
    fw = (double*) R_alloc(data->n_states * iter->length(), sizeof(double));
    bk = (double*) R_alloc(data->n_states * iter->length(), sizeof(double));
    PROTECT(result = allocMatrix(REALSXP, iter->length(), data->n_states));
    
    /* invoke forward, backward and posterior */
    double log_lik = data->hmm->forward((*iter), fw);
    log_lik = data->hmm->backward((*iter), bk);
    data->hmm->state_posterior((*iter), fw, bk, REAL(result));
    
    /* clean up */
    delete iter;
    
    /* prepare result */
    PROTECT(loglik = NEW_NUMERIC(1));
    REAL(loglik)[0] = log_lik;
    setAttrib(result, install("loglik"), loglik);
    
    UNPROTECT(3);
    
    return result;
  }
  
  SEXP rqhmm_em(SEXP rqhmm, SEXP emissions, SEXP covars, SEXP missing, SEXP tolerance) {
    SEXP result;
    SEXP ptr;
    RQHMMData * data;
    std::vector<Iter*> iterators;
    double loglik;

    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);

    /* create iterator vector */
    data->fill_iterator_list(iterators, emissions, covars, missing);

    /* invoke */
    try {
      loglik = data->hmm->em(iterators, REAL(tolerance)[0]);
    } catch (QHMMException & e) {
      REprintf("EM::RuntimeException::%s\n", e.what());
    }

    /* clean up */
    for (unsigned int i = 0; i < iterators.size(); ++i)
      delete iterators[i];

    /* prepare result */
    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = loglik;

    UNPROTECT(2);

    return result;
  }

  SEXP rqhmm_get_initial_probs(SEXP rqhmm, SEXP probs) {
    RQHMMData * data;
    SEXP ptr;
    SEXP result;

    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);

    PROTECT(result = NEW_NUMERIC(data->n_states));

    for (int i = 0; i < data->n_states; ++i)
      REAL(result)[i] = exp(data->init_log_probs[i]);

    UNPROTECT(2);

    return result;
   }

  SEXP rqhmm_set_initial_probs(SEXP rqhmm, SEXP probs) {
    RQHMMData * data;
    SEXP ptr;
    
    /* retrieve rqhmm pointer */
    PROTECT(ptr = GET_ATTR(rqhmm, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid rqhmm object");
    data = (RQHMMData*) R_ExternalPtrAddr(ptr);

    if (Rf_length(probs) != data->n_states)
      error("probability vector length must be equal to the number of states");
    
    data->hmm->set_initial_probs(REAL(probs));
    
    UNPROTECT(1);

    return R_NilValue;
  }
  
  SEXP rqhmm_path_blocks(SEXP path, SEXP states, SEXP in_sequence) {
    SEXP result;
    PROTECT(path = AS_INTEGER(path));
    PROTECT(states = AS_INTEGER(states));
    PROTECT(in_sequence = AS_LOGICAL(in_sequence));
    std::vector<block_t> * blocks;
    std::vector<int> state_vec = std::vector<int>(INTEGER(states), INTEGER(states) + Rf_length(states));
    std::vector<int> path_vec = std::vector<int>(INTEGER(path), INTEGER(path) + Rf_length(path));
    
    if (LOGICAL(in_sequence)[0] == TRUE)
      blocks = path_blocks_seq(&path_vec, &state_vec);
    else
      blocks = path_blocks(&path_vec, &state_vec);

    result = convert_block_vector(blocks);
    
    delete blocks;
    UNPROTECT(3);

    return result;
  }
  
  SEXP rqhmm_path_blocks_ext(SEXP path, SEXP start_states, SEXP middle_states, SEXP end_states) {
    SEXP result;
    PROTECT(path = AS_INTEGER(path));
    PROTECT(start_states = AS_INTEGER(start_states));
    PROTECT(middle_states = AS_INTEGER(middle_states));
    PROTECT(end_states = AS_INTEGER(end_states));
    std::vector<block_t> * blocks;
    std::vector<int> start_vec = std::vector<int>(INTEGER(start_states), INTEGER(start_states) + Rf_length(start_states));
    std::vector<int> middle_vec = std::vector<int>(INTEGER(middle_states), INTEGER(middle_states) + Rf_length(middle_states));
    std::vector<int> end_vec = std::vector<int>(INTEGER(end_states), INTEGER(end_states) + Rf_length(end_states));
    std::vector<int> path_vec = std::vector<int>(INTEGER(path), INTEGER(path) + Rf_length(path));
    
    blocks = path_blocks(&path_vec, &start_vec, &middle_vec, &end_vec);
    result = convert_block_vector(blocks);
    delete blocks;
    
    UNPROTECT(4);
    
    return result;
  }
  
  // R Entry points
  void attr_default R_init_rqhmm(DllInfo * info) {
    Rprintf("rqhmm init called\n");
    
    // register external routines
    RREGDEF(register_emission);
    RREGDEF(register_transition);
    RREGDEF(unregister_all);
    
    // add our basic transition functions
    register_transition(new TransitionEntry<Discrete>("discrete", "rqhmm_base", false));
    register_transition(new TransitionEntry<AutoCorr>("autocorr", "rqhmm_base", false));
    register_transition(new TransitionEntry<AutoCorrCovar>("autocorr_covar", "rqhmm_base", true));

    // add our basic emission functions
    register_emission(new EmissionEntry<Poisson>("poisson", "rqhmm_base", false));
    register_emission(new EmissionEntry<PoissonCovar>("poisson_covar", "rqhmm_base", true));
    register_emission(new EmissionEntry<PoissonScaledCovar>("poisson_scaled_covar", "rqhmm_base", true));
    register_emission(new EmissionEntry<DiscreteEmissions>("discrete", "rqhmm_base", false));
    register_emission(new EmissionEntry<Geometric>("geometric", "rqhmm_base", false));
    register_emission(new EmissionEntry<DirectEmission>("direct", "rqhmm_base", false));
    register_emission(new EmissionEntry<FixedEmission>("fixed", "rqhmm_base", false));
    register_emission(new EmissionEntry<DiscreteGamma>("dgamma", "rqhmm_base", false));
  }
  
  void attr_default R_unload_rqhmm(DllInfo * info) {
    Rprintf("rqhmm unload called\n");
    
    // empty emission & transition tables
    empty_tables();
  }
  
}
