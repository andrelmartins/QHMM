#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "func_entry.hpp"
#include <transitions/discrete.hpp>
#include <transitions/autocorr.hpp>
#include <emissions/poisson.hpp>
#include <hmm.hpp>
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

    emission_slots = length(emission_shape);
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
    
    covar_slots = length(covar_shape);
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
      m = length(data);
    }
  }
  
  Iter * create_iterator(SEXP emissions, SEXP covars) {
    /* validate emissions & covars */
    double * eptr;
    double * cptr = NULL;
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
    
    return new Iter(L, emission_slots, e_slot_dim, eptr,
                    covar_slots, c_slot_dim, cptr);
  }
  
  int emission_size;
  int emission_slots;
  int * e_slot_dim;
  int covar_size;
  int covar_slots;
  int * c_slot_dim;
  HMM * hmm;
  double * init_log_probs;
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
    
    ttable->insert(tfunc->create_transition_instance(n_states, n_targets, targets));
    
    delete[] targets;
  }
}

template<typename EType>
RQHMMData * _create_hmm_transitions(SEXP data_shape, EType * emissions, SEXP valid_transitions, SEXP transitions) {
  /* make choice about transition table */
  int n_states = length(transitions);
  bool needs_covars = false;
  RQHMMData * data = new RQHMMData(n_states, data_shape);
    
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
    
    data->hmm = HMM::create(ttable, emissions, data->init_log_probs);
  } else {
    HomogeneousTransitions * ttable = new HomogeneousTransitions(n_states);
    
    fill_transitions(ttable, n_states, valid_transitions, transitions);
    ttable->updateTransitions(); // default initialization ...
    
    data->hmm = HMM::create(ttable, emissions, data->init_log_probs);
  }
  
  return data;
}

RQHMMData * _create_hmm(SEXP data_shape, SEXP valid_transitions, SEXP transitions, SEXP emissions) {
  /* make choice about emission table */
  SEXP emission_shape = VECTOR_ELT(data_shape, 0);
  SEXP covar_shape = VECTOR_ELT(data_shape, 1);
  int n_emissions = length(emission_shape);
  int n_covars = length(covar_shape);
  int n_states = length(transitions);
  
  if (n_emissions == 1) {
    Emissions * etable = new Emissions(n_states);
    int dim = INTEGER(emission_shape)[0];
    
    for (int i = 0; i < n_states; ++i) {
      const char * efunc_name = CHAR(STRING_ELT(VECTOR_ELT(emissions, i), 0));
      FuncEntry * efunc = get_entry(__emissions, efunc_name);
      etable->insert(efunc->create_emission_instance(dim));
    }
    
    return _create_hmm_transitions(data_shape, etable, valid_transitions, transitions);
  } else {
    MultiEmissions * etable = new MultiEmissions(n_states, n_emissions);
    
    for (int i = 0; i < n_states; ++i) {
      std::vector<EmissionFunction *> funcs_i;
      SEXP emissions_i = VECTOR_ELT(emissions, i);
      
      for (int j = 0; j < n_emissions; ++j) {
        int dim = INTEGER(emission_shape)[j];
        const char * efunc_name = CHAR(STRING_ELT(emissions_i, j));
        FuncEntry * efunc = get_entry(__emissions, efunc_name);
        funcs_i.push_back(efunc->create_emission_instance(dim));
      }
      
      etable->insert(funcs_i);
    }
    
    return _create_hmm_transitions(data_shape, etable, valid_transitions, transitions);
  }
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
    while (it != __emissions.end()) {
      FuncEntry * entry = *it;
      if (!strcmp(entry->package, package))
        it = __emissions.erase(it);
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
  SEXP rqhmm_create_hmm(SEXP data_shape, SEXP valid_transitions, SEXP transitions, SEXP emissions) {
    RQHMMData * data = _create_hmm(data_shape, valid_transitions, transitions, emissions);
    SEXP ans;
    SEXP ptr;
    SEXP n_states;
    
    PROTECT(ans = NEW_LIST(1));
    ptr = R_MakeExternalPtr(data, install("RQHMM_struct"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, rqhmm_finalizer, (Rboolean) TRUE);
    setAttrib(ans, install("handle_ptr"), ptr);

    PROTECT(n_states = NEW_INTEGER(1));
    INTEGER(n_states)[0] = length(transitions);
    SET_VECTOR_ELT(ans, 0, n_states);
    
    UNPROTECT(3);
    
    return ans;
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

    // add our basic emission functions
    register_emission(new EmissionEntry<Poisson>("poisson", "rqhmm_base"));
  }
  
  void attr_default R_unload_rqhmm(DllInfo * info) {
    Rprintf("rqhmm unload called\n");
    
    // empty emission & transition tables
    empty_tables();
  }
  
}
