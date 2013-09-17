#ifndef FUNC_ENTRY_HPP
#define FUNC_ENTRY_HPP

#include <base_classes.hpp>

class FuncEntry {
public:
  FuncEntry(const char * fname, const char * pkg, const bool req_covars) : needs_covars(req_covars), package(pkg), name(fname) {}
  virtual EmissionFunction * create_emission_instance(int stateID, int slotID, int dim, bool with_missing = false, bool with_debug = false) const = 0;
  virtual TransitionFunction * create_transition_instance(int n_states, int stateID, int n_targets, int * targets, bool with_debug = false) const = 0;
  
  const bool needs_covars;
  const char * package;
  const char * name;
};

template<typename T>
class EmissionEntry : public FuncEntry {
public:
  EmissionEntry(const char * name, const char * package, const bool req_covars) : FuncEntry(name, package, req_covars) {
  }
  
  TransitionFunction * create_transition_instance(int n_states, int stateID, int n_targets, int * targets, bool with_debug = false) const {
    return NULL;
  }
  
  EmissionFunction * create_emission_instance(int stateID, int slotID, int dim, bool with_missing = false, bool with_debug = false) const {
    // TODO: do something with 'dim', like maybe pass it along or check it's valid
    EmissionFunction * func = new T(stateID, slotID);

    if (with_missing)
      func = new MissingEmissionFunction(func);

    if (with_debug)
      func = new DebugEmissionFunction(func);

    return func;
  }
};

template<typename T>
class TransitionEntry : public FuncEntry {
public:
  TransitionEntry(const char * name, const char * package, const bool req_covars) : FuncEntry(name, package, req_covars) {
  }
  
  TransitionFunction * create_transition_instance(int n_states, int stateID, int n_targets, int * targets, bool with_debug = false) const {
    return new T(n_states, stateID, n_targets, targets);
  }
  
  EmissionFunction * create_emission_instance(int stateID, int slotID, int dim, bool with_missing = false, bool with_debug = false) const {
    return NULL;
  }
};

typedef void (*reg_func_t)(FuncEntry*);
typedef void (*unreg_all_t)(const char*);

#endif
