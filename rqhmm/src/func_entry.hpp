#ifndef FUNC_ENTRY_HPP
#define FUNC_ENTRY_HPP

#include <base_classes.hpp>

class FuncEntry {
public:
  FuncEntry(const char * fname, const char * pkg) : package(pkg), name(fname) {}
  virtual EmissionFunction * create_emission_instance(int dim) const = 0;
  virtual TransitionFunction * create_transition_instance(int n_states, int n_targets, int * targets) const = 0;
  
  const char * package;
  const char * name;
};

template<typename T>
class EmissionEntry : public FuncEntry {
public:
  EmissionEntry(const char * name, const char * package) : FuncEntry(name, package) {
  }
  
  TransitionFunction * create_transition_instance(int n_states, int n_targets, int * targets) const {
    return NULL;
  }
  
  EmissionFunction * create_emission_instance(int dim) const {
    // TODO: do something with 'dim', like maybe pass it along or check it's valid
    return new T();
  }
};

template<typename T>
class TransitionEntry : public FuncEntry {
public:
  TransitionEntry(const char * name, const char * package) : FuncEntry(name, package) {
  }
  
  TransitionFunction * create_transition_instance(int n_states, int n_targets, int * targets) const {
    return new T(n_states, n_targets, targets);
  }
  
  EmissionFunction * create_emission_instance(int dim) const {
    return NULL;
  }
};

typedef void (*reg_func_t)(FuncEntry*);
typedef void (*unreg_all_t)(const char*);

#endif
