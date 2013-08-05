#ifndef FUNC_ENTRY_HPP
#define FUNC_ENTRY_HPP

#include <base_classes.hpp>

class FuncEntry {
public:
  
  virtual EmissionFunction * create_emission_instance(int dim) const = 0;
  virtual TransitionFunction * create_transition_instance(int n_states, int n_targets, int * targets) const = 0;
  
  const char * package;
  const char * name;
};


typedef void (*reg_func_t)(FuncEntry*);
typedef void (*unreg_all_t)(const char*);

#endif
