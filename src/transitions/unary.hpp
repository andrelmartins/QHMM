#ifndef UNARY_HPP
#define UNARY_HPP

#include "../base_classes.hpp"
#include "../em_base.hpp"

class Unary : public TransitionFunction {
public:
  Unary(int n_states, int stateID, int n_targets, int * targets) : TransitionFunction(n_states, stateID, n_targets, targets) {
    if (n_targets != 1)
      throw std::invalid_argument("Unary: requires one and only one transition target");
  }
  
  virtual double log_probability(int target) const {
    return 0;
  }
  
  virtual double log_probability(Iter const & iter, int target) const {
    return 0;
  }
  
};

#endif
