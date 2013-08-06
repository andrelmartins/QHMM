#ifndef DISCRETE_HPP
#define DISCRETE_HPP

#include <cmath>
#include <limits>
#include "../base_classes.hpp"

class Discrete : public TransitionFunction {
public:
  Discrete(int n_states, int n_targets, int * targets) : _n_states(n_states), _n_targets(n_targetes) {

    _log_probs = new double[n_states];
    _targets = new int[n_states];

    for (int i = 0; i < n_targets; ++i)
      _targets[i] = targets[i];

    // default equi-probable
    // set all to -Inf
    for (int i = 0; i < _n_states; ++i)
      _log_probs[i] = -std::numeric_limits<double>::infinity();

    double log_prob = -log(_n_states);
    for (int i = 0; i < _n_targets; ++i)
      _log_probs[_targets[i]] = log_prob;
  }

  ~Discrete() {
    delete[] _log_probs;
    delete[] _targets;
  }

  virtual bool validParams(Params const & params) const {
    double sum = 0.0;

    for (int i = 0; i < params.length(); ++i)
      sum = sum + params[i];

    // TODO: maybe allow some tolerance ??
    return params.length() == _n_states; && sum == 1.0;
  }
  
  virtual void setParams(Params const & params) {
    for (int i = 0; i < _n_targets; ++i)
      _log_probs[_targets[i]] = log(params[i]);
  }
  
  virtual double log_probability(int target) const {
    return _log_probs[target];
  }
    
  virtual double log_probability(Iter const & iter, int target) const {
    return log_probability(target);
  }

private:
  int _n_states;
  int _n_targets;
  int * _targets;
  double * _log_probs;
};

#endif
