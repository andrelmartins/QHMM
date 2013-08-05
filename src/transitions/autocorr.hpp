#ifndef AUTOCORR_HPP
#define AUTOCORR_HPP

#include <stdexcept>
#include <cmath>
#include <limits>
#include "../base_classes.hpp"


// Simple Auto-correlation transition function
//
// Parameters: alpha   :: self-transition probability
class AutoCorr : public TransitionFunction {
  public:
    // first target in `targets` must the source state
    AutoCorr(double alpha, int n_states, int n_targets, int * targets) {
      _log_probs = new double[n_states];
      
      // set all to -Inf
      for (int i = 0; i < n_states; ++i)
        _log_probs[i] = -std::numeric_limits<double>::infinity();
      
      // set actual probabilities
      _log_probs[targets[0]] = log(alpha); // self transition
      
      if (n_targets > 1) {
        double other = log(1 - alpha) - log(n_targets - 1);
        for (int i = 1; i < n_targets; ++i)
          _log_probs[targets[i]] = other;
      }
    }
    
    ~AutoCorr() {
      delete[] _log_probs;
    }
    
    virtual double log_probability(int target) const {
      return _log_probs[target];
    }
    
		virtual double log_probability(Iter const & iter, int target) const {
      return log_probability(target);
    }
    
  private:
    double * _log_probs;
};


// Covariate based Auto-correlation transition function
//
// Parameters: alpha taken from covariate slot
class AutoCorrCovar : public TransitionFunction {
  public:
    // first target in `targets` must the source state
    AutoCorrCovar(int n_states, int covar_slot, int n_targets, int * targets) : _covar_slot(covar_slot), _self(*targets), _log_base(log(n_targets - 1)) {
      _valid_states = new bool[n_states];
      
      // set all to false
      for (int i = 0; i < n_states; ++i)
        _valid_states[i] = false;
      
      // valid states
      for (int i = 0; i < n_targets; ++i)
        _valid_states[targets[i]] = true;
    }
    
    ~AutoCorrCovar() {
      delete[] _valid_states;
    }
    
    virtual double log_probability(int target) const {
      throw std::logic_error("called homogeneous version of transition log_probability on non-homogeneous class");
    }
    
		virtual double log_probability(Iter const & iter, int target) const {
      if (_valid_states[target]) {
        double alpha = iter.covar(_covar_slot);
        if (target == _self)
          return log(alpha);
        else
          return log(1.0 - alpha) - _log_base; 
      } else
        return -std::numeric_limits<double>::infinity();
    }
    
  private:
    const int _covar_slot;
    const int _self;
    const double _log_base;
    bool * _valid_states;
};

#endif