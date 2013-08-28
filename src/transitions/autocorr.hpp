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
  AutoCorr(int n_states, int n_targets, int * targets, double alpha = 0.5) : TransitionFunction(n_states, n_targets, targets) {
      _log_probs = new double[n_states];

      update_log_probs(alpha);
    }
    
    ~AutoCorr() {
      delete[] _log_probs;
    }
  
    virtual bool validParams(Params const & params) const {
      return params.length() == 1 && params[0] >= 0 && params[0] <= 1;
    }
  
    virtual Params * getParams() const {
      double alpha = exp(_log_probs[_targets[0]]);
      return new Params(1, &alpha);
    }

    virtual void setParams(Params const & params) {
      update_log_probs(params[0]);
    }
  
    virtual double log_probability(int target) const {
      return _log_probs[target];
    }
    
    virtual double log_probability(Iter const & iter, int target) const {
      return log_probability(target);
    }

  virtual void updateParams(EMSequences * sequences, std::vector<TransitionFunction*> * group) {
    // sufficient statistics are the self and total expected counts
    double expected_self_count = 0;
    double expected_total_count = 0;

    EMSequences::PosteriorTransitionIterators * siter = sequences->transition_iterators(*group);

    // sum expected counts
    do {
      PostIter & piter = siter->iter();

      piter.reset();
      do {
	for (unsigned int gidx = 0; gidx < group->size(); ++gidx) {
	  expected_self_count += piter.posterior(gidx, 0);

	  for (int tgt_idx = 1; tgt_idx < _n_targets; ++tgt_idx)
	    expected_total_count += piter.posterior(gidx, tgt_idx);
	}
      } while (piter.next());
    } while (siter->next());

    // estimate parameters
    double alpha = expected_self_count / expected_total_count;
    
    std::vector<TransitionFunction*>::iterator tf_it;
    for (tf_it = group->begin(); tf_it != group->end(); ++tf_it) {
      AutoCorr * tf = (AutoCorr*) *tf_it;

      tf->update_log_probs(alpha);
    }

    delete siter;
  }

    
  private:
    double * _log_probs;
  
    void update_log_probs(double alpha) {
      // set all to -Inf
      for (int i = 0; i < _n_states; ++i)
        _log_probs[i] = -std::numeric_limits<double>::infinity();
      
      // set actual probabilities
      _log_probs[_targets[0]] = log(alpha); // self transition
      
      if (_n_targets > 1) {
        double other = log(1 - alpha) - log(_n_targets - 1);
        for (int i = 1; i < _n_targets; ++i)
          _log_probs[_targets[i]] = other;
      }
    }
};


// Covariate based Auto-correlation transition function
//
// Parameters: alpha taken from covariate slot
class AutoCorrCovar : public TransitionFunction {
  public:
    // first target in `targets` must the source state
  AutoCorrCovar(int n_states, int n_targets, int * targets, int covar_slot = 0) : TransitionFunction(n_states, n_targets, targets), _covar_slot(covar_slot), _log_base(log(n_targets - 1)) {
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
        if (target == _stateID)
          return log(alpha);
        else
          return log(1.0 - alpha) - _log_base; 
      } else
        return -std::numeric_limits<double>::infinity();
    }
    
  private:
    const int _covar_slot;
    const double _log_base;
    bool * _valid_states;
};

#endif
