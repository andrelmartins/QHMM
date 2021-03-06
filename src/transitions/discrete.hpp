#ifndef DISCRETE_HPP
#define DISCRETE_HPP

#include <cmath>
#include <limits>
#include "../base_classes.hpp"
#include "../em_base.hpp"

class Discrete : public TransitionFunction {
public:
  Discrete(int n_states, int stateID, int n_targets, int * targets) : TransitionFunction(n_states, stateID, n_targets, targets), _pseudoCount(0.0) {

    _log_probs = new double[n_states];

    // default equi-probable
    // set all to -Inf
    for (int i = 0; i < _n_states; ++i)
      _log_probs[i] = -std::numeric_limits<double>::infinity();

    double log_prob = -log(n_targets);
    for (int i = 0; i < _n_targets; ++i)
      _log_probs[_targets[i]] = log_prob;

    _all_fixed = false;
    _is_fixed = new bool[n_targets];
    for (int i = 0; i < n_targets; ++i)
      _is_fixed[i] = false;
  }

  ~Discrete() {
    delete[] _log_probs;
    delete[] _is_fixed;
  }

  virtual bool validParams(Params const & params) const {
    double sum = 0.0;
    int fixedCount = 0;

    for (int i = 0; i < params.length(); ++i) {
      sum = sum + params[i];
      fixedCount += (params.isFixed(i) ? 1 : 0);
    }
    
    if (fixedCount == params.length() - 1 && params[0] != 1.0) { /* accept setting it to one */
      log_msg("it is not valid to have only one free parameter (zero degrees of freedom)\n");
      return false;
    }

    return params.length() == _n_targets && same_probability(sum, 1.0);
  }
  
  virtual Params * getParams() const {
    double * probs = new double[_n_targets];

    for (int i = 0; i < _n_targets; ++i)
      probs[i] = exp(_log_probs[_targets[i]]);

    Params * result = new Params(_n_targets, probs);

    for (int i = 0; i < _n_targets; ++i)
        result->setFixed(i, _is_fixed[i]);

    delete[] probs;
    return result;
  }

  virtual void setParams(Params const & params) {
    _fixedTotal = 0;
    for (int i = 0; i < _n_targets; ++i) {
      _log_probs[_targets[i]] = log(params[i]);
      if (params.isFixed(i)) {
        _is_fixed[i] = true;
        _fixedTotal += params[i];
      }
    }

    _all_fixed = params.isAllFixed();
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if(!strcmp(name, "pseudo_count")) {
      *out_value = _pseudoCount;
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "pseudo_count")) {
      if (value < 0) {
        log_msg("invalid pseudo_count: %g : shoud be >= 0\n",
                value);
        return false;
      }
      _pseudoCount = value;
      return true;
    }
    return false;
  }
  
  virtual double log_probability(int target) const {
    return _log_probs[target];
  }
    
  virtual double log_probability(Iter const & iter, int target) const {
    return log_probability(target);
  }

  virtual void updateParams(EMSequences * sequences, std::vector<TransitionFunction*> * group) {
    if (_all_fixed)
      return;

    // sufficient statistics are the per target expected counts
    double expected_counts[_n_targets];

    TransitionPosteriorIterator * piter = sequences->transition_iterator(*group);

    // initialize
    for (int i = 0; i < _n_targets; ++i)
      expected_counts[i] = _pseudoCount;

    // sum expected counts
    do {
      for (unsigned int gidx = 0; gidx < group->size(); ++gidx)
        for (int tgt_idx = 0; tgt_idx < _n_targets; ++tgt_idx)
          expected_counts[tgt_idx] += piter->posterior(gidx, tgt_idx);
    } while (piter->next());

    // estimate parameters
    double scaleFactor = 1.0 - _fixedTotal;
    double normalization = 0;
    for (int i = 0; i < _n_targets; ++i)
      if (!_is_fixed[i])
        normalization += expected_counts[i];
    for (int i = 0; i < _n_targets; ++i)
      if (!_is_fixed[i])
        _log_probs[_targets[i]] = log(expected_counts[i] / normalization * scaleFactor);
    
    // propagate to other elements in the group
    std::vector<TransitionFunction*>::iterator tf_it;
    for (tf_it = group->begin(); tf_it != group->end(); ++tf_it) {
      Discrete * tf = (Discrete*) (*tf_it)->inner();
      
      if (tf != this) {
        for (int i = 0; i < _n_targets; ++i)
          tf->_log_probs[tf->_targets[i]] = _log_probs[_targets[i]];
      }
    }

    delete piter;
  }

private:
  double * _log_probs;
  bool _all_fixed;
  bool * _is_fixed;
  double _fixedTotal;
  double _pseudoCount;
};

#endif
