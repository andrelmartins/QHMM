#ifndef LOGISTIC_HPP
#define LOGISTIC_HPP

#include <stdexcept>
#include <cmath>
#include <cstring>
#include <limits>
#include "../base_classes.hpp"
#include "../em_base.hpp"
#include "../math.hpp"
#include "../logsum.hpp"

class Logistic : public TransitionFunction {
public:
  
  Logistic(int n_states, int stateID, int n_targets, int * targets, int covar_slot = 0, int max_iters = 10, double tolerance = 1e-8) : TransitionFunction(n_states, stateID, n_targets, targets), _covar_slot(covar_slot), _max_iters(100), _tolerance(tolerance) {
    
    _logsum = LogSum::create(n_targets);
    
    // initial betas == 1
    int n_betas = 2*(_n_targets - 1);
    _betas = new double[n_betas];
    for (int i = 0; i < n_betas; ++i)
      _betas[i] = 1.0;
    
    // create reverse target map
    _rev_target_map = new int[n_states];
    // set all to -1
    for (int i = 0; i < n_states; ++i)
      _rev_target_map[i] = -1;
    // set targets
    for (int i = 0; i < n_targets; ++i)
      _rev_target_map[targets[i]] = i;
  }
  
  ~Logistic() {
    delete[] _betas;
    delete[] _rev_target_map;
    delete _logsum;
  }
  
  virtual bool validParams(Params const & params) const {
    return params.length() == 2*(_n_targets - 1);
  }
  
  virtual Params * getParams() const {
    Params * result = new Params(2*(_n_targets-1), _betas);
    if (_is_fixed)
      for (int i = 0; i < result->length(); ++i)
        result->setFixed(i, true);
    return result;
  }
  
  virtual void setParams(Params const & params) {
    _is_fixed = params.isAllFixed();
    
    for (int i = 0; i < params.length(); ++i)
      _betas[i] = params[i];
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "maxIters")) {
      *out_value = (double) _max_iters;
      return true;
    }
    if (!strcmp(name, "tolerance")) {
      *out_value = _tolerance;
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "maxIters")) {
      int iters = (int) value;
      if (iters <= 0) {
        log_msg("invalid maxIters: %d : should be > 0\n", iters);
        return false;
      }
      _max_iters = iters;
      return true;
    }
    if (!strcmp(name, "tolerance")) {
      if (value <= 0) {
        log_msg("invalid tolerance: %g : should be > 0\n", value);
        return false;
      }
      _tolerance = value;
      return true;
    }
    return false;
  }
  
  virtual double log_probability(int target) const {
    throw std::logic_error("called homogeneous version of transition log_probability on non-homogeneous class");
  }
  
  virtual double log_probability(Iter const & iter, int target) const {
    if (_rev_target_map[target] == -1)
      return -std::numeric_limits<double>::infinity();
    
    double x = iter.covar(_covar_slot);
    double log_num;
    
    _logsum->clear();
    _logsum->store(0);
    for (int i = 0; i < _n_targets - 1; ++i)
      _logsum->store(_betas[i*2] + _betas[i*2 + 1] * x);
    
    if (target == _targets[0])
      log_num = 0;
    else {
      int idx = _rev_target_map[target] - 1;
      log_num = _betas[idx*2] + _betas[idx*2 + 1] * x;
    }

    return log_num - _logsum->compute();
  }
  
  virtual void updateParams(EMSequences * sequences, std::vector<TransitionFunction*> * group) {
    if (_is_fixed)
      return;
    
    TransitionPosteriorIterator * piter = sequences->transition_iterator(*group);
    
    // optimize parameters
    int fail = 0;
    int n_betas = 2*(_n_targets - 1);
    struct opt_data udata;
    
    double * betas = new double[n_betas];
    for (int i = 0; i < n_betas; ++i)
      betas[i] = _betas[i];
    
    udata.group = group;
    udata.piter = piter;
    udata.logsum = _logsum;
    udata.covar_slot = _covar_slot;
    udata.n_targets = _n_targets;
    
    QHMM_fminimizer(optfunc, n_betas, betas, &udata, _max_iters, _tolerance, &fail);
    
    // TODO: do something if fail == 1
    
    // propagate parameters
    std::vector<TransitionFunction*>::iterator tf_it;
    for (tf_it = group->begin(); tf_it != group->end(); ++tf_it) {
      Logistic * tf = (Logistic*) (*tf_it)->inner();
      
      for (int i = 0; i < n_betas; ++i)
        tf->_betas[i] = betas[i];
    }
    
    delete piter;
    delete[] betas;
  }
  
  virtual bool setCovarSlots(int * slots, int length) {
    if (length != 1)
      return false;
    if (*slots < 0)
      return false;
    _covar_slot = *slots;
    return true;
  }

private:
  int _covar_slot;
  LogSum * _logsum;
  double * _betas;
  int * _rev_target_map;
  bool _is_fixed;
  int _max_iters;
  double _tolerance;
  
  struct opt_data {
    std::vector<TransitionFunction*> * group;
    TransitionPosteriorIterator * piter;
    LogSum * logsum;
    int covar_slot;
    int n_targets;
  };
  
  static double optfunc(int n, double * betas, void * udata) {
    struct opt_data * data = (struct opt_data*) udata;
    double result = 0;
    LogSum * logsum = data->logsum;
    unsigned int size = data->group->size();
    int n_targets = data->n_targets;
    TransitionPosteriorIterator * piter = data->piter;
    
    piter->reset();
    
    do {
      double x = piter->covar(data->covar_slot);
      double sum;
      
      logsum->clear();
      logsum->store(0);
      for (int i = 0; i < n_targets - 1; ++i)
        logsum->store(betas[i*2] + betas[i*2 + 1] * x);
      
      sum = logsum->compute();
      
      for (unsigned int gidx = 0; gidx < size; ++gidx) {
        for (int tgt_idx = 0; tgt_idx < n_targets; ++tgt_idx) {
          double post = piter->posterior(gidx, tgt_idx);
          result += post * ((*logsum)[tgt_idx] - sum);
        }
      }
    } while (piter->next());
    
    return -result;
  }
};

#endif
