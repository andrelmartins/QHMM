#ifndef NEGBINOMIAL_HPP
#define NEGBINOMIAL_HPP

#include <base_classes.hpp>
#include <em_base.hpp>

#include "../src/math.hpp"


class NegativeBinomial : public EmissionFunction {
public:
  NegativeBinomial(int stateID, int slotID, double mean = 1.0, double dispersion = 1.0) : EmissionFunction(stateID, slotID), _mean(mean), _dispersion(dispersion), _fixedParams(false), _offset(0), _tolerance(1e-6), _maxIter(100), _tblSize(64), _momInit(false) {

    _logp_tbl = new double[_tblSize];
    update_logp_tbl();
  }

  ~NegativeBinomial() {
    delete[] _logp_tbl;
  }
  
  
  virtual bool validParams(Params const & params) const {
    // supports two parameters: mean, dispersion
    // all must be > 0
    if (params.length() != 2)
      return false;
    
    for (int i = 0; i < params.length(); ++i)
      if (params[i] <= 0)
        return false;

    // TODO: check fixness
    return true;
  }
  
  virtual Params * getParams() const {
    Params * params;
    double pvals[2] = {_mean, _dispersion };
    params = new Params(2, pvals);

    if (_fixedParams) {
      params->setFixed(0, true);
      params->setFixed(1, true);
    }
    return params;
  }
  
  virtual void setParams(Params const & params) {
    _mean = params[0];
    _dispersion = params[1];
    _fixedParams = params.isFixed(0) || params.isFixed(1); // TODO: fix hack in validParams ...
    update_logp_tbl();
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "offset")) {
      *out_value = (double) _offset;
      return true;
    } else if(!strcmp(name, "maxIter")) {
      *out_value = _maxIter;
      return true;
    } else if (!strcmp(name, "tolerance")) {
      *out_value = _tolerance;
      return true;
    } else if (!strcmp(name, "tblSize")) {
      *out_value = _tblSize;
      return true;
    } else if (!strcmp(name, "momInit")) {
      *out_value = (_momInit ? 1 : 0);
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "offset")) {
      _offset = value;
      return true;
    } else if (!strcmp(name, "maxIter")) {
      int maxIter = (int) value;
      if (maxIter <= 0) {
        log_msg("maxIter must be > 0: %d\n", maxIter);
        return false;
      }
      _maxIter = maxIter;
      return true;
    } else if (!strcmp(name, "tolerance")) {
      if (value < 0) {
        log_msg("tolerance must be >= 0: %g\n", value);
        return false;
      }
      _tolerance = value;
      return true;
    } else if (!strcmp(name, "tblSize")) {
      int tblSize = (int)value;
      if (tblSize <= 0)
        _tblSize = tblSize; /* this just disables tbl use */
      else {
        _tblSize = tblSize;
        delete _logp_tbl;
        _logp_tbl = new double[_tblSize];
        update_logp_tbl();
      }
      return true;
    } else if (!strcmp(name, "momInit")) {
      _momInit = (value != 0.0);
      return true;
    }
    return false;
  }
  
  virtual double log_probability(Iter const & iter) const {
    int x = (int) (iter.emission(_slotID) + _offset);

    assert(x >= 0);

    if (x < _tblSize)
      return _logp_tbl[x];
    return logprob(x);
  }
  
  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (_fixedParams)
      return;

    // in these computations we revert back to standard (r, p) parameterization
    // as it makes the math easier
    double r = _dispersion;
    
    // sufficient statistics
    double sum_Pzi = 0;
    double sum_Pzi_xi = 0;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      NegativeBinomial * ef = (NegativeBinomial*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int x = (int) (iter.emission(ef->_slotID) + _offset);
          
          sum_Pzi += post_j[j];
          sum_Pzi_xi += post_j[j] * x;
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // update parameter
    // 1. estimate 'r' (dispersion)
    // 1.1 Apply Newton's method
    double r_prev = r_start_value(r, sum_Pzi, sum_Pzi_xi, sequences, group);
    double change;
    int i = 0;
    int reductionFactor = 2; /* how much to reduce the starting dispersion */
    do {
      ++i;
      r = r_prev - newton_ratio(sum_Pzi, sum_Pzi_xi, r_prev, sequences, group);
      
      /* test boundary conditions */
      if (QHMM_isinf(r) || QHMM_isnan(r)) {
        log_state_slot_msg(_stateID, _slotID, "dispersion update failed: %g (keeping old value: %g)\n", r, _dispersion);

        r = _dispersion;
        break;
      } else if (r <= 0) {
        /* we went to a very large value and then stepped back too much
         * instead try to step back to a fraction of the current parameter
         */
        if (r_prev > _dispersion) {
          log_state_slot_msg(_stateID, _slotID, "dispersion lower bound hit: %g (using %g)\n", r, _dispersion / reductionFactor);
          r = _dispersion / reductionFactor;
          r_prev = r;
          reductionFactor = reductionFactor * reductionFactor;
          continue;
        }
        
        /* otherwise just clamp value by tolerance */
        log_state_slot_msg(_stateID, _slotID, "dispersion lower bound hit: %g (using %g)\n", r, _tolerance);
        r = _tolerance;
        r_prev = _tolerance;
        continue;
      }

      change = fabs(r - r_prev);
      r_prev = r;
    } while (change > _tolerance && i < _maxIter);

    // 1.2 check for weirdness
    if (r > 1000 || QHMM_isnan(r) || QHMM_isinf(r)) {  // todo: fix magical value!!
      log_state_slot_msg(_stateID, _slotID, "dispersion update failed: %g (keeping old value: %g)\n", r, _dispersion);
      return;
    }

    // 2. estimate 'p'
    double p = sum_Pzi_xi / (sum_Pzi * r + sum_Pzi_xi);
    
    // accept update
    _mean = (p * r) /  (1.0 - p);
    _dispersion = r;
    
    update_logp_tbl();

    // propagate to other elements in the group
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      NegativeBinomial * ef = (NegativeBinomial*) (*ef_it)->inner();
      
      if (ef != this) {
        ef->_mean = _mean;
        ef->_dispersion = _dispersion;
        ef->copy_logp_tbl(this);
      }
    }
  }
  
private:
  double _mean; // := m
  double _dispersion; // := r
  bool _fixedParams;
  double _tolerance;
  double _offset;
  int _maxIter;
  int _tblSize;
  bool _momInit;
  
  // these are set in update_logp_tbl & copied by copy_logp_tbl
  double _A1; // := r [log(r) - log(r + m)]
  double _A2; // := [log(m) - log(r + m)]
  double _A3; // := log Gammafn(r)
  double * _logp_tbl;

  double logprob(int x) const {
    // TODO: check if computing log GammaFn[r + x] is faster or slower than:
    //       log GamamFn[r] + sum_{a=1}^x log(r + a - 1)
    return _A1 - _A3 + x * _A2 + QHMM_log_gamma(_dispersion + x) - QHMM_log_gamma(x + 1);// log(x!)
  }
  
  void update_logp_tbl() {
    _A1 = _dispersion * (log(_dispersion) - log(_dispersion + _mean));
    _A2 = log(_mean) - log(_dispersion + _mean);
    _A3 = QHMM_log_gamma(_dispersion);
    
    for (int i = 0; i < _tblSize; ++i)
      _logp_tbl[i] = logprob(i);
  }

  void copy_logp_tbl(NegativeBinomial * other) {
    other->_A1 = _A1;
    other->_A2 = _A2;
    other->_A3 = _A3;
    if (_tblSize > 0)
      memcpy(_logp_tbl, other->_logp_tbl, _tblSize * sizeof(double));
  }
  
  double r_start_value(double prev_r, double sum_Pzi, double sum_Pzi_xi, EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (!_momInit)
      return prev_r;
    
    // estimate variance
    double mean = sum_Pzi_xi / sum_Pzi;
    double sum_Pzi_sqdiff = 0.0;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      NegativeBinomial * ef = (NegativeBinomial*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int x = (int) (iter.emission(ef->_slotID) + _offset);
          
          sum_Pzi_sqdiff += post_j[j] * (x - mean) * (x - mean);
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    //
    double var = sum_Pzi_sqdiff / sum_Pzi;
    double r_est = fabs(mean / (var - mean));
    
    if (r_est > 1000) // todo: fix magical value!!
      return 500;
    return r_est;
  }
  
  double newton_ratio(double A, double B, double r, EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    
    // constant terms
    double const_num = 0;
    double const_denom = 0;
       
    const_num = -QHMM_digamma(r) + log(A * r) - log(A * r + B);
    const_denom = -QHMM_trigamma(r) + B / (r * (A * r + B));
    
    // data dependent terms
    double sum_num = 0;
    double sum_denom = 0;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      NegativeBinomial * ef = (NegativeBinomial*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          double x = (iter.emission(ef->_slotID) + _offset);
          
          sum_num += post_j[j] * QHMM_digamma(x + r);
          sum_denom += post_j[j] * QHMM_trigamma(x + r);
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // TODO: check if some trickery with the GammaFn can help here!
    
    // ratio
    double f_r = sum_num/A + const_num;
    double g_r = sum_denom/A + const_denom;
    
    return f_r / g_r;
  }
  
};

#endif
