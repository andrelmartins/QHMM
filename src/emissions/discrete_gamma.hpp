#ifndef DISCRETE_GAMMA_HPP
#define DISCRETE_GAMMA_HPP

#include <base_classes.hpp>
#include <em_base.hpp>

#include "../src/math.hpp"


class DiscreteGamma : public EmissionFunction {
public:
  DiscreteGamma(int stateID, int slotID, double shape = 1.0, double scale = 2.0) : EmissionFunction(stateID, slotID), _shape(shape), _scale(scale), _scale_private(1.0), _final_scale(scale), _fixedMean(false), _mean(shape * scale), _fixedParams(false), _offset(0), _shift(1.0), _tolerance(1e-6), _maxIter(100), _tblSize(64) {

    _logp_tbl = new double[_tblSize];
    update_logp_tbl();
  }

  ~DiscreteGamma() {
    delete[] _logp_tbl;
  }
  
  virtual bool validParams(Params const & params) const {
    // supports two or three parameters: shape, scale, mean
    // all three must be > 0
    if (params.length() < 2 || params.length() > 3)
      return false;
    
    for (int i = 0; i < params.length(); ++i)
      if (params[i] <= 0)
        return false;
    
    if (params.length() == 3) {
      double shape = params[0];
      double scale = params[1];
      double mean = params[2];
      double aux = shape * scale;
      
      if (fabs(aux - mean) > 10*std::numeric_limits<double>::epsilon()) {
        log_msg("shape * scale must be equal to mean: %g != %g\n", aux, mean);
        return false;
      }
    }
    return true;
  }
  
  virtual Params * getParams() const {
    Params * params;
    if (_fixedMean) {
      double pvals[3] = {_shape, _scale, _mean };
      params = new Params(3, pvals);
      params->setFixed(2, true); // mean is always fixed
    } else {
      double pvals[2] = {_shape, _scale };
      params = new Params(2, pvals);
    }
    if (_fixedParams) {
      params->setFixed(0, true);
      params->setFixed(1, true);
    }
    return params;
  }
  
  virtual void setParams(Params const & params) {
    if (params.length() == 3) {
      _mean = params[2];
      _scale_private = 1.0; // overide private scale
      
      _fixedMean = true;
    } else
      _fixedMean = false;
    
    _fixedParams = params.isFixed(0) || params.isFixed(1); // TODO: fix hack in validParams ...
    update_shape_scale(params[0], params[1]);
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "offset")) {
      *out_value = (double) _offset;
      return true;
    } else if(!strcmp(name, "shift")) {
      *out_value = _shift;
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
    } else if (!strcmp(name, "scale_private")) {
      *out_value = _scale_private;
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "offset")) {
      _offset = value;
      return true;
    } else if (!strcmp(name, "shift")) {
      if (value < 0) {
        log_msg("shift must be >= 0: %g\n", value);
        return false;
      }
      _shift = value;
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
    } else if (!strcmp(name, "scale_private")) {
      if (_fixedMean) {
        log_msg("scale_private is not a valid option with a fixed mean\n");
        return false;
      }
      if (value <= 0) {
        log_msg("scale_private must be > 0: %g\n", value);
        return false;
      }
      update_scale_private(value);
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
    
    // sufficient statistics
    double sum_Pzi = 0;
    double sum_Pzi_xi = 0;
    double sum_Pzi_log_xi = 0;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      DiscreteGamma * ef = (DiscreteGamma*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      double scale_i = ef->_scale_private;
      double log_scale_i = log(scale_i);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int x = (int) (iter.emission(ef->_slotID) + _offset);
          
          sum_Pzi += post_j[j];
          sum_Pzi_xi += post_j[j] * x / scale_i;
          sum_Pzi_log_xi += post_j[j] * (log(x) - log_scale_i);
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // update parameter
    // 1. estimate shape
    double mean = sum_Pzi_xi / sum_Pzi;
    double s = log(mean) - sum_Pzi_log_xi / sum_Pzi;

    //log_state_slot_msg(_stateID, _slotID, " U: mean = %g, sum_Pzi_xi = %g, sum_Pzi = %g sum_Pzi_log_xi = %g\n", mean, sum_Pzi_xi, sum_Pzi, sum_Pzi_log_xi);

    // 1.1 Initial guess
    double shape = (3 - s + sqrt(pow(s - 3, 2) + 24 * s)) / (12 * s);
    
    if (QHMM_isinf(shape) || QHMM_isnan(shape) || shape <= 0) {
      log_state_slot_msg(_stateID, _slotID, "initial shape guess failed: %g (starting with old value: %g)\n", shape, _shape);
      shape = _shape;
    }
    
    // 1.2 Apply Newton's method to refine estimate
    double ki;
    double knext = shape;
    int i = 0;
    do {
      ++i;
      ki = shape;
      /*log_state_slot_msg(_stateID, _slotID, " U[%d] s = %g, ki = %g\n",
			 i, s, ki);
      */
      
      /* update shape */
      //knext = ki - (log(ki) - digamma(ki) - s) / (1.0/ki - trigamma(ki));
      knext = ki - (log(ki) - QHMM_digamma(ki) - s) / (1.0 / ki - QHMM_trigamma(ki));
      
      /* test boundary conditions */
      if (QHMM_isinf(shape) || QHMM_isnan(shape) || shape <= 0) {
        log_state_slot_msg(_stateID, _slotID, "shape update failed: %g (keeping old value: %g)\n", knext, _shape);

        knext = shape;
        break;
      }
    } while (fabs(ki - shape) > _tolerance && i < _maxIter);

    shape = knext;
    
    // check for weirdness
    if (shape > 1000 || QHMM_isnan(shape) || QHMM_isinf(shape)) {
      log_state_slot_msg(_stateID, _slotID, "shape update failed: %g (keeping old value: %g)\n", shape, _shape);
      return;
    }
    
    // 1.3 Accept update
    _shape = shape;
    
    // 2. Update scale
    if (_fixedMean)
      _scale = _mean / _shape;
    else
      _scale = mean / _shape;
    
    _final_scale = _scale_private * _scale;
    update_logp_tbl();

    // propagate to other elements in the group
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      DiscreteGamma * ef = (DiscreteGamma*) (*ef_it)->inner();
      
      if (ef != this) {
        ef->_shape = _shape;
        ef->_scale = _scale;
        ef->_final_scale = _final_scale;
        ef->copy_logp_tbl(this);
      }
    }
  }
  
private:
  double _shape;
  double _scale;
  bool _fixedMean;
  double _mean;
  bool _fixedParams;
  double _offset;
  double _shift;
  double _tolerance;
  int _maxIter;
  int _tblSize;
  
  double _scale_private;
  double _final_scale;

  double * _logp_tbl;
  
  double logprob(int x) const {
    // TODO: double check this math: I think to be correct this need to have an open interval ...
    //       P(X = x) = p(x + (1 - shift) <= X < x + shift) = p(X < x + shift) - p(X < x + (shift - 1))
    //
    //       in other words, prob of X in [x + shift - 1, x + shift [
    //
    //       which differs the result below by p(X = x + shift) - p(X = x + (shift - 1)) when useLowerTail = TRUE
    //       and by p(X = x + (shift - 1)) - p(X = x + shift) when useLowerTail = FALSE
    //
    double x_low = x + (_shift - 1.0);
    double x_high = x + _shift;
    
    if (x_low <= 0)
      return QHMM_log_gamma_cdf_lower(x_high, _shape, _final_scale);

    // upper tail is more stable
    return QHMM_logdiff(QHMM_log_gamma_cdf_upper(x_low, _shape, _final_scale), QHMM_log_gamma_cdf_upper(x_high, _shape, _final_scale));
    // return logdiff(log_gamma_cdf_lower(x_high, _shape, _final_scale), log_gamma_cdf_lower(x_low, _shape, _final_scale));
  }
  
  void update_logp_tbl() {
    for (int i = 0; i < _tblSize; ++i)
      _logp_tbl[i] = logprob(i);
  }

  void copy_logp_tbl(DiscreteGamma * other) {
    if (_tblSize > 0)
      memcpy(_logp_tbl, other->_logp_tbl, _tblSize * sizeof(double));
  }
  
  void update_scale(double scale) {
    _scale = scale;
    _final_scale = scale * _scale_private;
    update_logp_tbl();
  }
  
  void update_scale_private(double value) {
    _scale_private = value;
    _final_scale = _scale * _scale_private;
    update_logp_tbl();
  }
  
  void update_shape_scale(double shape, double scale) {
    _shape = shape;
    update_scale(scale);
  }
};

#endif
