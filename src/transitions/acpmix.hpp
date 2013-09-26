#ifndef ACPMIX_HPP
#define ACPMIX_HPP

#include <stdexcept>
#include <cmath>
#include <cstring>
#include <limits>
#include "../base_classes.hpp"
#include "../em_base.hpp"
#include "../math.hpp"

// Mixture of Auto-correlation and prior (from covar)
//
// Parameters:
//   alpha: auto-correlation
//   gamma: auto-correlation mixture weight [always fixed]
//
// Covar slot must have dim == n_targets
// Covars are in log-space
class ACPMix : public TransitionFunction {
public:

  ACPMix(int n_states, int stateID, int n_targets, int * targets, double alpha = 0.5, double gamma = 0.5, int covar_slot = 0, int max_iters = 100, double tolerance = 1e-4) : TransitionFunction(n_states, stateID, n_targets, targets), _covar_slot(covar_slot), _alpha(alpha), _gamma(gamma), _is_fixed_alpha(false), _max_iters(100), _tolerance(tolerance) {
    _valid_states = new bool[n_states];
    
    // set all to false
    for (int i = 0; i < n_states; ++i)
      _valid_states[i] = false;
    
    // valid states
    for (int i = 0; i < n_targets; ++i)
      _valid_states[targets[i]] = true;

    // create reverse target map
    _rev_target_map = new int[n_states];
    // set all to -1
    for (int i = 0; i < n_states; ++i)
      _rev_target_map[i] = -1;
    // set targets
    for (int i = 0; i < n_targets; ++i)
      _rev_target_map[targets[i]] = i;

    update_log_probs();
  }

  ~ACPMix() {
    delete[] _valid_states;
    delete[] _rev_target_map;
  }

  virtual bool validParams(Params const & params) const {
    return params.length() == 2 && params[0] >= _tolerance && params[0] <= 1 - _tolerance && params[1] >= 0 && params[1] <= 1;
  }

  virtual Params * getParams() const {
    double pars[2] = { _alpha, _gamma };
    Params * result = new Params(2, pars);
    if (_is_fixed_alpha)
      result->setFixed(0, true);
    result->setFixed(1, true); // always fixed
    return result;
  }

  virtual void setParams(Params const & params) {
    _alpha = params[0];
    _gamma = params[1];
    update_log_probs();
    _is_fixed_alpha = params.isFixed(0) || _gamma == 0.0;
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
	log_msg("invalid maxIters: %d : should be > 0\n",
		iters);
	return false;
      }
      _max_iters = iters;
      return true;
    }
    if (!strcmp(name, "tolerance")) {
      if (value <= 0 || value > _alpha) {
	log_msg("invalid tolerance: %g : should be > 0 & <= alpha = %g\n",
		value, _alpha);
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
    if (_valid_states[target]) {
      int target_idx = _rev_target_map[target];
      double log_prior = iter.covar_i(_covar_slot, target_idx);

      if (target == _stateID)
	return QHMM_logsum(_log_acor_self, log_prior + _log_prior_weight);
      else
	return QHMM_logsum(_log_acor_other, log_prior + _log_prior_weight);
    } else
      return -std::numeric_limits<double>::infinity();
  }
  
  virtual void updateParams(EMSequences * sequences, std::vector<TransitionFunction*> * group) {
    if (_is_fixed_alpha)
      return;

    TransitionPosteriorIterator * piter = sequences->transition_iterator(*group);
    /* use Newton's method to fit alpha */
    double alpha = _alpha;
    bool hit_edge = false;
    for (int iters = _max_iters; iters > 0; --iters) {
      int sign;

      /* compute function and derivative */
      double fx = 0;
      double gx = 0;
      
      compute_fx_gx(alpha, piter, group, &fx, &gx);

      if (QHMM_isinf(gx) || QHMM_isinf(gx)) {
	log_state_msg(_stateID, "alpha update failed: iter alpha: %g prev alpha: %g\n", alpha, _alpha);
	delete piter;
	return;
      }

      if (fabs(fx) < _tolerance)
	break;

      /* update */
      sign = ratio_sign(fx, gx);
      alpha = alpha - sign * exp(log(fabs(fx)) - log(fabs(gx)));

      /* impose barriers */
      if (alpha < _tolerance) {
	hit_edge = true;
	alpha = _tolerance;
      }
      if (alpha > 1 - _tolerance) {
	hit_edge = true;
	alpha = 1 - _tolerance;
      }
    }

    if (hit_edge)
      log_state_msg(_stateID, "alpha update hit edge\n");

    // propagate parameters
    std::vector<TransitionFunction*>::iterator tf_it;
    for (tf_it = group->begin(); tf_it != group->end(); ++tf_it) {
      ACPMix * tf = (ACPMix*) (*tf_it)->inner();

      tf->update_log_probs(alpha);
    }

    delete piter;
  }


  // TODO: enforce covar slot dimension == n_targets
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
  double _alpha;
  double _gamma;
  bool _is_fixed_alpha;
  bool * _valid_states;
  int * _rev_target_map;
  double _log_acor_self;
  double _log_acor_other;
  double _log_prior_weight;
  int _max_iters;
  double _tolerance;

  void update_log_probs(double alpha) {
    _alpha = alpha;
    update_log_probs();
  }

  void update_log_probs() {
    _log_acor_self = log(_alpha) + log(_gamma);
    _log_acor_other = log(1.0 - _alpha) + log(_gamma) - log(_n_targets - 1);
    _log_prior_weight = log(1.0 - _gamma);
  }

  void compute_fx_gx(double alpha, TransitionPosteriorIterator * piter, std::vector<TransitionFunction*> * group, double * out_fx, double * out_gx) {
    double fx = 0;
    double gx = 0;

    do {
      // NOTE: assume all states in the group have the save _covar_slot

      for (int tgt_idx = 0; tgt_idx < _n_targets; ++tgt_idx) {
	double log_prior = piter->covar_i(_covar_slot, tgt_idx);
	double prior = exp(log_prior);

	for (unsigned int gidx = 0; gidx < group->size(); ++gidx) {
	  double post = piter->posterior(gidx, tgt_idx);

	  /* NOTE: Since this is intended to compute the ratio fx/gx
	   *       I'm simplifying the expressions by dividing both by gamma.
	   */
	  ACPMix * gState = (ACPMix*) (*group)[gidx]->inner();

	  bool is_self = gState->_stateID == gState->_targets[tgt_idx];
	  if (is_self) {
	    double denom = (_gamma * alpha + (1.0 - _gamma) * prior);
	    fx += post / denom;
	    gx -= post * _gamma / (denom * denom);
	  } else {
	    double denom = (_gamma * (1.0 - alpha) + (_n_targets - 1) * (1.0 - _gamma) * prior);
	    fx -= post / denom;
	    gx -= post * _gamma / (denom * denom);
	  }

	  if (QHMM_isinf(gx) || QHMM_isinf(gx)) {
	    log_state_msg(gState->_stateID, "alpha iter failed: fx: %g gx:%g prior: %g post: %g\n", fx, gx, prior, post);
	    *out_fx = fx;
	    *out_gx = gx;
	    
	    return;
	  }
	}
      }
    } while (piter->next());

    /* update output values */
    *out_fx = fx;
    *out_gx = gx;
  }

  int ratio_sign(double a, double b) {
    if (a > 0)
      return (b > 0 ? 1 : -1);
    return (b > 0 ? - 1: 1);
  }
};

#endif
