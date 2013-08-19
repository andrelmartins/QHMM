#ifndef BASE_CLASSES_HPP
#define BASE_CLASSES_HPP

#include "iter.hpp"
#include "params.hpp"

#include <cmath>
#include <limits>

enum SuffStatType {
  Fixed,       /* no EM work, parameters are fixed */
  SinglePass,  /* only a single pass over posteriors per EM step */
  MultiPass    /* multiple passes over posteriors per EM step */
};

class EmissionSuffStat {
public:
  virtual int state() const = 0;
  virtual void reset() = 0;
  virtual void accum(Iter & iter, int slot, double * posterior) = 0;

  virtual SuffStatType suff_stats_type() { return SinglePass; }
};

bool inline same_probability(double a, double b) {
  return fabs(a - b) < std::numeric_limits<double>::epsilon();
}

class TransitionFunction {
  public:
    virtual bool validParams(Params const & params) const { return true; }
    virtual Params * getParams() const { return NULL; }
    virtual void setParams(Params const & params) {};
    virtual double log_probability(int target) const = 0;
    virtual double log_probability(Iter const & iter, int target) const = 0;
    virtual ~TransitionFunction() {};
};

class EmissionFunction {
  public:
    virtual bool validParams(Params const & params) const { return true; }
    virtual Params * getParams() const { return NULL; }
    virtual void setParams(Params const & params) {};
    virtual double log_probability(Iter const & iter, int slot) const = 0;
    virtual ~EmissionFunction() {};

    virtual EmissionSuffStat * suff_stats_instance() const = 0;
};

class MissingEmissionFunction : EmissionFunction {
  public:
    MissingEmissionFunction(EmissionFunction * func) : _func(func) {}
    ~MissingEmissionFunction() { // TODO: review memory management responsabilities
      delete _func;
    }

    bool validParams(Params const & params) const {
      return _func->validParams(params);
    }
  
    void setParams(Params const & params) {
      _func->setParams(params);
    }
  
    double log_probability(Iter const & iter, int slot) const {
      if (iter.is_missing(slot))
        return 0.0; /* log(1) */
      return _func->log_probability(iter, slot);
    }
  
  private:
    EmissionFunction * _func;
};

#endif
