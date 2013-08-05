#ifndef BASE_CLASSES_HPP
#define BASE_CLASSES_HPP

#include "iter.hpp"
#include "params.hpp"

class TransitionFunction {
	public:
    virtual bool validParams(Params const & params) { return true; }
    virtual void setParams(Params const & params)  {}
    virtual double log_probability(int target) const = 0;
		virtual double log_probability(Iter const & iter, int target) const = 0;
    virtual ~TransitionFunction() {};
};

class EmissionFunction {
  public:
    virtual bool validParams(Params const & params) { return true; }
    virtual void setParams(Params const & params)  {}
    virtual double log_probability(Iter const & iter, int slot) const = 0;
    virtual ~EmissionFunction() {};
};

class MissingEmissionFunction : EmissionFunction {
  public:
    MissingEmissionFunction(EmissionFunction * func) : _func(func) {}
    ~MissingEmissionFunction() { // TODO: review memory management responsabilities
      delete _func;
    }

    bool validParams(Params const & params) {
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
