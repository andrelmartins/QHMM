#ifndef BASE_CLASSES_HPP
#define BASE_CLASSES_HPP

#include "iter.hpp"

class TransitionFunction {
	public:
	  virtual double log_probability(int target) const = 0;
		virtual double log_probability(Iter const & iter, int target) const = 0;
    virtual ~TransitionFunction() {};
};

class EmissionFunction {
  public:
    virtual double log_probability(Iter const & iter) const = 0;
    virtual ~EmissionFunction() {};
};

#endif
