#ifndef BASE_CLASSES_HPP
#define BASE_CLASSES_HPP

#include "iter.hpp"
#include "params.hpp"

#include <cmath>
#include <limits>

bool inline same_probability(double a, double b) {
  return fabs(a - b) < std::numeric_limits<double>::epsilon();
}

class TransitionFunction {
public:
  TransitionFunction(int stateID) : _stateID(stateID) {}

  virtual bool validParams(Params const & params) const { return true; }
  virtual Params * getParams() const { return NULL; }
  virtual void setParams(Params const & params) {};
  virtual double log_probability(int target) const = 0;
  virtual double log_probability(Iter const & iter, int target) const = 0;
  virtual ~TransitionFunction() {};

  int stateID() { return _stateID; }

protected:
  const int _stateID;
};

class EmissionFunction {
public:

  EmissionFunction(int stateID, int slotID) : _stateID(stateID), _slotID(slotID) {}

  virtual bool validParams(Params const & params) const { return true; }
  virtual Params * getParams() const { return NULL; }
  virtual void setParams(Params const & params) {};
  virtual double log_probability(Iter const & iter) const = 0;
  virtual ~EmissionFunction() {};

  int stateID() { return _stateID; }
  int slotID() { return _slotID; }

protected:
  const int _stateID;
  const int _slotID;
};

class MissingEmissionFunction : EmissionFunction {
  public:
  MissingEmissionFunction(EmissionFunction * func) : EmissionFunction(func->stateID(), func->slotID()), _func(func) {}
    ~MissingEmissionFunction() { // TODO: review memory management responsabilities
      delete _func;
    }

    bool validParams(Params const & params) const {
      return _func->validParams(params);
    }
  
    void setParams(Params const & params) {
      _func->setParams(params);
    }
  
    double log_probability(Iter const & iter) const {
      if (iter.is_missing(_slotID))
        return 0.0; /* log(1) */
      return _func->log_probability(iter);
    }
  
  private:
    EmissionFunction * _func;
};

#endif
