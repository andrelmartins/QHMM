#ifndef BASE_CLASSES_HPP
#define BASE_CLASSES_HPP

#include "iter.hpp"
#include "params.hpp"

#include <cmath>
#include <limits>

// forward declaration to avoid include loop
class EMSequences;

bool inline same_probability(double a, double b) {
  double diff = a - b;
  if (diff < 0.0)
    diff = -diff;
  return diff <= std::numeric_limits<double>::epsilon();
}

class TransitionFunction {
public:
  TransitionFunction(int n_states, int stateID, int n_targets, int * targets) : _stateID(stateID), _n_states(n_states), _n_targets(n_targets) {
    if (n_targets == 0)
      _targets = NULL;
    else {
      _targets = new int[_n_states];
    
      for (int i = 0; i < _n_targets; ++i)
	_targets[i] = targets[i];
    }
  }

  virtual ~TransitionFunction() {
    if (_targets != NULL)
      delete[] _targets;
  }

  virtual bool validParams(Params const & params) const { return true; }
  virtual Params * getParams() const { return NULL; }
  virtual void setParams(Params const & params) {};
  virtual bool setCovarSlots(int * slots, int length) { return false; } // default is no slots
  virtual bool getOption(const char * name, double * out_value) { return false; } // no options
  virtual bool setOption(const char * name, double value) { return false; } // no options
  
  virtual double log_probability(int target) const = 0;
  virtual double log_probability(Iter const & iter, int target) const = 0;

  int stateID() const { return _stateID; }
  int n_targets() const { return _n_targets; }
  const int * targets() const { return _targets; }

  virtual void updateParams(EMSequences * sequences, std::vector<TransitionFunction*> * group) {}

protected:
  const int _stateID;
  const int _n_states;
  const int _n_targets;
  int * _targets;
};

class EmissionFunction {
public:

  EmissionFunction(int stateID, int slotID) : _stateID(stateID), _slotID(slotID) {}

  virtual bool validParams(Params const & params) const { return true; }
  virtual Params * getParams() const { return NULL; }
  virtual void setParams(Params const & params) {};
  virtual bool setCovarSlots(int * slots, int length) { return false; } // default is no slots
  virtual bool getOption(const char * name, double * out_value) { return false; } // no options
  virtual bool setOption(const char * name, double value) { return false; } // no options
  
  virtual double log_probability(Iter const & iter) const = 0;
  virtual ~EmissionFunction() {};

  int stateID() { return _stateID; }
  int slotID() { return _slotID; }
  
  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {}

protected:
  const int _stateID;
  const int _slotID;
};

class MissingEmissionFunction : public EmissionFunction {
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
  
    bool setCovarSlots(int * slots, int length) {
      return _func->setCovarSlots(slots, length);
    }
  
    bool getOption(const char * name, double * out_value) {
      return _func->getOption(name, out_value);
    }
  
    bool setOption(const char * name, double value) {
      return _func->setOption(name, value);
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
