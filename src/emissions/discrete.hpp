#ifndef DISCRETE_EMISSION_HPP
#define DISCRETE_EMISSION_HPP

#include <cmath>
#include <limits>
#include "../base_classes.hpp"

class DiscreteEmissions : public EmissionFunction {
public:
  DiscreteEmissions(int stateID, int slotID, int offset = 1) : EmissionFunction(stateID, slotID), _offset(offset), _log_probs(NULL), _alphabetSize(0) {}
  ~DiscreteEmissions() {
    if (_log_probs)
      delete[] _log_probs;
  }

  virtual bool validParams(Params const & params) const {
    double sum = 0.0;

    for (int i = 0; i < params.length(); ++i)
      sum += params[i];

    return params.length() > 0 && same_probability(sum, 1.0);
  }

  virtual Params * getParams() const {
    double * probs = new double[_alphabetSize];

    for (int i = 0; i < _alphabetSize; ++i)
      probs[i] = exp(_log_probs[i]);

    Params * result = new Params(_alphabetSize, probs);

    delete probs;

    return result;
  }

  virtual void setParams(Params const & params) {
    if (_alphabetSize != params.length()) {
      delete[] _log_probs;
      _alphabetSize = params.length();
      _log_probs = new double[_alphabetSize];
    }

    for (int i = 0; i < _alphabetSize; ++i)
      _log_probs[i] = log(params[i]);
  }

  virtual double log_probability(Iter const & iter) const {
    int x = (int) iter.emission(_slotID); // cast to integer
    int y = x - _offset;
    
    if (y < 0 || y >= _alphabetSize)
      return -std::numeric_limits<double>::infinity();

    return _log_probs[y];
  }

private:
  int _offset;
  int _alphabetSize;
  double * _log_probs;

};

#endif
