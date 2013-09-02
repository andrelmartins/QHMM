#ifndef DISCRETE_EMISSION_HPP
#define DISCRETE_EMISSION_HPP

#include <cmath>
#include <limits>
#include <cstring>
#include "../base_classes.hpp"
#include "../em_base.hpp"

class DiscreteEmissions : public EmissionFunction {
public:
  DiscreteEmissions(int stateID, int slotID, int offset = 1) : EmissionFunction(stateID, slotID), _offset(offset), _alphabetSize(0), _log_probs(NULL) {}
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
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "offset")) {
      *out_value = (double) _offset;
      return true;
    }
    return false;
  }
  
  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "offset")) {
      _offset = (int) value;
      return true;
    }
    return false;
  }

  virtual double log_probability(Iter const & iter) const {
    int x = (int) iter.emission(_slotID); // cast to integer
    int y = x - _offset;
    
    if (y < 0 || y >= _alphabetSize)
      return -std::numeric_limits<double>::infinity();

    return _log_probs[y];
  }

  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    // sufficient statistics are the per symbol expected counts
    double expected_counts[_alphabetSize];
    
    std::vector<EmissionFunction*>::iterator ef_it;
    
    for (int i = 0; i < _alphabetSize; ++i)
      expected_counts[i] = 0;
    
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      DiscreteEmissions * ef = (DiscreteEmissions*) *ef_it;
      EMSequences::PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);

      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int symbol = (int) iter.emission(ef->_slotID) - _offset;
          
          expected_counts[symbol] += post_j[j];
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // use expected counts to estimate parameter values
    double normalization = 0;
    for (int i = 0; i < _alphabetSize; ++i)
      normalization += expected_counts[i];
    for (int i = 0; i < _alphabetSize; ++i)
      _log_probs[i] = log(expected_counts[i] / normalization);
    
    // propagate to other elements in the group
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      DiscreteEmissions * ef = (DiscreteEmissions*) *ef_it;
      
      if (ef != this)
        memcpy(ef->_log_probs, _log_probs, _alphabetSize * sizeof(double));
    }
  }
  
private:
  int _offset;
  int _alphabetSize;
  double * _log_probs;

};

#endif
