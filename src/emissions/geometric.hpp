#ifndef GEOMETRIC_H
#define GEOMETRIC_H

#include "../base_classes.hpp"
#include "../em_base.hpp"

class Geometric : public EmissionFunction {
public:
  Geometric(int stateID, int slotID, int base = 0, double prob = 0.5) : EmissionFunction(stateID, slotID), _base(base), _log_prob(log(prob)), _log_1_prob(log(1 - prob)), _is_fixed(false) {}

  virtual bool validParams(Params const & params) const {
    return params.length() == 1 && params[0] > 0 && params[0] <= 1;
  }

  virtual Params * getParams() const {
    double prob = exp(_log_prob);
    Params * params = new Params(1, &prob);
    params->setFixed(0, _is_fixed);
    return params;
  }
  
  virtual void setParams(Params const & params) {
    double prob = params[0];
    _log_prob = log(prob);
    _log_1_prob = log(1 - prob);
    _is_fixed = params.isFixed(0);
  }
  
  virtual bool getOption(const char * name, double * out_value) {
    if (!strcmp(name, "base")) {
      *out_value = _base;
      return true;
    }
    return false;
  }

  virtual bool setOption(const char * name, double value) {
    if (!strcmp(name, "base")) {
      int base = (int) value;
      if (base != 0 && base != 1) // TODO: add warning message
	return false;
      _base = base;
    }
    return false;
  }

  virtual double log_probability(Iter const & iter) const {
    int x = (int) iter.emission(_slotID); // cast to integer
    
    // log prob(x) = log( (1 - prob)^(x - base) * prob )
    //             = (x - base) * log(1 - prob) + log(prob)
    return (x - _base) * _log_1_prob + _log_prob;
  }

  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (_is_fixed)
      return;

    // sufficient statistics are the sum of the state posteriors and the sum of the posterior times
    // the observations (less the base)
    double sum_Pzi = 0;
    double sum_Pzi_xi = 0;
      
    std::vector<EmissionFunction*>::iterator ef_it;
      
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      Geometric * ef = (Geometric*) (*ef_it)->inner();
      PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
      
      do {
        const double * post_j = post_it->posterior();
        Iter & iter = post_it->iter();
        iter.resetFirst();
        
        for (int j = 0; j < iter.length(); iter.next(), ++j) {
          int x = (int) iter.emission(ef->_slotID) - _base;
          
          sum_Pzi += post_j[j];
          sum_Pzi_xi += post_j[j] * x;
        }
      } while (post_it->next());
      
      delete post_it;
    }
    
    // use expected counts to estimate parameter value
    double mean_xi = sum_Pzi_xi/sum_Pzi;
    double prob = 1.0/(1.0 + mean_xi);
    _log_prob = log(prob);
    _log_1_prob = log(1 - prob);

    // propagate to other elements in the group
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      Geometric * ef = (Geometric*) (*ef_it)->inner();
        
      if (ef != this) {
        ef->_log_prob = _log_prob;
        ef->_log_1_prob = _log_1_prob;
      }
    }
  }
    
private:
  int _base;
  double _log_prob;
  double _log_1_prob;
  bool _is_fixed;
};

#endif
