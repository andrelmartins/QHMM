#ifndef NORMAL_HPP
#define NORMAL_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include "../base_classes.hpp"
#include "../em_base.hpp"

class Normal : public EmissionFunction {
public:

  Normal(int stateID, int slotID, double mean = 0.0, double var = 1.0) : EmissionFunction(stateID, slotID), _mean(mean), _var(var), _is_fixed_mean(false), _is_fixed_var(false) {}
 
  virtual bool validParams(Params const & params) const {
    return params.length() == 2 && params[1] > 0;
  }
  
  virtual Params * getParams() const {
    double * values = new double[2];
    
    values[0] = _mean;
    values[1] = _var;
    
    Params * result = new Params(2, values);
    result->setFixed(0, _is_fixed_mean);
    result->setFixed(1, _is_fixed_var);
    
    delete[] values;
    
    return result;
  }
  
  virtual void setParams(Params const & params) {
    _mean = params[0];
    _var = params[1];
    _is_fixed_mean = params.isFixed(0);
    _is_fixed_var = params.isFixed(1);
    
    update_consts();
  }
  
  virtual double log_probability(Iter const & iter) const {
    // log P(x | mu, var) = -log(sqrt(2 pi var)) - (x - mu)^2 / (2 var)
    int x = iter.emission(_slotID);
    double diff = x - _mean;
    
    return _A - (diff * diff) / (2 * _var);
  }
  
  virtual void updateParams(EMSequences * sequences, std::vector<EmissionFunction*> * group) {
    if (_is_fixed_mean && _is_fixed_var)
      return;
    
    // suff stats
    double sum_Pzi = 0;
    double sum_Pzi_xi = 0;
    double mu = _mean;
    double sig2 = _var;
    
    std::vector<EmissionFunction*>::iterator ef_it;
    
    if (!_is_fixed_mean) {
      for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
        Normal * ef = (Normal*) (*ef_it)->inner();
        PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
        
        do {
          const double * post_j = post_it->posterior();
          Iter & iter = post_it->iter();
          iter.resetFirst();
          
          for (int j = 0; j < iter.length(); iter.next(), ++j) {
            double x = iter.emission(ef->_slotID);
            
            sum_Pzi += post_j[j];
            sum_Pzi_xi += post_j[j] * x;
          }
        } while (post_it->next());
        
        delete post_it;
      }
      
      mu = sum_Pzi_xi /  sum_Pzi;
    }
    
    // for variance
    if (!_is_fixed_var) {
      double sum_Pzi_sdiff = 0;
      
      if (_is_fixed_mean) {
        for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
          Normal * ef = (Normal*) (*ef_it)->inner();
          PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
          
          do {
            const double * post_j = post_it->posterior();
            Iter & iter = post_it->iter();
            iter.resetFirst();
            
            for (int j = 0; j < iter.length(); iter.next(), ++j) {
              double x = iter.emission(ef->_slotID);
              
              sum_Pzi += post_j[j];
              sum_Pzi_sdiff += post_j[j] * (x - mu) * (x - mu);
            }
          } while (post_it->next());
          
          delete post_it;
        }
      } else {
        for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
          Normal * ef = (Normal*) (*ef_it)->inner();
          PosteriorIterator * post_it = sequences->iterator(ef->_stateID, ef->_slotID);
          
          do {
            const double * post_j = post_it->posterior();
            Iter & iter = post_it->iter();
            iter.resetFirst();
            
            for (int j = 0; j < iter.length(); iter.next(), ++j) {
              double x = iter.emission(ef->_slotID);
              
              sum_Pzi_sdiff += post_j[j] * (x - mu) * (x - mu);
            }
          } while (post_it->next());
          
          delete post_it;
        }
      }
      
      sig2 = sum_Pzi_sdiff / sum_Pzi;
    }
    
    // update values
    _mean = mu;
    _var = sig2;
    update_consts();
    
    // propagate to other elements in the group
    for (ef_it = group->begin(); ef_it != group->end(); ++ef_it) {
      Normal * ef = (Normal*) (*ef_it)->inner();
      
      if (ef != this) {
        ef->_mean = _mean;
        ef->_var = _var;
        ef->update_consts();
      }
    }
  }

private:
  double _mean;
  double _var;
  bool _is_fixed_mean;
  bool _is_fixed_var;
  
  double _A; // -log(sqrt(2 pi var))
  
  void update_consts() {
    _A = -log(sqrt(2.0 * M_PI * _var));
  }
  
};

#endif
