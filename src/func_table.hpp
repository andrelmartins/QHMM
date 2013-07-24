#ifndef FUNC_TABLE_HPP
#define FUNC_TABLE_HPP

#include <vector>
#include <cassert>
#include "base_classes.hpp"

template<typename T>
class FunctionTable {
  protected:
    FunctionTable(int n_states) : _n_states(n_states) {
      _funcs.reserve(n_states);
    }
    
    void insert(T * func) {
      assert(_funcs.size() < _n_states);
      _funcs.push_back(func);
    }
  
    const int _n_states;
    std::vector<T *> _funcs;
    
  public:
    int n_states() const { return _n_states; }
};

class HomogeneousTransitions : public FunctionTable<TransitionFunction> {
	public:
		HomogeneousTransitions(int n_states) : FunctionTable<TransitionFunction>(n_states) {
		  _m = new double*[n_states];
		  for (int i = 0; i < n_states; ++i)
		    _m[i] = new double[n_states];
		  updateTransitions();
		}
		
		~HomogeneousTransitions() {
		  for (int i = 0; i < _n_states; ++i)
		    delete[] _m[i];
		  delete[] _m;
		}
		
		void updateTransitions() {
		  for (int i = 0; i < _n_states; ++i) {
		    double * row = _m[i];
		    for (int j = 0; j < _n_states; ++j)
		      row[j] = _funcs[i]->log_probability(j);
		  }
		}
		
		double operator() (Iter const & iter, int i, int j) const {
			return _m[i][j];
		}
		
	private:
		double ** _m;
};

class NonHomogeneousTransitions : public FunctionTable<TransitionFunction> {
	public:
		NonHomogeneousTransitions(int n_states) : FunctionTable<TransitionFunction>(n_states) {}
		
		double operator() (Iter const & iter, int i, int j) const {
			return _funcs[i]->log_probability(iter, j);
		}
};

class Emissions : public FunctionTable<EmissionFunction> {
  Emissions(int n_states) : FunctionTable<EmissionFunction>(n_states) {}
  
  double operator() (Iter const & iter, int i) const {
    return _funcs[i]->log_probability(iter, 0);
  }
};

class MultiEmissions {
  MultiEmissions(int n_states, int n_slots) : _n_states(n_states), _n_slots(n_slots) {
    _funcs.reserve(n_states);
  }
  
  void insert(std::vector<EmissionFunction *> funcs) {
      assert(_funcs.size() < _n_states);
      assert(funcs.size() == _n_slots);
      _funcs.push_back(funcs);
    }
  
  double operator() (Iter const & iter, int i) const {
    double log_prob = 0;
    
    for (int slot = 0; slot < _n_slots; ++slot)
      log_prob += _funcs[i][slot]->log_probability(iter, slot);
    
    return log_prob;
  }
  
  int n_states() { return _n_states; }
  int n_slots() { return _n_slots; }
  
  private:
    const int _n_states;
    const int _n_slots;
    std::vector<std::vector<EmissionFunction *> > _funcs;
};

#endif
