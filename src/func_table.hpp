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
    return _funcs[i]->log_probability(iter);
  }
};

#endif
