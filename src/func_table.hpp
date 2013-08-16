#ifndef FUNC_TABLE_HPP
#define FUNC_TABLE_HPP

#include <stdexcept>
#include <limits>
#include <vector>
#include <cassert>
#include "base_classes.hpp"

template<typename T>
class FunctionTable {
  protected:
    FunctionTable(int n_states) : _n_states(n_states) {
      _funcs.reserve(n_states);
    }

    ~FunctionTable() { // TODO: review memory management responsabilities
      for (unsigned int i = 0; i < _funcs.size(); ++i)
        delete _funcs[i];
    }
    
    const int _n_states;
    std::vector<T *> _funcs;
    
  public:
    void insert(T * func) {
      assert((int) _funcs.size() < _n_states);
      _funcs.push_back(func);
    }
  
    bool validParams(int state, Params const & params) const {
      return _funcs[state]->validParams(params);
    }

    virtual Params * getParams(int state) const {
      return _funcs[state]->getParams();
    }

    virtual void setParams(int state, Params const & params) {
      _funcs[state]->setParams(params);
    }
  
    int n_states() const { return _n_states; }
};

class HomogeneousTransitions : public FunctionTable<TransitionFunction> {
	public:
		HomogeneousTransitions(int n_states) : FunctionTable<TransitionFunction>(n_states) {
		  _m = new double*[n_states];
		  for (int i = 0; i < n_states; ++i)
		    _m[i] = new double[n_states];
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
  
    virtual void setParams(int state, Params const & params) {
      FunctionTable<TransitionFunction>::setParams(state, params);
      
      /* partial update */
      double * row = _m[state];
      for (int j = 0; j < _n_states; ++j)
        row[j] = _funcs[state]->log_probability(j);
    }
		
		bool isSparse() {
		  int invalid_count = 0;
		  for (int i = 0; i < _n_states; ++i)
		    for (int j = 0; j < _n_states; ++j)
		      if (_m[i][j] == -std::numeric_limits<double>::infinity())
		        ++invalid_count;
		  return invalid_count >= (_n_states * _n_states / 2); // heuristic threshold
		}
		
		int ** previousStates() {
		  int ** previous = new int*[_n_states];
		  
		  for (int j = 0; j < _n_states; ++j) {
		    int length = 0;
		    
		    // count valid transitions
		    for (int i = 0; i < _n_states; ++i)
		      if (_m[i][j] != -std::numeric_limits<double>::infinity())
		        ++length;
		    
		    // record valid transitions
		    previous[j] = new int[length + 1];
		    previous[j][length] = -1; // termination mark
		    int k = 0;
		    for (int i = 0; i < _n_states; ++i)
		      if (_m[i][j] != -std::numeric_limits<double>::infinity())
		        previous[j][k++] = i;
		  }
		  
		  return previous;
		}
		
		int ** nextStates() {
		  int ** next = new int*[_n_states];
		  
		  for (int i = 0; i < _n_states; ++i) {
		    int length = 0;
		    
		    // count valid transitions
		    for (int j = 0; j < _n_states; ++j)
		      if (_m[i][j] != -std::numeric_limits<double>::infinity())
		        ++length;
		    
		    // record valid transitions
		    next[i] = new int[length + 1];
		    next[i][length] = -1; // termination mark
		    int k = 0;
		    for (int j = 0; j < _n_states; ++j)
		      if (_m[i][j] != -std::numeric_limits<double>::infinity())
		        next[i][k++] = i;
		  }
		  
		  return next;
		}
		
		// TODO: Check if it's worth it to transpose this matrix, since we'll be accessing
		//       it column by columns
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
		
		bool isSparse() {
		  return false; // For now assume non-homogeneous means not-sparse
		                // this is not strictly true, we could have constraints on valid
		                // transitions that make things sparse ...
		}
  
    int ** previousStates() {
      throw std::logic_error("called previousStates on non-homogeneous class (not supported)");
    }
  
    int ** nextStates() {
      throw std::logic_error("called nextStates on non-homogeneous class (not supported)");      
    }
};

class Emissions : public FunctionTable<EmissionFunction> {
public:
  Emissions(int n_states) : FunctionTable<EmissionFunction>(n_states) {}

  bool validSlotParams(int state, int slot, Params const & params) const {
    return validParams(state, params);
  }
  
  void setSlotParams(int state, int slot, Params const & params) {
    setParams(state, params);
  }
  
  double operator() (Iter const & iter, int i) const {
    return _funcs[i]->log_probability(iter, 0);
  }
};

class MultiEmissions {
public:
  MultiEmissions(int n_states, int n_slots) : _n_states(n_states), _n_slots(n_slots) {
    _funcs.reserve(n_states);
  }
  
  ~MultiEmissions() { // TODO: review memory management responsabilities
    for (unsigned int i = 0; i < _funcs.size(); ++i) {
      std::vector<EmissionFunction *> funcs_i = _funcs[i];
      
      for (unsigned int j = 0; j < funcs_i.size(); ++j)
        delete funcs_i[j];
    }
  }
  
  bool validSlotParams(int state, int slot, Params const & params) const {
    return _funcs[state][slot]->validParams(params);
  }
  
  void setSlotParams(int state, int slot, Params const & params) {
    _funcs[state][slot]->setParams(params);
  }
  
  void insert(std::vector<EmissionFunction *> funcs) {
    assert((int)_funcs.size() < _n_states);
    assert((int) funcs.size() == _n_slots);
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
