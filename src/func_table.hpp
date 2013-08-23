#ifndef FUNC_TABLE_HPP
#define FUNC_TABLE_HPP

#include <algorithm>
#include <stdexcept>
#include <limits>
#include <vector>
#include <cassert>
#include <cstring>
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
  
    std::vector<T *> _singletons;
    std::vector<std::vector<T *> > _groups;
  
  public:
    void insert(T * func) {
      assert((int) _funcs.size() < _n_states);
      _funcs.push_back(func);
      _singletons.push_back(func);
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
  
    virtual void makeGroup(int * idxs, int length, int slot = 0) {
      std::vector<T * > group;
      
      for (int i = 0; i < length; ++i) {
        T * f_i = _funcs[idxs[i]];
        group.push_back(f_i);
        
        // remove from singletons
        typename std::vector<T*>::iterator it = find(_singletons.begin(), _singletons.end(), f_i);
        _singletons.erase(it);
      }
      
      _groups.push_back(group);
    }
  
    // turn remaining singletons to groups
    virtual void commitGroups() {
      typename std::vector<T*>::iterator it;
      
      for (it = _singletons.begin(); it != _singletons.end(); ++it) {
        std::vector<T*> group;
        group.push_back(*it);
        _groups.push_back(group);
      }
    }
  
    virtual const std::vector<std::vector<T*> > & groups() {
      return _groups;
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
      TransitionFunction * f_i = _funcs[i];
      int length = f_i->n_targets() + 1;
      next[i] = new int[length];
      next[i][f_i->n_targets()] = -1; // termination mark

      memcpy(next[i], f_i->targets(), f_i->n_targets() * sizeof(int));
      // sort to keep memory access sequential
      std::sort(next[i], next[i] + length);
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
  
  Params * getSlotParams(int state, int slot) const {
    return getParams(state);
  }

  void setSlotParams(int state, int slot, Params const & params) {
    setParams(state, params);
  }
  
  double operator() (Iter const & iter, int i) const {
    return _funcs[i]->log_probability(iter);
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
  
  Params * getSlotParams(int state, int slot) const {
    return _funcs[state][slot]->getParams();
  }

  void setSlotParams(int state, int slot, Params const & params) {
    _funcs[state][slot]->setParams(params);
  }
  
  void insert(std::vector<EmissionFunction *> funcs) {
    assert((int)_funcs.size() < _n_states);
    assert((int) funcs.size() == _n_slots);
    _funcs.push_back(funcs);
    
    std::vector<EmissionFunction*>::iterator it;
    int slot = 0;
    for (it = funcs.begin(); it != funcs.end(); ++it, ++slot) {
      _singletons[slot].push_back(*it);
    }
  }
  
  double operator() (Iter const & iter, int i) const {
    double log_prob = 0;
    
    for (int slot = 0; slot < _n_slots; ++slot)
      log_prob += _funcs[i][slot]->log_probability(iter);
    
    return log_prob;
  }
  
  int n_states() { return _n_states; }
  int n_slots() { return _n_slots; }
  
  virtual void makeGroup(int * idxs, int length, int slot) {
    std::vector<EmissionFunction * > group;
    
    for (int i = 0; i < length; ++i) {
      EmissionFunction * f_i = _funcs[idxs[i]][slot];
      group.push_back(f_i);
      
      // remove from singletons
      std::vector<EmissionFunction*>::iterator it = find(_singletons[slot].begin(), _singletons[slot].end(), f_i);
      _singletons[slot].erase(it);
    }
    
    _groups.push_back(group);
  }
  
  // turn remaining singletons to groups
  virtual void commitGroups() {
    std::vector<EmissionFunction*>::iterator it;
    
    for (int slot = 0; slot < _n_slots; ++slot) {
      for (it = _singletons[slot].begin(); it != _singletons[slot].end(); ++it) {
        std::vector<EmissionFunction*> group;
        group.push_back(*it);
        _groups.push_back(group);
      }
    }
  }
  
  virtual const std::vector<std::vector<EmissionFunction*> > & groups() {
    return _groups;
  }
  
  private:
    const int _n_states;
    const int _n_slots;
    std::vector<std::vector<EmissionFunction *> > _funcs;
  
    std::vector<std::vector<EmissionFunction *> > _singletons;
    std::vector<std::vector<EmissionFunction *> > _groups;

};

#endif
