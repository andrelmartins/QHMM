#ifndef BASE_FUNC_TABLE_HPP
#define BASE_FUNC_TABLE_HPP

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
  virtual void insert(T * func) {
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
  
  virtual bool setCovars(int state, int * idxs, int length) {
    return _funcs[state]->setCovarSlots(idxs, length);
  }
  
  virtual bool getOption(int state, const char * name, double * out_value) const {
    return _funcs[state]->getOption(name, out_value);
  }
  
  virtual bool setOption(int state, const char * name, double value) {
    return _funcs[state]->setOption(name, value);
  }
  
  virtual void makeGroup(int * idxs, int length) {
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
  
  const T * function(int index) {
    assert(index >= 0 && index < _n_states);
    return _funcs[index];
  }
};

class TransitionTable : public FunctionTable<TransitionFunction> {
public:
  TransitionTable(int n_states) : FunctionTable<TransitionFunction>(n_states) {}
  
  bool isSparse() {
    int invalid_count = 0;
    for (int i = 0; i < _n_states; ++i)
      invalid_count += (_n_states - _funcs[i]->n_targets());
    return invalid_count >= (_n_states * _n_states / 2); // heuristic threshold
  }
  
  int ** previousStates() {
    int ** previous = new int*[_n_states];
    
    // count valid transitions
    int lengths[_n_states];
    
    for (int i = 0; i < _n_states; ++i)
      lengths[i] = 0;
    
    for (int i = 0; i < _n_states; ++i) {
      TransitionFunction * f_i = _funcs[i];
      int n = f_i->n_targets();
      const int * targets = f_i->targets();
      
      for (int j = 0; j < n; ++j)
        lengths[targets[j]]++;
    }
    
    // allocate space
    for (int i = 0; i < _n_states; ++i) {
      previous[i] = new int[lengths[i] + 1];
      previous[i][lengths[i]] = - 1; // termination mark
    }
    
    // record valid transitions
    for (int i = 0; i < _n_states; ++i) {
      TransitionFunction * f_i = _funcs[i];
      int n = f_i->n_targets();
      const int * targets = f_i->targets();
      
      for (int j = 0; j < n; ++j) {
        int k = targets[j];
        --lengths[k];
        previous[k][lengths[k]] = i;
      }
    }
    
    // sort transitions
    for (int i = 0; i < _n_states; ++i)
      std::sort(previous[i], previous[i] + lengths[i]);
    
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
      std::sort(next[i], next[i] + length - 1); // exclude termination mark
    }
    
    return next;
  }
  
  virtual void refresh() {}
};

class EmissionTable {

public:
  virtual bool validSlotParams(int state, int slot, Params const & params) const  = 0;
  
  virtual Params * getSlotParams(int state, int slot) const = 0;
  
  virtual void setSlotParams(int state, int slot, Params const & params) = 0;
  
  virtual bool setSlotCovars(int state, int slot, int * idxs, int length) = 0;
  
  virtual bool getSlotOption(int state, int slot, const char * name, double * out_value) const = 0;
  
  virtual bool setSlotOption(int state, int slot, const char * name, double value) = 0;
  
  // groups
  virtual void insert(std::vector<EmissionFunction *> funcs) = 0;
  
  virtual void makeGroup(int * idxs, int length, int * slots = NULL, int n_slots = 0) = 0;
  
  virtual void makeGroupExt(int length, int * idxs, int * slots = NULL) = 0;
  
  virtual void commitGroups() = 0;
  
  virtual const std::vector<std::vector<EmissionFunction*> > & groups() = 0;
  
  virtual int n_states() const = 0;
};

#endif
