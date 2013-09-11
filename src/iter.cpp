#include "iter.hpp"
#include <cstdio>

Iter::Iter(int length, int emission_slots, int * e_slot_dim, double * emissions,
         int covar_slots, int * c_slot_dim, double * covars, int * missing) {
  assert(length > 0);
  assert(emission_slots > 0);
  assert(e_slot_dim != NULL);
  assert(emissions != NULL);
  
  _is_subiterator = false;
  _offset = 0;
  _length = length;
  _index = 0;
  _emission_slot_count = emission_slots;
  _covar_slot_count = covar_slots;
  
  _emission_step = 0;
  _emission_offsets = new int[emission_slots];
  for (int i = 0; i < emission_slots; ++i) {
    _emission_offsets[i] = _emission_step;
    _emission_step += e_slot_dim[i];
  }
  _emission_ptr = emissions;
  _emission_start = emissions;
  _emission_end = emissions + (length - 1) * _emission_step;
  
  _covar_step = 0;
  if (covar_slots == 0)
    _covar_offsets = NULL;
  else
    _covar_offsets = new int[covar_slots];
  for (int i = 0; i < covar_slots; ++i) {
    _covar_offsets[i] = _covar_step;
    _covar_step += c_slot_dim[i];
  }
  _covar_ptr = covars;
  _covar_start = covars;
  _covar_end = covars + (length - 1) * _covar_step;

  _missing_step = (missing == NULL ? 0 : emission_slots);
  _missing_ptr = missing;
  _missing_start = missing;
  _missing_end = missing + (length - 1) * _missing_step;
}

Iter::~Iter() {
  if (!_is_subiterator) {
    delete[] _emission_offsets;
    if (_covar_offsets != NULL)
      delete[] _covar_offsets;
  }
}

Iter::Iter(Iter * parent, int start, int end) : _is_subiterator(true), _missing_ptr(NULL), _missing_start(NULL), _missing_end(NULL), _missing_step(0) {
  // set length
  _length = end - start + 1;
  _index = 0;
  _offset = start;
  
  // copy shared information
  _emission_slot_count = parent->_emission_slot_count;
  _covar_slot_count = parent->_covar_slot_count;
  _emission_offsets = parent->_emission_offsets;
  _covar_offsets = parent->_covar_offsets;
  _emission_step = parent->_emission_step;
  _covar_step = parent->_covar_step;

  // start/end/ptr pointers
  _emission_start = parent->_emission_start + parent->_emission_step * start;
  _emission_end = parent->_emission_start + parent->_emission_step * end;
  _emission_ptr = _emission_start;
  
  _covar_start = parent->_covar_start + parent->_covar_step * start;
  _covar_end = parent->_covar_start + parent->_covar_step * end;
  _covar_ptr = _covar_start;
}

std::vector<Iter> * Iter::sub_iterators(int slot) {
  assert(!_is_subiterator);
  assert(slot >= 0 && slot < _emission_slot_count);
  std::vector<Iter> * result = new std::vector<Iter>();
  
  if (has_missing()) {
    int start = -1;
    int * mptr = _missing_start + slot;
    
    for (int i = 0; i < _length; ++i, mptr += _missing_step) {
      if (*mptr != 0 && start >= 0) {
        result->push_back(Iter(this, start, i - 1));
        start = -1;
      } else if (start < 0)
        start = i;
    }
    
    // close last, if any
    if (start >= 0)
      result->push_back(Iter(this, start, _length - 1));
  } else
    result->push_back(Iter(this, 0, _length - 1));

  return result;
}
