#include "iter.hpp"

Iter::Iter(int length, int emission_slots, int * e_slot_dim, double * emissions,
         int covar_slots, int * c_slot_dim, double * covars, int * missing) {
  assert(length > 0);
  assert(emission_slots > 0);
  assert(e_slot_dim != NULL);
  assert(emissions != NULL);
  
  _length = length;
  
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
  delete[] _emission_offsets;
  delete[] _covar_offsets;
}
