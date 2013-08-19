#ifndef ITER_HPP
#define ITER_HPP

#include <cstdlib> // for NULL
#include <cassert>

class Iter {

  public:
    // Missing data must follow the same slot count as emissions
    Iter(int length, int emission_slots, int * e_slot_dim, double * emissions,
         int covar_slots, int * c_slot_dim, double * covars, int * missing = NULL);
    virtual ~Iter();
  
    // control ops
    void resetFirst() {
      _emission_ptr = _emission_start;
      _covar_ptr = _covar_start;
      _missing_ptr = _missing_start;
    }
    
    void resetLast() {
      _emission_ptr = _emission_end;
      _covar_ptr = _covar_end;
      _missing_ptr = _missing_end;
    }
    
    bool next() {
      if (_emission_ptr == _emission_end)
        return false;
      _emission_ptr += _emission_step;
      _covar_ptr += _covar_step;
      _missing_ptr += _missing_step;
      return true;
    }
    
    bool prev() {
      if (_emission_ptr == _emission_start)
        return false;
      _emission_ptr -= _emission_step;
      _covar_ptr -= _covar_step;
      _missing_ptr -= _missing_step;
      return true;
    }
        
    // data ops
    double emission(const int slot) const {
      return _emission_ptr[_emission_offsets[slot]];
    }
    
    double emission_i(const int slot, const int i) const {
      return _emission_ptr[_emission_offsets[slot] + i];
    }
    
    double covar(const int slot) const {
      assert(_covar_start != NULL);
      return _covar_ptr[_covar_offsets[slot]];
    }
    
    double covar_i(const int slot, const int i) const {
      assert(_covar_start != NULL);
      return _covar_ptr[_covar_offsets[slot] + i];
    }
    
    // offset with respect to current iterator sequence position
    double covar_ext(const int slot, const int i, const int offset) const {
      assert(_covar_ptr != NULL);
      return _covar_ptr[offset * _covar_step + _covar_offsets[slot]];
    }
    
    int length() const { return _length; }
  
    bool is_missing(int slot) const {
      if (_missing_ptr == NULL)
        return false;
      return (_missing_ptr[slot] != 0);
    }
    
  protected:
    int _length;
  
    double * _emission_ptr;
    double * _emission_start;
    double * _emission_end;
    int _emission_step;
    int * _emission_offsets;
    
    double * _covar_ptr;
    double * _covar_start;
    double * _covar_end;
    int _covar_step;
    int * _covar_offsets;
  
    int * _missing_ptr;
    int * _missing_start;
    int * _missing_end;
    int _missing_step;
  
    // could also store slot descriptions and add asserts to data ops
    // ideally, there should be a non-performance penalty way to add the
    // checks at runtime for debug purposes ... (exceptions ?)
};

#endif
