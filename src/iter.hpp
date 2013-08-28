#ifndef ITER_HPP
#define ITER_HPP

#include <cstdlib> // for NULL
#include <cassert>
#include <vector>

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
      _index = 0;
    }
    
    void resetLast() {
      _emission_ptr = _emission_end;
      _covar_ptr = _covar_end;
      _missing_ptr = _missing_end;
      _index = _length - 1;
    }
    
    bool next() {
      if (_emission_ptr == _emission_end)
        return false;
      _emission_ptr += _emission_step;
      _covar_ptr += _covar_step;
      _missing_ptr += _missing_step;
      ++_index;
      return true;
    }
    
    bool prev() {
      if (_emission_ptr == _emission_start)
        return false;
      _emission_ptr -= _emission_step;
      _covar_ptr -= _covar_step;
      _missing_ptr -= _missing_step;
      --_index;
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

    int index() const { return _index; }
  
    bool is_missing(int slot) const {
      if (_missing_ptr == NULL)
        return false;
      return (_missing_ptr[slot] != 0);
    }

    bool has_missing() const { return _missing_ptr != NULL; }

    // Partition sequence into blocks with no missing data (for a given slot)
    // NOTES: 1. These wull share data with the parent iterator and, as such, are not valid
    //           after the parent iterator is destroyed.
    //        2. Sub-iterators always report no missing data. That is only valid from the
    //           perspective of the indicated slot (if missing data does not align across
    //           emission slots, then a separate set of sub-iterators must be created per
    //           emission slot).
    std::vector<Iter> * sub_iterators(int slot);

    int emission_slot_count() { return _emission_slot_count; }
    int iter_offset() const { return _offset; }
  
  protected:
    Iter(Iter * parent, int start, int end); // constructor for sub_iterator() function
  
    bool _is_subiterator;
    int _offset;
    int _length;
    int _index;
  
    int _emission_slot_count;
    double * _emission_ptr;
    double * _emission_start;
    double * _emission_end;
    int _emission_step;
    int * _emission_offsets;
  
    int _covar_slot_count;
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
