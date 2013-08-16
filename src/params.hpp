#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <stdexcept>
#define INDEX_CHECK(index) if (index < 0 || index >= _length) throw std::out_of_range("index");

class Params {

public:
  
  Params(unsigned int length, const double * const values) : _length(length) {
    _allFixed = false;
    _fixed = new bool[length];
    _values = new double[length];
    for (unsigned int i = 0; i < length; ++i) {
      _fixed[i] = false;
      _values[i] = values[i];
    }
  }
  
  ~Params() {
    delete[] _values;
    delete[] _fixed;
  }
  
  bool isAllFixed() const {
    return _allFixed;
  }
  
  bool isFixed(unsigned int index) const {
    INDEX_CHECK(index);
    return _fixed[index];
  }
  
  void setFixed(unsigned int index, bool value) {
    INDEX_CHECK(index);
    _fixed[index] = value;
    
    if (_allFixed)
      _allFixed = value;
    else {
      _allFixed = value;
      for (unsigned int i = 0; i < _length && _allFixed; ++i)
        _allFixed = _allFixed && _fixed[i];
    }
  }
  
  double & operator[](unsigned int index) {
    INDEX_CHECK(index);
    return _values[index];
  }

  const double & operator[](unsigned int index) const {
    INDEX_CHECK(index);
    return _values[index];
  }
  
  int length() const {
    return _length;
  }
  
private:
  const unsigned int _length;
  bool * _fixed;
  bool _allFixed;
  double * _values;
  
};

#endif
