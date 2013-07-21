#ifndef LOGSUM_H
#define LOGSUM_H

#include <cassert>

class LogSum {
  public:
    static const double SUM_LOG_THRESHOLD;

    static LogSum * create(const int capacity);
    // for benchmarks
    static LogSum * createType(const int type, const int capacity);

    virtual ~LogSum() { delete _values; }
    
    void store(const double logValue) {
      assert(_count < _length);
      _values[_count++] = logValue;
    }
    
    void clear() {
      _count = 0;
    }
    
    virtual double compute();
    
  protected:
    LogSum(const int capacity) : _values(new double[capacity]), _length(capacity), _count(0) {}
    
    double * _values;
    int _length;
    int _count;
};

#endif
