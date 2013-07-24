#ifndef LOGSUM_HPP
#define LOGSUM_HPP

#include <cassert>
#include <limits>

class LogSum {
  public:
    static const double SUM_LOG_THRESHOLD;

    static LogSum * create(const int capacity, const bool optimize = true);
    // for benchmarks
    static LogSum * createType(const int type, const int capacity, const bool optimize);
    
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
    LogSum(const int capacity, const bool optimize) : _threshold(optimize ? SUM_LOG_THRESHOLD : -std::numeric_limits<double>::infinity()), _values(new double[capacity]), _length(capacity), _count(0) {}
    
    const double _threshold;
    double * _values;
    int _length;
    int _count;
};

#endif
