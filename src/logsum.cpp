#include "logsum.hpp"
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>
const double LogSum::SUM_LOG_THRESHOLD = -10;

class LogSum2 : public LogSum {
  public:
    LogSum2(const int capacity) : LogSum(capacity) {}
    
    virtual double compute() {
      assert(_count == 2);
      
      double a = _values[0];
      double b = _values[1];
      
      if (a == -std::numeric_limits<double>::infinity())
        return b;
      
      if (b == -std::numeric_limits<double>::infinity())
        return a;
      
      if (a > b) {
        if (b - a < SUM_LOG_THRESHOLD)
          return a;
        else
          return a + log(exp(b - a));
      } else {
        if (a - b < SUM_LOG_THRESHOLD)
          return b;
        else
          return b + log(exp(a - b));
      }
    }
};

class LogSumSmall : public LogSum {
  public:
    LogSumSmall(const int capacity) : LogSum(capacity) {}
    
    virtual double compute() {
      assert(_count > 0);
    
      double max = *(std::max_element(_values, _values + _count));
      int i;
      double * dptr;
      double expsum = 1;
      
      if (max == -std::numeric_limits<double>::infinity())
        return max;
      
      for (i = 0, dptr = _values; i < _count; ++i, ++dptr) {
        double logdiff = (*dptr - max);
        if (logdiff > SUM_LOG_THRESHOLD)
          expsum += exp(logdiff);
      }
      
      return max + log(expsum);
    }
};

LogSum * LogSum::create(const int capacity) {
  // special case implementations
  if (capacity == 2)
    return new LogSum2(capacity);
  if (capacity < 5)
    return new LogSumSmall(capacity);

  // generic implementation
  return new LogSum(capacity);
}

LogSum * LogSum::createType(const int type, const int capacity) {
  switch (type) {
    case 1: return new LogSum2(capacity);
    case 2: return new LogSumSmall(capacity);
    default: return new LogSum(capacity);
  }
}

// Inspired by code in Adam Siepel's Phast Library
double LogSum::compute() {
  assert(_count > 0);

  double max, expsum;

  // sort in descending order
  std::sort(_values, _values + _count, std::greater<double>());

  // sum_i exp(x_i - max(x))
  int i;
  double * dptr;
  double log_diff;
  
  max = _values[0];
  if (max == -std::numeric_limits<double>::infinity())
    return max;
  
  expsum = 1;
  for (i = 1, dptr = _values + 1, log_diff = 0;
       i < _count && (log_diff = *dptr - max) > SUM_LOG_THRESHOLD;
       ++i, ++dptr) {
    expsum += exp(log_diff);
  }
  
  return max + log(expsum);
}
