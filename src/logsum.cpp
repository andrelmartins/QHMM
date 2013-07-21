#include "logsum.hpp"
#include <algorithm>
#include <functional>
#include <cmath>
const double LogSum::SUM_LOG_THRESHOLD = -10;

class LogSum2 : public LogSum {
  public:
    LogSum2(const int capacity, const bool optimize) : LogSum(capacity, optimize) {}
    
    virtual double compute() {
      assert(_count == 2);
      
      double a = _values[0];
      double b = _values[1];
      
      if (a == -std::numeric_limits<double>::infinity())
        return b;
      
      if (b == -std::numeric_limits<double>::infinity())
        return a;
      
      if (a > b) {
        if (b - a < _threshold)
          return a;
        else
          return a + log(exp(b - a));
      } else {
        if (a - b < _threshold)
          return b;
        else
          return b + log(exp(a - b));
      }
    }
};

LogSum * LogSum::create(const int capacity, const bool optimize) {
  // special case implementations
  if (capacity == 2)
    return new LogSum2(capacity, optimize);

  return new LogSum(capacity, optimize);
}

LogSum * LogSum::createType(const int type, const int capacity, const bool optimize) {
  if (capacity == 2 && type == 1)
    return new LogSum2(capacity, optimize);

  return new LogSum(capacity, optimize);
}
   
double LogSum::compute() {
  assert(_count > 0);
    
  double max = *(std::max_element(_values, _values + _count));
  int i;
  double * dptr;
  double expsum = 1;

  if (max == -std::numeric_limits<double>::infinity())
    return max;

  for (i = 0, dptr = _values; i < _count; ++i, ++dptr) {
    double logdiff = (*dptr - max);
    if (logdiff > _threshold)
      expsum += exp(logdiff);
  }

  return max + log(expsum);
}      
