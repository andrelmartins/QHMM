#ifndef QHMMTHREADHELPER_HPP
#define QHMMTHREADHELPER_HPP

#include "QHMMException.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class QHMMThreadHelper {
  QHMMException copy;
  bool has_exception;

public:
  QHMMThreadHelper() { 
    has_exception = false; 
  }

  ~QHMMThreadHelper(){
    this->rethrow(); 
  }

  void rethrow(){
    if(has_exception)
      throw QHMMException(copy);
  }

  void captureException(const QHMMException & other) {
    #pragma omp critical
    {
      if (!has_exception) { // capture only first exception
	has_exception = true;
	copy = other;
      }
    }
  }
};

#endif
