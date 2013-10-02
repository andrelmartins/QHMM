#ifndef QHMMTHREADEXCEPTION_HPP
#define QHMMTHREADEXCEPTION_HPP

#include "QHMMException.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class QHMMThreadException {
  QHMMException copy;
  bool has_exception;

public:
  QHMMThreadException() { 
    has_exception = false; 
  }

  ~QHMMThreadException(){
    this->rethrow(); 
  }

  void rethrow(){
    if(has_exception)
      throw QHMMException(copy);
  }

  void captureException(const QHMMException & other) {
    #pragma omp critical
    {
      has_exception = true;
      copy = other;
    }
  }
};

#endif
