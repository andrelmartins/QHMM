#ifndef INNER_TMPL_HPP
#define INNER_TMPL_HPP

#include "iter.hpp"
#include "logsum.hpp"

template<typename FuncType>
class InnerFwdDense {
public:
  double operator() (const int & n_states, double const * const m_col_prev, int l, Iter const & iter, FuncType logAkl, LogSum * lg) {
    lg->clear();
    
    for (int k = 0; k < n_states; ++k)
      lg->store(m_col_prev[k] + (*logAkl)(iter, k, l));

    return lg->compute();
  }
};

template<typename FuncType>
class InnerFwdSparse {
public:
  InnerFwdSparse(FuncType transitions) : _previous(transitions->previousState()), _n_states(transitions->n_states()) {}
  ~InnerFwdSparse() { // TODO: clean up memory management responsibilities
    for (int i = 0; i < _n_states; ++i)
      delete[] _previous[i];
    delete[] _previous;
  }

  double operator() (const int & n_states, double const * const  m_col_prev, int l, Iter const & iter, FuncType logAkl, LogSum * lg) {
    lg->clear();
  
    for (int * ptr = _previous[l]; *ptr >= 0; ++ptr)
      lg->store(m_col_prev[*ptr] + (*logAkl)(iter, *ptr, l));

    return lg->compute();
  }

private:
  const int ** _previous;
  const int _n_states;
};


#endif
