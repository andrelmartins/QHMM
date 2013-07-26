#ifndef INNER_TMPL_HPP
#define INNER_TMPL_HPP

#include "iter.hpp"
#include "logsum.hpp"

//
// Forward Inner Loop
//

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

//
// Backward Inner Loop
//

template<typename FuncAkl, typename FuncEkb>
class InnerBckDense {
public:
  double operator() (const int & n_states, double const * const m_col_next, int k, Iter const & iter, FuncAkl logAkl, FuncEkb logEkb, LogSum * lg) {
    lg->clear();

    for (int l = 0; l < n_states; ++l) {
      double value = m_col_next[l] + (*logAkl)(iter, k, l) + (*logEkb)(iter, l);
      lg->store(value);
    }
    return lg->compute();
  }
};

template<typename FuncAkl, typename FuncEkb>
class InnerBckSparse {
public:
  InnerBckSparse(FuncAkl transitions) : _next(transitions->nextState()), _n_states(transitions->n_states()) {}
  ~InnerBckSparse() { // TODO: clean up memory management responsibilities
    for (int i = 0; i < _n_states; ++i)
      delete[] _next[i];
    delete[] _next;
  }

  double operator() (const int & n_states, double const * const m_col_next, int k, Iter const & iter, FuncAkl logAkl, FuncEkb logEkb, LogSum * lg) {
    lg->clear();

    for (int * ptr = _next[k]; *ptr >= 0; ++ptr) {
      int l = *ptr;
      double value = m_col_next[l] + (*logAkl)(iter, k, l) + (*logEkb)(iter, l);
      lg->store(value);
    }
    return lg->compute();
  }

private:
  const int ** _next;
  const int _n_states;
};


#endif
