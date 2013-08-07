#ifndef HMMTMPL_HPP
#define HMMTMPL_HPP

#include "logsum.hpp"
#include "inner_tmpl.hpp"
#include "func_table.hpp"
#include "hmm.hpp"

#include <cmath>

template <typename InnerFwd, typename InnerBck, typename FuncAkl, typename FuncEkb>
class HMMImpl : public HMM {
  private:
    const int _n_states;
    const FuncAkl _logAkl;
    const FuncEkb _logEkb;
    const InnerFwd _innerFwd;
    const InnerBck _innerBck;
    
    double * _init_log_probs;
    
  public:
    HMMImpl(InnerFwd innerFwd, InnerBck innerBck, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) : _n_states(logAkl->n_states()), _logAkl(logAkl), _logEkb(logEkb), _innerFwd(innerFwd), _innerBck(innerBck), _init_log_probs(init_log_probs) { }

    ~HMMImpl() {
      delete _innerFwd; // TODO: clean up memory management responsibilities
      delete _innerBck; // ""
    }
  
    virtual bool valid_transition_params(int state, Params const & params) const {
      return _logAkl->validParams(state, params);
    }
  
    virtual bool valid_emission_params(int state, int slot, Params const & params) const {
      return _logEkb->validSlotParams(state, slot, params);
    }

    virtual void set_transition_params(int state, Params const & params) {
      _logAkl->setParams(state, params);
    }
  
    virtual void set_emission_params(int state, int slot, Params const & params) {
      _logEkb->setSlotParams(state, slot, params);
    }

    virtual void set_initial_probs(double * probs) {
      for (int i = 0; i < _n_states; ++i)
        _init_log_probs[i] = log(probs[i]);
    }
  
    double forward(Iter & iter, double * matrix) {
      double * m_col, * m_col_prev;
      LogSum * logsum = LogSum::create(_n_states);
      iter.resetFirst();
    
      /* border conditions - position i = 0
       * f_k(0) = e_k(0) * a0k
       * log f_k(0) = log e_k(0) + log a0k
       */
      m_col = matrix;
      for (int k = 0; k < _n_states; ++k)
        m_col[k] = (*_logEkb)(iter, k) + _init_log_probs[k];

      /* inner cells */
      m_col_prev = m_col;
      m_col += _n_states;
      for (; iter.next(); m_col_prev += _n_states, m_col += _n_states) {
        for (int l = 0; l < _n_states; ++l)
          m_col[l] = (*_logEkb)(iter, l) + (*_innerFwd)(_n_states, m_col_prev, l, iter, _logAkl, logsum);
      }
      
      /* log-likelihood */
      logsum->clear();
      m_col = matrix + (iter.length() - 1)*_n_states;
      for (int i = 0; i < _n_states; ++i)
        logsum->store(m_col[i]);
      double loglik = logsum->compute();

      // clean up
      delete logsum;
      
      return loglik;
    }

    double backward(Iter & iter, double * matrix) {
      double * m_col, * m_col_next;
      LogSum * logsum = LogSum::create(_n_states);
      
      /* border conditions @ position = N - 1*/
      m_col = matrix + (iter.length() - 1)*_n_states;
      for (int k = 0; k < _n_states; ++k)
        m_col[k] = 0; /* log(1) */

      /* inner cells */
      iter.resetLast();
      m_col_next = m_col;
      m_col -= _n_states;
      for (; m_col >= matrix; m_col_next -= _n_states, m_col -= _n_states, iter.prev()) {

        for (int k = 0; k < _n_states; ++k)
          m_col[k] = (*_innerBck)(_n_states, m_col_next, k, iter, _logAkl, _logEkb, logsum);
      }

      /* log-likelihood */
      m_col = matrix;
      logsum->clear();
      iter.resetFirst();
      for (int k = 0; k < _n_states; ++k) {
        double value = m_col[k] + _init_log_probs[k] + (*_logEkb)(iter, k);
        logsum->store(value);
      }
      double loglik = logsum->compute();

      // clean-up
      delete logsum;
      
      return loglik;
    }

#define AT(M, I, J) M[(I) + (J)*rows]

    void viterbi(Iter & iter, int * path) {
      int rows = _n_states; /* needed by AT macro */
      int cols = iter.length();
      double * m_col, * m_col_prev;
      int * b_col;
      int * backptr;
      int * pptr;
      double * matrix;

      /* setup matrices */
      matrix = new double[rows*cols];
      backptr = new int[rows*cols];

      /* fill first column */
      for (int l = 0; l < _n_states; ++l) {
        AT(matrix, l, 0) = (*_logEkb)(iter, l) + _init_log_probs[l];
        AT(backptr, l, 0) = -1; /* stop */
      }

      /* inner columns */
      m_col_prev = matrix;
      m_col = matrix + rows;
      b_col = backptr + rows;
      iter.resetFirst();
      for ( ; iter.next(); m_col += rows, m_col_prev += rows, b_col += rows) {

        for (int l = 0; l < _n_states; ++l) {
          double max = -std::numeric_limits<double>::infinity();
          int argmax = -1;

          for (int k = 0; k < _n_states; ++k) {
            double value = m_col_prev[k] + (*_logAkl)(iter, k, l);
            
            if (value > max) {
              max = value;
              argmax = k;
            }
          }

          /* assert(argmax != -1); */
          m_col[l] = (*_logEkb)(iter, l) + max;
          b_col[l] = argmax;
        }
      }

      /* backtrace */

      /* last state */
      iter.resetLast();
      pptr = path + iter.length() - 1;
      {
        double max = -std::numeric_limits<double>::infinity();
        int argmax = -1;
        m_col = matrix + (iter.length() - 1)*rows;

        for (int k = 0; k < _n_states; ++k) { // TODO: Consider inner loop sparse optimization
          double value = m_col[k];
          if (value > max) {
            max = value;
            argmax = k;
          }
        }
        /* assert(argmax != -1); */
        *pptr = argmax;
      }

      /* other states */
      int z = *pptr;
      --pptr;
      for (int l = iter.length() - 1; iter.prev() ; --pptr) {
        z = AT(backptr, z, l);
        *pptr = z;
        /* assert(prev >= 0); */
      }
      
      // clean up
      delete[] matrix;
      delete[] backptr;
    }
};

// auxiliary function to enable type inference
template <typename InnerFwd, typename InnerBck, typename FuncAkl, typename FuncEkb>
HMM * new_hmm_instance(InnerFwd innerFwd, InnerBck innerBck, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) {
  return new HMMImpl<InnerFwd, InnerBck, FuncAkl, FuncEkb>(innerFwd, innerBck, logAkl, logEkb, init_log_probs);
}

#endif
