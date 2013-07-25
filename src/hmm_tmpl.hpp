#ifndef HMMTMPL_HPP
#define HMMTMPL_HPP

#include "hmm.hpp"
#include "logsum.hpp"
#include "inner_tmpl.hpp"
#include "func_table.hpp"

template <typename InnerFwd, typename FuncAkl, typename FuncEkb>
class HMMImpl : public HMM {
  private:
    const int _n_states;
    const FuncAkl _logAkl;
    const FuncEkb _logEkb;
    const InnerFwd _innerFwd;
    
    double * _init_log_probs;
    
  public:
    HMMImpl(InnerFwd innerFwd, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) : _n_states(logAkl->n_states()), _logAkl(logAkl), _logEkb(logEkb), _innerFwd(innerFwd), _init_log_probs(init_log_probs) { }

    ~HMMImpl() {
      delete _innerFwd; // TODO: clean up memory management responsibilities
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

};

// auxiliary function to enable type inference
template <typename InnerFwd, typename FuncAkl, typename FuncEkb>
HMM * new_hmm_instante(InnerFwd innerFwd, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) {
  return new HMMImpl<InnerFwd, FuncAkl, FuncEkb>(innerFwd, logAkl, logEkb, init_log_probs);
}

#endif
