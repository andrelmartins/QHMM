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
  protected:
  
    virtual const std::vector<std::vector<EmissionFunction*> > & emission_groups() const {
      return _logEkb->groups();
    }

    virtual const std::vector<std::vector<TransitionFunction*> > & transition_groups() const {
      return _logAkl->groups();
    }

    virtual void refresh_transition_table() {
      return _logAkl->refresh();
    }
  
  public:
    HMMImpl(InnerFwd innerFwd, InnerBck innerBck, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) : _n_states(logAkl->n_states()), _logAkl(logAkl), _logEkb(logEkb), _innerFwd(innerFwd), _innerBck(innerBck), _init_log_probs(init_log_probs) { }

    ~HMMImpl() {
      delete _innerFwd; // TODO: clean up memory management responsibilities
      delete _innerBck; // ""
    }

    virtual int state_count() const {
      return _n_states;
    }

    virtual int state_n_targets(int state) const {
      return _logAkl->function(state)->n_targets();
    }

    virtual const int * state_targets(int state) const {
      return _logAkl->function(state)->targets();
    }
  
    virtual bool valid_transition_params(int state, Params const & params) const {
      return _logAkl->validParams(state, params);
    }
  
    virtual bool valid_emission_params(int state, int slot, Params const & params) const {
      return _logEkb->validSlotParams(state, slot, params);
    }

    virtual Params * get_transition_params(int state) const {
      return _logAkl->getParams(state);
    }

    virtual void set_transition_params(int state, Params const & params) {
      _logAkl->setParams(state, params);
    }

    virtual Params * get_emission_params(int state, int slot) const {
      return _logEkb->getSlotParams(state, slot);
    }
  
    virtual void set_emission_params(int state, int slot, Params const & params) {
      _logEkb->setSlotParams(state, slot, params);
    }

    virtual void set_initial_probs(double * probs) {
      for (int i = 0; i < _n_states; ++i)
        _init_log_probs[i] = log(probs[i]);
    }
  
    virtual bool set_transition_covars(int state, int * idxs, int length) {
      return _logAkl->setCovars(state, idxs, length);
    }
  
    virtual bool set_emission_covars(int state, int slot, int * idxs, int length) {
      return _logEkb->setSlotCovars(state, slot, idxs, length);
    }
  
    virtual bool get_emission_option(int state, int slot, const char * name, double * out_value) const {
      return _logEkb->getSlotOption(state, slot, name, out_value);
    }
  
    virtual bool set_emission_option(int state, int slot, const char * name, double value) {
      return _logEkb->setSlotOption(state, slot, name, value);
    }
  
    virtual bool get_transition_option(int state, const char * name, double * out_value) const {
      return _logAkl->getOption(state, name, out_value);
    }
  
    virtual bool set_transition_option(int state, const char * name, double value) {
      return _logAkl->setOption(state, name, value);
    }
  
    double forward(Iter & iter, double * matrix) const {
      double * m_col, * m_col_prev;
      LogSum * logsum = LogSum::create(_n_states);
      iter.resetFirst();
    
      try {
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
        
      } catch (QHMMException & e) {
        // clean up
        delete logsum;
        
        e.stack.push_back("forward");
        throw;
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

    double backward(Iter & iter, double * matrix) const {
      double * m_col, * m_col_next;
      LogSum * logsum = LogSum::create(_n_states);
      
      /* border conditions @ position = N - 1*/
      m_col = matrix + (iter.length() - 1)*_n_states;
      for (int k = 0; k < _n_states; ++k)
        m_col[k] = 0; /* log(1) */

      try {
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
      } catch (QHMMException & e) {
        // clean up
        delete logsum;
        
        e.stack.push_back("backward");
        throw;
      }
      double loglik = logsum->compute();

      // clean-up
      delete logsum;
      
      return loglik;
    }

#define AT(M, I, J) M[(I) + (J)*rows]

    void viterbi(Iter & iter, int * path) const {
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
      iter.resetFirst();
      for (int l = 0; l < _n_states; ++l) {
        AT(matrix, l, 0) = (*_logEkb)(iter, l) + _init_log_probs[l];
        AT(backptr, l, 0) = -1; /* stop */
      }

      try {
        /* inner columns */
        m_col_prev = matrix;
        m_col = matrix + rows;
        b_col = backptr + rows;
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
      } catch (QHMMException & e) {
        // clean up
        delete[] matrix;
        delete[] backptr;
        
        e.stack.push_back("viterbi");
        throw;
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
      for (int l = iter.length() - 1; l > 0 ; --pptr, --l) {
        z = AT(backptr, z, l);
        *pptr = z;
        /* assert(prev >= 0); */
      }
      
      // clean up
      delete[] matrix;
      delete[] backptr;
    }
  
    void state_posterior(Iter & iter, const double * const fw, const double * const bk, double * matrix) const {
      /* posterior matrix is filled, state by state */
      LogSum * logsum = LogSum::create(_n_states);
      
      /* posterior_i,k = exp(fw[i,k] + bk[i,k] - logPx_i)
       
         logPx_i = log sum_k exp(fw[i, k] + bk[i,k])
       */
      
      /* TODO: optimize this */
      for (int i = 0; i < iter.length(); ++i) {
        
        /* compute local log-lik */
        logsum->clear();
        for (int j = 0; j < _n_states; ++j)
          logsum->store(fw[i*_n_states + j] + bk[i*_n_states + j]);
        double logPx = logsum->compute();
        
        /* fill posterior */
        for (int j = 0; j < _n_states; ++j)
          matrix[j*iter.length() + i] = exp(fw[i*_n_states + j] + bk[i*_n_states + j] - logPx);
      }
    }

    void local_loglik(Iter & iter, const double * const fw, const double * const bk, double * result) const {
      LogSum * logsum = LogSum::create(_n_states);

      for (int i = 0; i < iter.length(); ++i) {
        logsum->clear();
        for (int j = 0; j < _n_states; ++j)
          logsum->store(fw[i*_n_states + j] + bk[i*_n_states + j]);
        result[i] = logsum->compute();
      }

      delete logsum;
    }

    void transition_posterior(Iter & iter_at_target, const double * const fw, const double * const bk, double loglik, int n_src, const int * const src, int n_tgt, double * result) const {

      int index_tgt = iter_at_target.index();
      const double * const fw_src = fw + _n_states * (index_tgt - 1);
      const double * const bk_tgt = bk + _n_states * index_tgt;
      double * rptr = result;

      for (int isrc = 0; isrc < n_src; ++isrc) {
        int k = src[isrc];
        const int * const tgt = _logAkl->function(k)->targets();
        
        for (int itgt = 0; itgt < n_tgt; ++itgt, ++rptr) {
          int l = tgt[itgt];
          double log_emission = (*_logEkb)(iter_at_target, l);
          double log_trans = (*_logAkl)(iter_at_target, k, l);
          
          *rptr = exp(fw_src[k] + log_trans + log_emission + bk_tgt[l] - loglik);
        }
      }
    }
};

// auxiliary function to enable type inference
template <typename InnerFwd, typename InnerBck, typename FuncAkl, typename FuncEkb>
HMM * new_hmm_instance(InnerFwd innerFwd, InnerBck innerBck, FuncAkl logAkl, FuncEkb logEkb, double * init_log_probs) {
  return new HMMImpl<InnerFwd, InnerBck, FuncAkl, FuncEkb>(innerFwd, innerBck, logAkl, logEkb, init_log_probs);
}

#endif
