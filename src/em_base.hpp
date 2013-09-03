#ifndef EM_BASE_HPP
#define EM_BASE_HPP

#include "hmm.hpp"
#include "iter.hpp"
#include <vector>

class EMSequences;

class EMSequence {
public:
  EMSequence(HMM * hmm, Iter * iter);
  ~EMSequence();

  // returns sequence log-likelihood
  double updateFwBk();
  
  // accessors
  const double * forward() { return _forward; }
  const double * backward() { return _backward; }
  Iter & iter() { return *_iter; }
  const HMM * hmm() { return _hmm; }
  const double * local_loglik();
  
  friend class EMSequences;

private:
  Iter * _iter;
  HMM * _hmm;
  
  bool _posterior_dirty;
  double * _forward;
  double * _backward;
  double * _posterior;
  std::vector<std::vector<Iter>* > * _slot_subiters;

  void update_posterior();

  bool _local_loglik_dirty;
  double * _local_loglik;
};

#include "trans_post_iter.hpp"

class EMSequences {
public:
  EMSequences(HMM * hmm, std::vector<Iter*> & iters);
  ~EMSequences();

  class PosteriorIterator {
  public:
    PosteriorIterator(int state, int slot, const std::vector<EMSequence*> * seqs) : _state(state), _slot(slot), _seqs(seqs) {
      reset();
    }

    void reset() {
      _outer_iter = _seqs->begin();
      _inner_vec = (*((*_outer_iter)->_slot_subiters))[_slot];
      _inner_iter = _inner_vec->begin();
      (*_outer_iter)->update_posterior();
    }

    bool next() {
      /* try to move inner iterator */
      ++_inner_iter;
      if (_inner_iter == _inner_vec->end()) {
        /* try to move outer iterator */
        ++_outer_iter;
        if (_outer_iter == _seqs->end())
          return false;

        /* update posterior */
        (*_outer_iter)->update_posterior();

        /* update inner iterator */
        _inner_vec = (*((*_outer_iter)->_slot_subiters))[_slot];
        _inner_iter = _inner_vec->begin();
      }

      return true;
    }

    Iter & iter() {
      return *_inner_iter;
    }

    const double * posterior() {
      EMSequence * em_seq = (*_outer_iter);
      double * post_state = em_seq->_posterior + (_state * em_seq->_iter->length());
      return post_state + (*_inner_iter).iter_offset();
    }

  private:
    const int _state;
    const int _slot;

    std::vector<EMSequence*>::const_iterator _outer_iter;
    std::vector<Iter> * _inner_vec;
    std::vector<Iter>::iterator _inner_iter;

    const std::vector<EMSequence*> * _seqs;
  };

  PosteriorIterator * iterator(int state, int slot) {
    return new PosteriorIterator(state, slot, &_em_seqs);
  }

  TransitionPosteriorIterator * transition_iterator(std::vector<TransitionFunction*> & group) {
    return new TransitionPosteriorIterator(group, &_em_seqs);
  }

  // returns sequence log-likelihood
  double updateFwBk() {
    std::vector<EMSequence*>::iterator it;
    double loglik = 0;
    
    for (it = _em_seqs.begin(); it != _em_seqs.end(); ++it)
      loglik += (*it)->updateFwBk();
    
    return loglik;
  }

private:
  std::vector<EMSequence*> _em_seqs;
};

#endif
