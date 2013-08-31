#ifndef POST_ITER_HPP
#define POST_ITER_HPP

#include "hmm.hpp"

// Iterator for posterior transitions
class TransitionPosteriorIterator {
public:
  TransitionPosteriorIterator(std::vector<TransitionFunction*> & group, const std::vector<EMSequence*> * seqs) : _hmm((*seqs)[0]->hmm()), _iter((*seqs)[0]->iter()) {
    
    // initialize sequence iterator
    _seqs = seqs;
    _seq_iter = seqs->begin();
    
    // initialize group information
    _group_size = group.size();
    _group_ids = new int[_group_size];
    for (unsigned int i = 0; i < _group_size; ++i)
      _group_ids[i] = group[i]->stateID();

    // initialize target information
    TransitionFunction * head = group[0];
    _n_targets = head->n_targets(); /* all functions in a group must share
                                     the same number of targets */
    
    // internal memory
    _trans_post = new double[_n_targets * _group_size];
    
    // initialize to first position
    reset();
  }

  ~TransitionPosteriorIterator() {
    delete[] _trans_post;
    delete[] _group_ids;
  }

  void reset() {
    _seq_iter = _seqs->begin();
    changed_sequence();
    next(); /* will update values for first transition
	       and cause iterator to be over second position
	       as transitions take values at target
	     */
  }

  double posterior(int grp_idx, int tgt_idx) {
    return _trans_post[grp_idx * _n_targets + tgt_idx];
  }

  bool next() {
    bool res = _iter.next();
    
    // try to move to next sequence(s)
    while (!res && (++_seq_iter) != _seqs->end()) {
      changed_sequence();
      res = _iter.next();
    }

    if (res) {
      double logPxi = _local_logPx[_iter.index()]; // NOTE: RHMM used local Px at src not target ...
      _hmm->transition_posterior(_iter, _fw, _bk, logPxi, _group_size, _group_ids, _n_targets, _trans_post);

    }
    return res;
  }

  /* implement covariate interface if needed */
  double covar(const int slot) const {
    return _iter.covar(slot);
  }
    
  double covar_i(const int slot, const int i) const {
    return _iter.covar_i(slot, i);
  }
  
  // offset with respect to current iterator sequence position
  double covar_ext(const int slot, const int i, const int offset) const {
    return _iter.covar_ext(slot, i, offset);
  }

  int length() const { return _iter.length(); }

private:
  const double * _fw;
  const double * _bk;
  const double * _local_logPx;
  const HMM * _hmm;
  Iter & _iter;
  unsigned int _group_size;
  int * _group_ids;
  int _n_targets;
  double * _trans_post;
  
  const std::vector<EMSequence*> * _seqs;
  std::vector<EMSequence*>::const_iterator _seq_iter;
  
  void changed_sequence() {
    _iter = (*_seq_iter)->iter();
    _iter.resetFirst();
    _fw = (*_seq_iter)->forward();
    _bk = (*_seq_iter)->backward();
    _local_logPx = (*_seq_iter)->local_loglik();
    _hmm = (*_seq_iter)->hmm();
  }
};

#endif
