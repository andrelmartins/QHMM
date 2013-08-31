#ifndef POST_ITER_HPP
#define POST_ITER_HPP

#include "hmm.hpp"

// Iterator for posterior transitions
class PostIter {
public:
  PostIter(const double * forward, const double * backward, 
	   const double * local_loglik,
	   std::vector<TransitionFunction*> & group, HMM * hmm,
	   Iter & iter) : _fw(forward), _bk(backward), 
			  _local_logPx(local_loglik), _hmm(hmm),
        _iter(iter) {
    _group_size = group.size();
    _group_ids = new int[_group_size];
    for (unsigned int i = 0; i < _group_size; ++i)
      _group_ids[i] = group[i]->stateID();

    TransitionFunction * head = group[0];

    _n_targets = head->n_targets(); /* all functions in a group must share
				       the same number of targets */

    _trans_post = new double[_n_targets * _group_size];
    reset();
  }

  ~PostIter() {
    delete[] _trans_post;
    delete[] _group_ids;
  }

  void reset() {
    _iter.resetFirst();
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
  unsigned int _group_size;
  int * _group_ids;
  int _n_targets;
  Iter & _iter;
  double * _trans_post;
};

#endif
