#include "trans_post_iter.hpp"
#include "em_seq.hpp"

TransitionPosteriorIterator::TransitionPosteriorIterator(std::vector<TransitionFunction*> & group, const std::vector<EMSequence*> * seqs) : _iter((*seqs)[0]->iter()) {
  
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

TransitionPosteriorIterator::~TransitionPosteriorIterator() {
  delete[] _trans_post;
  delete[] _group_ids;
}

void TransitionPosteriorIterator::reset() {
  _seq_iter = _seqs->begin();
  changed_sequence();
  next(); /* will update values for first transition
           and cause iterator to be over second position
           as transitions take values at target
           */
}

bool TransitionPosteriorIterator::next() {
  bool res = _iter.next();
  
  // try to move to next sequence(s)
  while (!res && (++_seq_iter) != _seqs->end()) {
    changed_sequence();
    res = _iter.next();
  }
  
  if (res) {
    double logPxi = _local_logPx[_iter.index()]; // NOTE: RHMM used local Px at src not target ...
    (*_seq_iter)->hmm()->transition_posterior(_iter, _fw, _bk, logPxi, _group_size, _group_ids, _n_targets, _trans_post);
    
  }
  return res;
}

void TransitionPosteriorIterator::changed_sequence() {
  _iter = (*_seq_iter)->iter();
  _iter.resetFirst();
  _fw = (*_seq_iter)->forward();
  _bk = (*_seq_iter)->backward();
  _local_logPx = (*_seq_iter)->local_loglik();
}
