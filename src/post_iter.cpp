#include "post_iter.hpp"
#include "em_seq.hpp"

PosteriorIterator::PosteriorIterator(int state, int slot, const std::vector<EMSequence*> * seqs) : _state(state), _slot(slot), _seqs(seqs) {
  reset();
}

void PosteriorIterator::reset() {
  _outer_iter = _seqs->begin();
  _inner_vec = (*((*_outer_iter)->_slot_subiters))[_slot];
  _inner_iter = _inner_vec->begin();
  (*_outer_iter)->update_posterior();
}

bool PosteriorIterator::next() {
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

Iter & PosteriorIterator::iter() {
  return *_inner_iter;
}

const double * PosteriorIterator::posterior() {
  EMSequence * em_seq = (*_outer_iter);
  double * post_state = em_seq->_posterior + (_state * em_seq->_iter->length());
  return post_state + (*_inner_iter).iter_offset();
}
