#include "em_base.hpp"

EMSequence::EMSequence(HMM * hmm, Iter * iter) {
  /* keep pointer to main iterator and HMM */
  _iter = iter;
  _hmm = hmm;

  /* initialize sub-iterators */
  int n_slots = iter->emission_slot_count();

  _slot_subiters = new std::vector<std::vector<Iter>* >();
  for (int i = 0; i < n_slots; ++i)
    _slot_subiters->push_back(iter->sub_iterators(i));

  /* allocate space for forward, backward */
  int n_states = hmm->state_count();
  _forward = new double[n_states * iter->length()];
  _backward = new double[n_states * iter->length()];

  _posterior = NULL; /* only allocate posterior on first use */
  _posterior_dirty = true; /* needs update */  
  _local_loglik = NULL; /* only allocate on first use */
  _local_loglik_dirty = true; /* needs update */
}

EMSequence::~EMSequence() {
  for (unsigned int i = 0; i < _slot_subiters->size(); ++i)
    delete (*_slot_subiters)[i];
  delete _slot_subiters;

  delete[] _forward;
  delete[] _backward;
  if (_posterior != NULL)
    delete[] _posterior;
  if (_local_loglik != NULL)
    delete[] _local_loglik;
}

double EMSequence::updateFwBk() {
  double loglik = 0;
  
  loglik = _hmm->forward(*_iter, _forward);
  loglik = _hmm->backward(*_iter, _backward); // ideally should check both are eq

  /* we don't update the posterior here just in case all emissions
     are fixed and we don't actually need it
  */
  _posterior_dirty = true; /* needs update */
  _local_loglik_dirty = true; /* needs update */
  return loglik;
}

void EMSequence::update_posterior() {
  /* update posterior if needed */
  if (_posterior_dirty) {
    /* allocate posterior if needed */
    if (_posterior == NULL) {
      int n_states = _hmm->state_count();
      _posterior = new double[n_states * _iter->length()];
    }

    _hmm->state_posterior(*_iter, _forward, _backward, _posterior);
    _posterior_dirty = false;
  }
}

PostIter * EMSequence::transition_iterator(std::vector<TransitionFunction*> & group) {
  return new PostIter(_forward, _backward, local_loglik(), group, _hmm, *_iter);
}

double * EMSequence::local_loglik() {
  if (_local_loglik_dirty) {
    if (_local_loglik == NULL)
      _local_loglik = new double[_iter->length()];

    _hmm->local_loglik(*_iter, _forward, _backward, _local_loglik);
    _local_loglik_dirty = false;
  }

  return _local_loglik;
}

EMSequences::EMSequences(HMM * hmm, std::vector<Iter*> & iters) {
  std::vector<Iter*>::iterator it;

  for (it = iters.begin(); it != iters.end(); ++it) {
    EMSequence * seq = new EMSequence(hmm, (*it));
    _em_seqs.push_back(seq);
  }
}

EMSequences::~EMSequences() {
  std::vector<EMSequence*>::iterator it;
  
  for (it = _em_seqs.begin(); it != _em_seqs.end(); ++it)
    delete (*it);
}
