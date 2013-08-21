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
}

EMSequence::~EMSequence() {
  for (int i = 0; i < _slot_subiters->size(); ++i)
    delete (*_slot_subiters)[i];
  delete _slot_subiters;

  delete[] _forward;
  delete[] _backward;
  if (_posterior != NULL)
    delete[] _posterior;
}

double EMSequence::updateFwBk() {
  double loglik = 0;
  
  loglik = _hmm->forward(*_iter, _forward);
  loglik = _hmm->backward(*_iter, _backward); // ideally should check both are eq

  /* we don't update the posterior here just in case all emissions
     are fixed and we don't actually need it
  */
  _posterior_dirty = true; /* needs update */
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
