#include "em_base.hpp"
#include "em_seq.hpp"

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

PosteriorIterator * EMSequences::iterator(int state, int slot) {
  return new PosteriorIterator(state, slot, &_em_seqs);
}

TransitionPosteriorIterator * EMSequences::transition_iterator(std::vector<TransitionFunction*> & group) {
  return new TransitionPosteriorIterator(group, &_em_seqs);
}

double EMSequences::updateFwBk() {
  std::vector<EMSequence*>::iterator it;
  double loglik = 0;
  
  for (it = _em_seqs.begin(); it != _em_seqs.end(); ++it)
    loglik += (*it)->updateFwBk();
  
  return loglik;
}
