#include "em_base.hpp"
#include "em_seq.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

EMSequences::EMSequences(HMM * hmm, std::vector<Iter*> & iters) {
  std::vector<Iter*>::iterator it;
  _unitarySequences = true;

  for (it = iters.begin(); it != iters.end(); ++it) {
    EMSequence * seq = new EMSequence(hmm, (*it));
    _em_seqs.push_back(seq);

    if ((*it)->length() > 1)
      _unitarySequences = false;
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
  QHMMThreadHelper helper;

  #pragma omp parallel shared(helper, loglik)
  {
    #pragma omp for
    for (unsigned int i = 0; i < _em_seqs.size(); ++i)
      (_em_seqs[i])->updateFwBk(i, helper, loglik);
  }

  // helper will rethrow exceptions on exit
  return loglik;
}
