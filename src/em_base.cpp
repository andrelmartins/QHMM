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
  int n_threads = 1;

  #ifdef _OPENMP  
  omp_set_nested(1); // Allow nested OMP directives.
  n_threads = omp_get_num_threads(); // Get number of threads (set elsewhere).
  n_threads = (n_threads==1)?1:(n_threads/2); // Threads split between forward and backward during updates.
  #endif

  #pragma omp parallel for num_threads(n_threads)
  for (int i=0;i<_em_seqs.size();i++) {
//    try {
      loglik += (_em_seqs[i])->updateFwBk();
//    } catch (QHMMException & e) {
//      int seq_id = it - _em_seqs.begin();
//      e.sequence_id = seq_id;
//      throw;
//    }
  }
  
  return loglik;
}
