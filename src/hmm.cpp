#include "hmm_tmpl.hpp"

template<typename TransTableT, typename EmissionTableT>
HMM * HMM::Create(TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs) {

  // determine appropriate inner loop type (Sparse vs Dense)
  if (transitions->isSparse()) {
    int ** previous = transitions->previousStates();
    InnerFwdSparse<TransTableT *> * innerFwd = new InnerFwdSparse<TransTableT *>(previous, transitions->n_states());
    
    return new_hmm_instance(innerFwd, transitions, emissions, init_log_probs);
  }
  
  // dense
  return new_hmm_instance(new InnerFwdDense<TransTableT *>(), transitions, emissions, init_log_probs);
}
