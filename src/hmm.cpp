#include "hmm_tmpl.hpp"

template<typename TransTableT, typename EmissionTableT>
HMM * HMM::Create(TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs) {

  // determine appropriate inner loop type (Sparse vs Dense)
  if (transitions->isSparse()) {
    InnerFwdSparse<TransTableT *> * innerFwd = new InnerFwdSparse<TransTableT *>(transitions);
    InnerBckSparse<TransTableT *, EmissionTableT *> * innerBck = new InnerBckSparse<TransTableT *, EmissionTableT *>(transitions);
    
    return new_hmm_instance(innerFwd, innerBck, transitions, emissions, init_log_probs);
  }
  
  // dense
  InnerFwdDense<TransTableT *> * innerFwd = new InnerFwdDense<TransTableT *>();
  InnerBckDense<TransTableT *, EmissionTableT *> * innerBck = new InnerBckDense<TransTableT *, EmissionTableT *>();
  
  return new_hmm_instance(innerFwd, innerBck, transitions, emissions, init_log_probs);
}
