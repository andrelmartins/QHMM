#ifndef EM_BASE_HPP
#define EM_BASE_HPP

#include "hmm.hpp"
#include "iter.hpp"
#include <vector>

class EMSequence {
public:
  EMSequence(HMM * hmm, Iter * iter);
  ~EMSequence();

  // returns sequence log-likelihood
  double updateFwBk();

  // call suff stat objects (for given slot) with posteriors for
  // non-missing emission data segments of this sequence
  void emission_suff_stats(int slot, std::vector<EmissionSuffStat*> * ssobjs);

private:
  Iter * _iter;
  HMM * _hmm;
  
  bool _posterior_dirty;
  double * _forward;
  double * _backward;
  double * _posterior;
  std::vector<std::vector<Iter>* > * _slot_subiters;

};

#endif
