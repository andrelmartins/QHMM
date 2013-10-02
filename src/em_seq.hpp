#ifndef EM_SEQ_HPP
#define EM_SEQ_HPP

#include "QHMMThreadHelper.hpp"
#include "hmm.hpp"
#include "iter.hpp"
#include <vector>

class PosteriorIterator;

class EMSequence {
public:
  EMSequence(HMM * hmm, Iter * iter);
  ~EMSequence();
  
  // returns sequence log-likelihood
  void updateFwBk(int seq_id, QHMMThreadHelper & helper, double & loglik);
  
  // accessors
  const double * forward() { return _forward; }
  const double * backward() { return _backward; }
  Iter & iter() { return *_iter; }
  const HMM * hmm() { return _hmm; }
  const double * local_loglik();
  
  friend class PosteriorIterator;
  
private:
  Iter * _iter;
  Iter * _iterCopy;
  HMM * _hmm;
  
  bool _posterior_dirty;
  double * _forward;
  double * _backward;
  double * _posterior;
  std::vector<std::vector<Iter>* > * _slot_subiters;
  
  void update_posterior();
  
  bool _local_loglik_dirty;
  double * _local_loglik;
};

#endif
