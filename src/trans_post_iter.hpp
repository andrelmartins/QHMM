#ifndef TRANS_POST_ITER_HPP
#define TRANS_POST_ITER_HPP

#include "iter.hpp"
#include "base_classes.hpp"

class EMSequence;

// Iterator for posterior transitions
class TransitionPosteriorIterator {
public:
  TransitionPosteriorIterator(std::vector<TransitionFunction*> & group, const std::vector<EMSequence*> * seqs);
  ~TransitionPosteriorIterator();

  bool next();
  
  void reset();
  
  double posterior(int grp_idx, int tgt_idx) {
    return _trans_post[grp_idx * _n_targets + tgt_idx];
  }

  /* implement covariate interface if needed */
  double covar(const int slot) const {
    return _iter.covar(slot);
  }
    
  double covar_i(const int slot, const int i) const {
    return _iter.covar_i(slot, i);
  }
  
  // offset with respect to current iterator sequence position
  double covar_ext(const int slot, const int i, const int offset) const {
    return _iter.covar_ext(slot, i, offset);
  }

  int length() const { return _iter.length(); }

private:
  const double * _fw;
  const double * _bk;
  const double * _local_logPx;
  Iter & _iter;
  unsigned int _group_size;
  int * _group_ids;
  int _n_targets;
  double * _trans_post;
  
  const std::vector<EMSequence*> * _seqs;
  std::vector<EMSequence*>::const_iterator _seq_iter;
  
  void changed_sequence();
};

#endif
