#ifndef EM_BASE_HPP
#define EM_BASE_HPP

#include <vector>

class HMM;
class EMSequence;

#include "post_iter.hpp"
#include "trans_post_iter.hpp"

class EMSequences {
public:
  EMSequences(HMM * hmm, std::vector<Iter*> & iters);
  ~EMSequences();

  PosteriorIterator * iterator(int state, int slot);
  
  TransitionPosteriorIterator * transition_iterator(std::vector<TransitionFunction*> & group);
  
  // returns sequence set log-likelihood
  double updateFwBk();
  
private:
  std::vector<EMSequence*> _em_seqs;
};

#endif
