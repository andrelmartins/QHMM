#ifndef POST_ITER_HPP
#define POST_ITER_HPP

#include "iter.hpp"

class EMSequence;

class PosteriorIterator {
public:
  PosteriorIterator(int state, int slot, const std::vector<EMSequence*> * seqs);
  void reset();
  
  bool next();
  
  Iter & iter();
  
  const double * posterior();
  
private:
  const int _state;
  const int _slot;
  
  std::vector<EMSequence*>::const_iterator _outer_iter;
  std::vector<Iter> * _inner_vec;
  std::vector<Iter>::iterator _inner_iter;
  
  const std::vector<EMSequence*> * _seqs;
};

#endif
