#ifndef HMM_HPP
#define HMM_HPP

#include "iter.hpp"
#include "func_table.hpp"

class HMM {
  public:
    virtual ~HMM() {}

    virtual double forward(Iter & iter, double * matrix) = 0;
    virtual double backward(Iter & iter, double * matrix) = 0;
    virtual void viterbi(Iter & iter, int * path) = 0;

    template<typename TransTableT, typename EmissionTableT>
    static HMM * create(TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs);

};

#include "hmm.cpp"

#endif
