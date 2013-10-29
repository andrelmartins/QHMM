#ifndef HMM_HPP
#define HMM_HPP

#include <vector>
#include "iter.hpp"
#include "base_func_table.hpp"
#include "param_record.hpp"

typedef struct EMResult {
  std::vector<double> * log_likelihood;
  std::vector<ParamRecord*> * param_trace;
} EMResult;

class HMM {
  public:
    virtual ~HMM() {}

    virtual TransitionTable * transitions() const = 0;
    virtual EmissionTable * emissions() const = 0;
  
    virtual void set_initial_probs(double * probs) = 0;

    virtual double forward(Iter & iter, double * matrix) const = 0;
    virtual double backward(Iter & iter, double * matrix) const = 0;
    virtual void viterbi(Iter & iter, int * path) const = 0;
    virtual void state_posterior(Iter & iter, const double * const fw, const double * const bk, double * matrix) const = 0;
    virtual void local_loglik(Iter & iter, const double * const fw, const double * const bk, double * result) const = 0;
    virtual void transition_posterior(Iter & iter_at_target, const double * const fw, const double * const bk, double loglik, int n_src, const int * const src, int n_tgt, double * result) const = 0;

    virtual struct EMResult em(std::vector<Iter*> & iters, double tolerance);

    template<typename TransTableT, typename EmissionTableT>
    static HMM * create(TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs);

    // properties
    virtual int state_count() const = 0;
  
    static void delete_records(std::vector<ParamRecord*> * ptr);
protected:
    virtual const std::vector<std::vector<EmissionFunction*> > & emission_groups() const = 0;
    virtual const std::vector<std::vector<TransitionFunction*> > & transition_groups() const = 0;
    virtual void refresh_transition_table() = 0;

    std::vector<ParamRecord*> * init_records() const;
    void update_records(std::vector<ParamRecord*> * ptr) const;
};

#include "hmm.cpp"

#endif
