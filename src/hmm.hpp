#ifndef HMM_HPP
#define HMM_HPP

#include <vector>
#include "iter.hpp"
#include "func_table.hpp"

class HMM {
  public:
    virtual ~HMM() {}

    virtual bool valid_transition_params(int state, Params const & params) const = 0;
    virtual bool valid_emission_params(int state, int slot, Params const & params) const = 0;
  
    virtual Params * get_transition_params(int state) const = 0;
    virtual void set_transition_params(int state, Params const & params) = 0;
    virtual Params * get_emission_params(int state, int slot) const = 0;
    virtual void set_emission_params(int state, int slot, Params const & params) = 0;
    virtual void set_initial_probs(double * probs) = 0;
    virtual bool set_transition_covars(int state, int * idxs, int length) = 0;
    virtual bool set_emission_covars(int state, int slot, int * idxs, int length) = 0;
    virtual bool get_emission_option(int state, int slot, const char * name, double * out_value) = 0;
    virtual bool set_emission_option(int state, int slot, const char * name, double value) = 0;
    virtual bool get_transition_option(int state, const char * name, double * out_value) = 0;
    virtual bool set_transition_option(int state, const char * name, double value) = 0;

    virtual double forward(Iter & iter, double * matrix) = 0;
    virtual double backward(Iter & iter, double * matrix) = 0;
    virtual void viterbi(Iter & iter, int * path) = 0;
    virtual void state_posterior(Iter & iter, const double * const fw, const double * const bk, double * matrix) = 0;
    virtual void local_loglik(Iter & iter, const double * const fw, const double * const bk, double * result) const = 0;
    virtual void transition_posterior(Iter & iter_at_target, const double * const fw, const double * const bk, double loglik, int n_src, const int * const src, int n_tgt, double * result) const = 0;

    virtual double em(std::vector<Iter*> & iters, double tolerance);

    template<typename TransTableT, typename EmissionTableT>
    static HMM * create(TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs);

    // properties
    virtual int state_count() const = 0;
  
protected:
    virtual const std::vector<std::vector<EmissionFunction*> > & emission_groups() const = 0;
    virtual const std::vector<std::vector<TransitionFunction*> > & transition_groups() const = 0;
};

#include "hmm.cpp"

#endif
