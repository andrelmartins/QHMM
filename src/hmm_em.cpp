#include "hmm.hpp"
#include "em_base.hpp"
#include <limits>

double HMM::em(std::vector<Iter*> & iters, double tolerance) {
  int iter_count = 0;
  double cur_loglik, prev_loglik;

  /* initialize sequences & fw/bk memory
     (handles spliting by missing data)
   */
  std::vector<EMSequence*> sequences;
  for (std::vector<Iter*>::iterator it = iters.begin();
       it != iters.end();
       ++it)
    sequences.push_back(new EMSequence(this, *it));

  /* initialize sufficient statistic instances */

  /* . emissions */
  std::vector<std::vector<EmissionSuffStat*> > single_pass_emissions;

  for (int i = 0; i < this->n_slots(); ++i) {
    std::vector<EmissionSuffStat*> slot_i_single;

    for (int k = 0; k < this->state_count(); ++k) {
      EmissionSuffStat * ss_i = this->emission_suff_stats_instance(k, i);

      if (ss_i != NULL) {
	if (ss_i->suff_stats_type() == SinglePass)
	  slot_i_single.push_back(ss_i);
	else if (ss_i->suff_stats_type() == MultiPass) {
	  /* TODO: implement this */
	}
      }
    }

    single_pass_emissions.push_back(slot_i_single);
  }

  /* main EM loop */
  prev_loglik = -std::numeric_limits<double>::infinity();
  while (1) {
    ++iter_count;

    /* reset sufficient statistics */

    /* compute forward/backward per sequence => get log-lik */
    cur_loglik = 0;
    for (std::vector<EMSequence*>::iterator it = sequences.begin();
       it != sequences.end();
       ++it)
      cur_loglik += (*it)->updateFwBk();

    /* output cur_loglik & current parameters */

    /* check log-lik */
    if (cur_loglik < prev_loglik ||
        cur_loglik - prev_loglik < tolerance)
      break;

    /* update parameters */

    /* */
    prev_loglik = cur_loglik;
  }

  /* clean up */
  for (std::vector<EMSequence*>::iterator it = sequences.begin();
       it != sequences.end();
       ++it)
     delete *it;

  return cur_loglik;
}
