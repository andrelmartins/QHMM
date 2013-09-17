#include "hmm.hpp"
#include "em_base.hpp"
#include <limits>
#include <cstdio>
#include <cmath>

double HMM::em(std::vector<Iter*> & iters, double tolerance) {
  int iter_count = 0;
  double cur_loglik, prev_loglik;

  /* initialize sequences & fw/bk memory
     (handles spliting by missing data)
   */
  EMSequences * sequences = new EMSequences(this, iters);

  /* main EM loop */
  prev_loglik = -std::numeric_limits<double>::infinity();
  while (1) {
    ++iter_count;

    /* compute forward/backward per sequence => get log-lik */
    cur_loglik = sequences->updateFwBk();

    /* output cur_loglik & current parameters */
    printf("[%d] loglik: %g\n", iter_count, cur_loglik);

    /* check log-lik */
    if (cur_loglik < prev_loglik ||
        cur_loglik - prev_loglik < tolerance)
      break;

    /* update parameters */

    /* - transition functions */
    std::vector<std::vector<TransitionFunction*> > tgroups = transition_groups();
    std::vector<std::vector<TransitionFunction*> >::iterator tit;
    for (tit = tgroups.begin(); tit != tgroups.end(); ++tit) {
      std::vector<TransitionFunction*> group_i = *tit;
      TransitionFunction * head = group_i[0];
      
      head->updateParams(sequences, &group_i);
    }
    refresh_transition_table(); // refresh internal caches

    /* - emission functions */
    std::vector<std::vector<EmissionFunction*> > groups = emission_groups();
    std::vector<std::vector<EmissionFunction*> >::iterator it;
    for (it = groups.begin(); it != groups.end(); ++it) {
      std::vector<EmissionFunction*> group_i = *it;
      EmissionFunction * head = group_i[0];
      
      head->updateParams(sequences, &group_i);
    }
   
    /* */
    prev_loglik = cur_loglik;
    
    /* check things went ok */
    if (prev_loglik == -std::numeric_limits<double>::infinity()) {
      printf("-Inf log likelihood! Aborted!\n");
      break;
    }
    if (std::isnan(prev_loglik)) {
      printf("NaN log likelihood! Aborted!\n");
      break;
    }
  }

  /* clean up */
  delete sequences;

  return cur_loglik;
}
