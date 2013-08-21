#include "hmm.hpp"
#include "em_base.hpp"
#include <limits>

double HMM::em(std::vector<Iter*> & iters, double tolerance) {
  int iter_count = 0;
  double cur_loglik, prev_loglik;

  /* initialize sequences & fw/bk memory
     (handles spliting by missing data)
   */
  EMSequences * sequences = new EMSequences(this, iters);

  /* initialize sufficient statistic instances */


  /* main EM loop */
  prev_loglik = -std::numeric_limits<double>::infinity();
  while (1) {
    ++iter_count;

    /* reset sufficient statistics */

    /* compute forward/backward per sequence => get log-lik */
    cur_loglik = sequences->updateFwBk();

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
  delete sequences;

  return cur_loglik;
}
