#include "hmm.hpp"
#include <limits>

double HMM::em(std::vector<Iter*> iters, double tolerance) {
  int iter_count = 0;
  double cur_loglik, prev_loglik;

  /* if there is any missing data, split sequences for posterior/parameter updates ??? */

  /* allocate memory for forward/backward & posterior */

  /* initialize sufficient statistic instances */


  /* main EM loop */
  prev_loglik = -std::numeric_limits<double>::infinity();
  while (1) {
    ++iter_count;

    /* reset sufficient statistics */

    /* compute forward/backward & posterior per sequence => get log-lik */

    /* output cur_loglik & current parameters */

    /* check log-lik */
    if (cur_loglik < prev_loglik ||
        cur_loglik - prev_loglik < tolerance)
      break;

    /* update parameters */

    /* */
    prev_loglik = cur_loglik;
  }

  return cur_loglik;
}
