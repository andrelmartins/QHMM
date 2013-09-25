#include "hmm.hpp"
#include "em_base.hpp"
#include <limits>
#include <cstdio>
#include <cmath>

std::vector<ParamRecord*> * HMM::init_records() const {
  std::vector<ParamRecord*> * result = new std::vector<ParamRecord*>();

  /* add records for transitions
     (only need to look at head of each group)
  */

  std::vector<std::vector<TransitionFunction*> > tgroups = transition_groups();
  std::vector<std::vector<TransitionFunction*> >::iterator tit;

  for (tit = tgroups.begin(); tit != tgroups.end(); ++tit) {
    std::vector<TransitionFunction*> group_i = *tit;
    TransitionFunction * head = group_i[0];

    Params * par = head->getParams();
    if (par != NULL && !par->isAllFixed())
      result->push_back(new ParamRecord(head));

    if (par != NULL)
      delete par;
  }

  /* add records for emissions */
  std::vector<std::vector<EmissionFunction*> > groups = emission_groups();
  std::vector<std::vector<EmissionFunction*> >::iterator it;

  for (it = groups.begin(); it != groups.end(); ++it) {
    std::vector<EmissionFunction*> group_i = *it;
    EmissionFunction * head = group_i[0];
    
    Params * par = head->getParams();
    if (par != NULL && !par->isAllFixed())
      result->push_back(new ParamRecord(head));

    if (par != NULL)
      delete par;
  }    

  return result;
}

void HMM::update_records(std::vector<ParamRecord*> * ptr) const {
  for (unsigned int i = 0; i < ptr->size(); ++i)
    (*ptr)[i]->collect();
}

void HMM::delete_records(std::vector<ParamRecord*> * ptr) {
  for (unsigned int i = 0; i < ptr->size(); ++i)
    delete (*ptr)[i];

  delete ptr;
}

EMResult HMM::em(std::vector<Iter*> & iters, double tolerance) {
  int iter_count = 0;
  double cur_loglik, prev_loglik;
  EMResult result;
  bool skip_transitions;

  /* initialize result trace */
  result.param_trace = init_records();

  /* initialize sequences & fw/bk memory
     (handles spliting by missing data)
   */
  EMSequences * sequences = new EMSequences(this, iters);
  skip_transitions = sequences->unitarySequences();

  /* main EM loop */
  prev_loglik = -std::numeric_limits<double>::infinity();
  while (1) {
    ++iter_count;

    /* compute forward/backward per sequence => get log-lik */
    cur_loglik = sequences->updateFwBk();

    /* output cur_loglik & store current parameters */
    update_records(result.param_trace);
    printf("[%d] loglik: %g\n", iter_count, cur_loglik);

    /* check log-lik */
    if (cur_loglik < prev_loglik ||
        cur_loglik - prev_loglik < tolerance)
      break;

    /* update parameters */

    /* - transition functions */
    if (!skip_transitions) {
      std::vector<std::vector<TransitionFunction*> > tgroups = transition_groups();
      std::vector<std::vector<TransitionFunction*> >::iterator tit;
      for (tit = tgroups.begin(); tit != tgroups.end(); ++tit) {
	std::vector<TransitionFunction*> group_i = *tit;
	TransitionFunction * head = group_i[0];
	
	head->updateParams(sequences, &group_i);
      }
      refresh_transition_table(); // refresh internal caches
    }

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

  result.log_likelihood = cur_loglik;

  return result;
}
