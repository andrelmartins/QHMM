/*
  Benchmark for forward algorithm and different choices of HMM
  implementation details.
*/

#include <iostream>
#include <cstdlib>
#include "../src/hmm.hpp"
#include "../src/transitions/autocorr.hpp"
#include "../src/emissions/poisson.hpp"

using namespace std;

void fill_sequence(double * array, int length, double low, double high, double step) {
  double value = low;

  for (int i = 0; i < length; ++i, ++array) {
    *array = value;

    value += step;
    if (value > high)
      value = low;
  }
}

double * seq_merge(double * array1, double * array2, int length) {
  int res_len = 2 * length;
  double * result = new double[2*length];

  for (int i = 0, j = 0; i < res_len; i+=2, ++j) {
    result[i] = array1[j];
    result[i + 1] = array2[j];
  }
   
  return result;
}

/* homogeneous transition table */
HomogeneousTransitions * create_homogeneous_transitions(int n_states, double alpha, int k) {
  HomogeneousTransitions * result = new HomogeneousTransitions(n_states);
  
  int n_targets = 2 + k; /* k = 0 :: self + next */

  if (n_targets > n_states)
    n_targets = n_states;

  if (n_states == 1) {
    n_targets = 1;
    int target = 0;

    AutoCorr * state = new AutoCorr(alpha, n_states, n_targets, &target);
    result->insert(state);
    result->updateTransitions();
    return result;
  }

  int * targets = new int[n_targets];

  /* add states */
  for (int i = 0; i < n_states; ++i) {
    AutoCorr * state_i;

    /* fill target vector */
    for (int j = 0; j < n_targets; ++j)
      targets[j] = (i + j) % n_states;
    
    state_i = new AutoCorr(alpha, n_states, n_targets, targets);
    result->insert(state_i);
  }

  delete[] targets;

  /* fill in cached transition table */
  result->updateTransitions();

  return result;
}

/* non-homogeneous transition table */
NonHomogeneousTransitions * create_nonhomogeneous_transitions(int n_states, int covar_slot, int k) {
  NonHomogeneousTransitions * result = new NonHomogeneousTransitions(n_states);
  
  int n_targets = 2 + k; /* k = 0 :: self + next */
  
  if (n_targets > n_states)
    n_targets = n_states;
  
  if (n_states == 1) {
    n_targets = 1;
    int target = 0;
    
    AutoCorrCovar * state = new AutoCorrCovar(n_states, covar_slot, n_targets, &target);
    result->insert(state);
    return result;
  }
  
  int * targets = new int[n_targets];
  
  /* add states */
  for (int i = 0; i < n_states; ++i) {
    AutoCorrCovar * state_i;
    
    /* fill target vector */
    for (int j = 0; j < n_targets; ++j)
      targets[j] = (i + j) % n_states;
    
    state_i = new AutoCorrCovar(n_states, covar_slot, n_targets, targets);
    result->insert(state_i);
  }
  
  delete[] targets;
  return result;
}


Emissions * create_poisson_emissions(int n_states, double lambda) {
  Emissions * result = new Emissions(n_states);

  for (int i = 0; i < n_states; ++i)
    result->insert(new Poisson(lambda));

  return result;
}

Emissions * create_poisson_covar_emissions(int n_states, int covar_slot) {
  Emissions * result = new Emissions(n_states);
  
  for (int i = 0; i < n_states; ++i)
    result->insert(new PoissonCovar(covar_slot));
  
  return result;
}

MultiEmissions * create_poisson_both_emissions(int n_states, double lambda, int covar_slot) {
  MultiEmissions * result = new MultiEmissions(n_states, 2);
  
  for (int i = 0; i < n_states; ++i) {
    std::vector<EmissionFunction *> funcs;
    funcs.push_back(new Poisson(lambda));
    funcs.push_back(new PoissonCovar(covar_slot));
    
    result->insert(funcs);
  }
  
  return result;
}

void print_matrix(double * ptr, int n, int m) {
  for (int j = 0; j < m; ++j) {
    cout << *ptr;
    ++ptr;
    for (int i = 1; i < n; ++i, ++ptr) {
      cout << " " << *ptr;
    }
    cout << endl;
  }
}

template<typename TransTableT, typename EmissionTableT>
void run_test(Iter & iter, TransTableT * transitions, EmissionTableT * emissions, double * init_log_probs, double * fwd, int repeats, int test_number, bool try_sparse = true) {

  double loglik = 0;
  clock_t start, end;

  //
  // HMM with default settings (decided by HMM::Create)
  //
  HMM * hmm_auto = HMM::create(transitions, emissions, init_log_probs);
  
  
  start = clock();
  for (int r = 0; r < repeats; ++r)
    loglik += hmm_auto->forward(iter, fwd);
  
  end = clock();
  
  cout << "T" << test_number << ":A\t" << loglik/repeats << "\t" << (end - start)*1000.0/CLOCKS_PER_SEC << endl;
  
  delete hmm_auto;
  
  if (try_sparse) {
    //
    // HMM with explicit sparse mode
    //
    HMM * hmm_sparse = new_hmm_instance(new InnerFwdSparse<TransTableT *>(transitions),
                                        new InnerBckSparse<TransTableT *, EmissionTableT *>(transitions),
                                        transitions,
                                        emissions,
                                        init_log_probs);
    loglik = 0;
    start = clock();
    for (int r = 0; r < repeats; ++r)
      loglik += hmm_sparse->forward(iter, fwd);
    
    end = clock();
    cout << "T" << test_number << ":S\t" << loglik/repeats << "\t" << (end - start)*1000.0/CLOCKS_PER_SEC << endl;
    
    delete hmm_sparse;
    
    //
    // HMM with explicit dense mode
    //
    HMM * hmm_dense = new_hmm_instance(new InnerFwdDense<TransTableT *>(),
                                       new InnerBckDense<TransTableT *, EmissionTableT *>(),
                                       transitions,
                                       emissions,
                                       init_log_probs);
    loglik = 0;
    
    start = clock();
    for (int r = 0; r < repeats; ++r)
      loglik += hmm_dense->forward(iter, fwd);
    
    end = clock();
    cout << "T" << test_number << ":D\t" << loglik/repeats << "\t" << (end - start)*1000.0/CLOCKS_PER_SEC << endl;
    
    delete hmm_dense;
  }
}

int main(int argc, char ** argv) {
  int seq_len;
  int n_states;
  int sparseness;
  int repeats;

  /* data */
  double * data;
  double * auto_corr_covar;
  double * lambda_covar;
  double * covar_2d;
  double * emission_2d;

  /* params */
  double auto_corr = 0.2;
  double lambda = 0.2;

  if (argc != 5) {
    cout << "Usage: " << argv[0] << " <seq lengh> <n. states> <sparseness> <repeats>" << endl;
    return EXIT_FAILURE;
  }

  seq_len = atoi(argv[1]);
  n_states = atoi(argv[2]);
  sparseness = atoi(argv[3]);
  repeats = atoi(argv[4]);

  /* fill in dataset */
  data = new double[seq_len];
  fill_sequence(data, seq_len, 0, 3, 1);
  auto_corr_covar = new double[seq_len];
  fill_sequence(auto_corr_covar, seq_len, 0.1, 0.3, 0.1);
  lambda_covar = new double[seq_len];
  fill_sequence(lambda_covar, seq_len, 0.1, 0.3, 0.1);
  covar_2d = seq_merge(auto_corr_covar, lambda_covar, seq_len);
  emission_2d = seq_merge(data, data, seq_len);

  double * init_log_probs = new double[n_states];
  init_log_probs[0] = 0; /* log(1) */
  for (int i = 1; i < n_states; ++i)
    init_log_probs[i] = -std::numeric_limits<double>::infinity();

  double * fwd = new double[seq_len * n_states];
  
  /* Test 1: Poisson emission, AutoCorr transition
   *
   */
  int eslot_dim = 1;
  Iter * iter1 = new Iter(seq_len, 1, &eslot_dim, data, 0, NULL, NULL);

  /* create states */
  HomogeneousTransitions * trans = create_homogeneous_transitions(n_states, auto_corr, sparseness);
  Emissions * emissions = create_poisson_emissions(n_states, lambda);

  run_test(*iter1, trans, emissions, init_log_probs, fwd, repeats, 1);
  
  /* clean up */
  delete iter1;
  delete trans;
  delete emissions;

  /* Test 2: PoissonCovar emission, AutoCorr transition
   *
   */
  int cslot_dim = 1;
  Iter * iter2 = new Iter(seq_len, 1, &eslot_dim, data, 1, &cslot_dim, lambda_covar);
  
  /* create states */
  trans = create_homogeneous_transitions(n_states, auto_corr, sparseness);
  emissions = create_poisson_covar_emissions(n_states, 0);
  
  run_test(*iter2, trans, emissions, init_log_probs, fwd, repeats, 2);
  
  /* clean up */
  delete iter2;
  delete trans;
  delete emissions;
  
  /* Test 3: (Poisson, PoissonCovar emission), AutoCorr transition
   *
   */
  int eslot_dim2[2] = {1, 1};
  Iter * iter3 = new Iter(seq_len, 2, eslot_dim2, emission_2d, 1, &cslot_dim, lambda_covar);
  
  /* create states */
  trans = create_homogeneous_transitions(n_states, auto_corr, sparseness);
  MultiEmissions * m_emissions = create_poisson_both_emissions(n_states, lambda, 0);
  
  run_test(*iter3, trans, m_emissions, init_log_probs, fwd, repeats, 3);
  
  /* clean up */
  delete iter3;
  delete trans;
  delete m_emissions;

  /* Test 4: Poisson emission, AutoCorrCovar transition
   *
   */
  Iter * iter4 = new Iter(seq_len, 1, &eslot_dim, data, 1, &cslot_dim, auto_corr_covar);
  
  /* create states */
  NonHomogeneousTransitions * nh_trans = create_nonhomogeneous_transitions(n_states, 0, sparseness);
  emissions = create_poisson_emissions(n_states, lambda);
  
  run_test(*iter4, nh_trans, emissions, init_log_probs, fwd, repeats, 4, false);
  
  /* clean up */
  delete iter4;
  delete nh_trans;
  delete emissions;

  /* Test 5: PoissonCovar emission, AutoCorrCovar transition
   *
   */
  int cslot_dim2[2] = {1, 1};
  Iter * iter5 = new Iter(seq_len, 1, &eslot_dim, data, 2, cslot_dim2, covar_2d);
  
  /* create states */
  nh_trans = create_nonhomogeneous_transitions(n_states, 0, sparseness);
  emissions = create_poisson_covar_emissions(n_states, 1);
  
  run_test(*iter5, nh_trans, emissions, init_log_probs, fwd, repeats, 5, false);
  
  /* clean up */
  delete iter5;
  delete nh_trans;
  delete emissions;

  /* clean up */
  delete[] fwd;
  delete[] data;
  delete[] auto_corr_covar;
  delete[] lambda_covar;
  delete[] covar_2d;
  delete[] emission_2d;
  delete[] init_log_probs;

  return EXIT_SUCCESS;
}
