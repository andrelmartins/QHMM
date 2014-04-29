#include "math.hpp"
#include <cmath>

#if defined(USE_RMATH)
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Random.h>
#elif defined(USE_GSL)
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#endif

bool QHMM_isnan(const double x) {
  return std::isnan(x);
}

bool QHMM_isinf(const double x) {
  return std::isinf(x);
}

#if defined(USE_RMATH)
double QHMM_digamma(const double x) {
  return Rf_digamma(x);
}

double QHMM_trigamma(const double x) {
  return Rf_trigamma(x);
}

double QHMM_log_gamma(const double x) {
  return Rf_lgammafn(x);
}

double QHMM_logdiff(const double ln_x1, const double ln_x2) {
  return(ln_x1 + log1p(-exp(ln_x2 - ln_x1)));
}

double QHMM_logsum(const double ln_x1, const double ln_x2) {
  if (ln_x1 > ln_x2)
    return(ln_x1 + log1p(exp(ln_x2 - ln_x1)));
  return(ln_x2 + log1p(exp(ln_x1 - ln_x2)));
}

double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale) {
  return pgamma(x, shape, scale, TRUE, TRUE);
}

double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale) {
  return pgamma(x, shape, scale, FALSE, TRUE);
}

double QHMM_fminimizer(qhmmfn func, int n, double * x0, void * params, int maxit, double tol, int * out_fail) {
  double * xout = new double[n];
  double fout;
  double intol = 1e-8;
  double nm_alpha = 1;
  double nm_beta = 0.5;
  double nm_gamma = 2;
  int fncount = 0;

  nmmin(n, x0, xout, &fout, func, out_fail, tol, intol, params,
        nm_alpha, nm_beta, nm_gamma, 0,
        &fncount, maxit);
  
  for (int i = 0; i < n; ++i)
    x0[i] = xout[i];
  
  return fout;
}

void QHMM_rnd_prepare(void) {
  GetRNGstate();
}

void QHMM_rnd_cleanup(void) {
  PutRNGstate();
}

double QHMM_runif(void) {
  return unif_rand();
}

#elif defined(USE_GSL)

double QHMM_digamma(const double x) {
  return gsl_sf_psi(x):
}

double QHMM_trigamma(const double x) {
  return gsl_sf_psi_1(x);
}

double QHMM_logdiff(const double ln_x1, const double ln_x2) {
  return(ln_x1 + gsl_log1p(-exp(ln_x2 - ln_x1)));
}

double QHMM_logsum(const double ln_x1, const double ln_x2) {
  if (ln_x1 > ln_x2)
    return(ln_x1 + gsl_log1p(exp(ln_x2 - ln_x1)));
  return(ln_x2 + gsl_log1p(exp(ln_x1 - ln_x2)));
}

double QHMM_log_gamma(const double x) {
  return gsl_sf_lngamma(x);
}

double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale) {
  return log(gsl_sf_gamma_inc_P(shape, x / scale));
}

double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale) {
  return log(gsl_sf_gamma_inc_Q(shape, x / scale));
}

double QHMM_fminimizer(qhmmfn func, int n, double * x0, void * params, int maxit, double tol, int * out_fail) {

  // TODO: Implement GSL version ...
  /* Implement this using GSL's nmsimplex2

     - intermediate function that maps qhmmfn to the function form needed by GSL
     - main optimization loop ...
  */
}

void QHMM_rnd_prepare(void) {
  // TODO: seed random number generator
  // TODO: initialize GSL random number generator
}

void QHMM_rnd_cleanup(void) {
  // TODO: free GSL random number generator
}

double QHMM_runif(void) {
  // TODO: Implement GSL version (use gsl's random number generators ...)
  // TODO: fix this! Should be number in [0, 1] */
  return drand48(); /* random number in [0, 1) */
}

#else

/* Thanks to Charles Danko for the di/tri gamma code */

double QHMM_digamma(const double k) {
  if(k < 8) {
    return(QHMM_digamma(k+1)-1/k);
  }
  else {
    return(log(k)-(1+(1- (0.1-1/(21*k*k)) /(k*k))/(6*k))/(2*k));
  }
}

double QHMM_trigamma(const double k) {
  if(k < 8) {
    return(QHMM_trigamma(k+1)+1/(k*k));
  }
  else {
    return((1+(1+(1-(1/5-1/(7*k*k))/(k*k))/(3*k))/(2*k))/k);
  }
}

double QHMM_logdiff(const double ln_x1, const double ln_x2) {
  return(ln_x1 + log(1 - exp(ln_x2 - ln_x1)));
}

double QHMM_logsum(const double ln_x1, const double ln_x2) {
  if (ln_x1 > ln_x2)
    return(ln_x1 + log(1 + exp(ln_x2 - ln_x1)));
  return(ln_x2 + log(1 + exp(ln_x1 - ln_x2)));
}


double QHMM_log_gamma(const double x) {
  return std::lgamma(x);
}

/*
double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale) {
1/\Gamma(a) \int_0^x dt t^{a-1} \exp(-t) for a > 0, x >= 0.

  return log(gsl_sf_gamma_inc_P(shape, x / scale));
}

double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale) {
  return log(gsl_sf_gamma_inc_Q(shape, x / scale));
}
*/

void QHMM_rnd_prepare(void) {
  // TODO: seed random number generator
}

void QHMM_rnd_cleanup(void) {
  // blank
}

double QHMM_runif(void) {
  // TODO: fix this! Should be number in [0, 1] */
  return drand48(); /* random number in [0, 1) */
}

#endif
