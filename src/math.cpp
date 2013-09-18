#include "math.hpp"
#include <cmath>

#if defined(USE_RMATH)
#include <Rmath.h>
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

double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale) {
  return pgamma(x, shape, scale, TRUE, TRUE);
}

double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale) {
  return pgamma(x, shape, scale, FALSE, TRUE);
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

double QHMM_log_gamma(const double x) {
  return gsl_sf_lngamma(x);
}

double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale) {
  return log(gsl_sf_gamma_inc_P(shape, x / scale));
}

double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale) {
  return log(gsl_sf_gamma_inc_Q(shape, x / scale));
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

#endif
