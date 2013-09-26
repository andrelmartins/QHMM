#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>

bool QHMM_isnan(const double x);
bool QHMM_isinf(const double x);

double QHMM_digamma(const double x);
double QHMM_trigamma(const double x);

double QHMM_logdiff(const double ln_x1, const double ln_x2);
double QHMM_logsum(const double ln_x1, const double ln_x2);
double QHMM_log_gamma(const double x);

double QHMM_log_gamma_cdf_lower(const double x, const double shape, const double scale);
double QHMM_log_gamma_cdf_upper(const double x, const double shape, const double scale);


#endif
