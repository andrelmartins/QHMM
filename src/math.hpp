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

typedef double qhmmfn(int n, double *x, void *params);

double QHMM_fminimizer(qhmmfn func, int n, double * x0, void * params, int maxit, double tol, int * out_fail);

void QHMM_rnd_prepare(void);
void QHMM_rnd_cleanup(void);
double QHMM_runif(void);

typedef double optimfn(int, double *, void *);
typedef void optimgr(int, double *, double *, void *);

void QHMM_lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
            double *Fmin, optimfn fn, optimgr gr, int *fail, void *ex,
            double factr, double pgtol, int *fncount, int *grcount,
            int maxit, char *msg, int trace, int nREPORT);
#endif
