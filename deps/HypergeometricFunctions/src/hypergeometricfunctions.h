#ifndef HYPERGEOMETRICFUNCTIONS_H
#define HYPERGEOMETRICFUNCTIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpfr.h>

#define MAX(a,b) ((a) > (b) ? a : b)
#define M_EPSf         0x1p-23f               /* powf(2.0f, -23)    */
#define M_EPS          0x1p-52                /* pow(2.0, -52)      */
#define M_EPSl         0x1p-64l               /* powl(2.0l, -64)    */

#define ZERO(FLT) ((FLT) 0)
#define ONE(FLT) ((FLT) 1)
#define TWO(FLT) ((FLT) 2)

static inline float epsf(void) {return M_EPSf;}
static inline double eps(void) {return M_EPS;}
static inline long double epsl(void) {return M_EPSl;}
static inline void eps_mpfr(mpfr_t epsilon, mpfr_prec_t prec, mpfr_rnd_t rnd) {
    mpfr_t nxt;
    mpfr_init2(nxt, prec);
    mpfr_set_d(nxt, 1.0, rnd);
    mpfr_nextabove(nxt);
    mpfr_sub_d(epsilon, nxt, 1.0, rnd);
    mpfr_clear(nxt);
    return;
}
static inline int errcheck(const double x, const double y, const double tol) {return (isfinite(x) && isfinite(y) && (fabs(x-y) > MAX(fabs(x), fabs(y))*tol));}
static inline int isfinite_mpfr(const mpfr_t x) {
    return !(mpfr_inf_p(x) != 0) && !(mpfr_nan_p(x) != 0);
}
static inline int errcheck_mpfr(const mpfr_t x, const mpfr_t y, const mpfr_t tol, mpfr_t t[6], mpfr_prec_t prec, mpfr_rnd_t rnd) {
    mpfr_abs(t[0], x, rnd);
    mpfr_abs(t[1], y, rnd);
    mpfr_sub(t[2], x, y, rnd);
    mpfr_abs(t[3], t[2], rnd);
    mpfr_max(t[4], t[0], t[1], rnd);
    mpfr_mul(t[5], t[4], tol, rnd);
    return isfinite_mpfr(x) && isfinite_mpfr(y) && mpfr_greater_p(t[3], t[5]);
}

double hf_drummond0f0(const double z, const int kmax);
void hf_drummond0f0_mpfr(mpfr_t ret, const mpfr_t z, const int kmax, mpfr_prec_t prec, mpfr_rnd_t rnd);

double hf_weniger0f0(const double z, const int kmax);
double hf_weniger1f0(const double alpha, const double z, const int kmax);
double hf_weniger2f0(const double alpha, const double beta, const double z, const int kmax);
double hf_weniger2f1(const double alpha, const double beta, const double gamma, const double z, const int kmax);

#endif // HYPERGEOMETRICFUNCTIONS_H
