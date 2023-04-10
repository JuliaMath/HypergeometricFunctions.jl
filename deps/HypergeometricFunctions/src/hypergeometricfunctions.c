#include "hypergeometricfunctions.h"

static inline void hf_update2(double (*N)[2], double (*D)[2], double (*T)[2], const double t[2]) {
    (*N)[0] = (*N)[1];
    (*N)[1] = t[0];
    (*D)[0] = (*D)[1];
    (*D)[1] = t[1];
    (*T)[0] = (*T)[1];
    (*T)[1] = (*N)[1]/(*D)[1];
}

static inline void hf_update3(double (*N)[3], double (*D)[3], double (*T)[3], const double t[2]) {
    (*N)[0] = (*N)[1];
    (*N)[1] = (*N)[2];
    (*N)[2] = t[0];
    (*D)[0] = (*D)[1];
    (*D)[1] = (*D)[2];
    (*D)[2] = t[1];
    (*T)[0] = (*T)[1];
    (*T)[1] = (*T)[2];
    (*T)[2] = (*N)[2]/(*D)[2];
}

double hf_drummond0f0(const double z, const int kmax) {
    if (fabs(z) < eps()) {return 1.0;}
    double zeta = 1.0/z;
    double N[2] = {zeta, (2.0*zeta-1.0)*zeta + 2.0*zeta};
    double D[2] = {zeta, (2.0*zeta-1.0)*zeta};
    double T[2] = {N[0]/D[0], N[1]/D[1]};
    int k = 1;
    double t[2] = {((k+2.0)*zeta-1.0)*N[1] + k*zeta*N[0] + zeta, ((k+2.0)*zeta-1.0)*D[1] + k*zeta*D[0]};
    hf_update2(&N, &D, &T, t);
    k++;
    while (k < kmax && errcheck(T[0], T[1], 8.0*eps())) {
        t[0] = ((k+2.0)*zeta-1.0)*N[1] + k*zeta*N[0];
        t[1] = ((k+2.0)*zeta-1.0)*D[1] + k*zeta*D[0];
        hf_update2(&N, &D, &T, t);
        k++;
    }
    return isfinite(T[1]) ? T[1] : T[0];
}

void hf_drummond0f0_mpfr(mpfr_t ret, const mpfr_t z, const int kmax, mpfr_prec_t prec, mpfr_rnd_t rnd) {
    mpfr_t epsilon, tol;
    mpfr_init2(epsilon, prec);
    mpfr_init2(tol, prec);
    eps_mpfr(epsilon, prec, rnd);
    mpfr_mul_d(tol, epsilon, 8.0, rnd);
    mpfr_t absz;
    mpfr_init2(absz, prec);
    mpfr_abs(absz, z, rnd);

    if (mpfr_less_p(absz, epsilon)) {
        mpfr_set_d(ret, 1.0, rnd);
        return;
    }
    //double zeta = 1.0/z;
    mpfr_t zeta;
    mpfr_init2(zeta, prec);
    mpfr_d_div(zeta, 1.0, z, rnd);

    mpfr_t N[2], D[2], T[2];
    mpfr_init2(N[0], prec);
    mpfr_init2(N[1], prec);
    mpfr_init2(D[0], prec);
    mpfr_init2(D[1], prec);
    mpfr_init2(T[0], prec);
    mpfr_init2(T[1], prec);

    mpfr_t t1, t2, t3, t4, t5, t6;
    mpfr_init2(t1, prec);
    mpfr_init2(t2, prec);
    mpfr_init2(t3, prec);
    mpfr_init2(t4, prec);
    mpfr_init2(t5, prec);
    mpfr_init2(t6, prec);

    // N[0] = zeta;
    // D[0] = zeta;
    mpfr_set(N[0], zeta, rnd);
    mpfr_set(D[0], zeta, rnd);
    // N[1] = (2.0*zeta-1.0)*zeta + 2.0*zeta
    // D[1] = (2.0*zeta-1.0)*zeta
    mpfr_set_d(t1, 2.0, rnd);
    mpfr_set_d(t2, 1.0, rnd);
    mpfr_fms(t3, t1, zeta, t2, rnd);
    mpfr_mul_d(t4, zeta, 2.0, rnd);
    mpfr_fma(N[1], t3, zeta, t4, rnd);
    mpfr_mul(D[1], t3, zeta, rnd);
    mpfr_div(T[0], N[0], D[0], rnd);
    mpfr_div(T[1], N[1], D[1], rnd);
    /*
    int k = 1;
    double t[2] = {((k+2.0)*zeta-1.0)*N[1] + k*zeta*N[0] + zeta, ((k+2.0)*zeta-1.0)*D[1] + k*zeta*D[0]};
    hf_update2(&N, &D, &T, t);
    k++;
    */
    int k = 1;
    mpfr_set_d(t1, k+2, rnd);
    mpfr_fms(t3, t1, zeta, t2, rnd);
    mpfr_mul_d(t4, zeta, k, rnd);
    mpfr_fma(t5, t4, N[0], zeta, rnd);
    mpfr_fma(t6, t3, N[1], t5, rnd);
    mpfr_set(N[0], N[1], rnd);
    mpfr_set(N[1], t6, rnd);
    mpfr_mul(t5, t4, D[0], rnd);
    mpfr_fma(t6, t3, D[1], t5, rnd);
    mpfr_set(D[0], D[1], rnd);
    mpfr_set(D[1], t6, rnd);
    mpfr_set(T[0], T[1], rnd);
    mpfr_div(T[1], N[1], D[1], rnd);
    k++;
    /*
    while (k < kmax && errcheck(T[0], T[1], 8.0*eps())) {
        t[0] = ((k+2.0)*zeta-1.0)*N[1] + k*zeta*N[0];
        t[1] = ((k+2.0)*zeta-1.0)*D[1] + k*zeta*D[0];
        hf_update2(&N, &D, &T, t);
        k++;
    }
    */
    mpfr_t errtemps[6];
    mpfr_init2(errtemps[0], prec);
    mpfr_init2(errtemps[1], prec);
    mpfr_init2(errtemps[2], prec);
    mpfr_init2(errtemps[3], prec);
    mpfr_init2(errtemps[4], prec);
    mpfr_init2(errtemps[5], prec);
    while (k < kmax && errcheck_mpfr(T[0], T[1], tol, errtemps, prec, rnd)) {
        mpfr_set_d(t1, k+2, rnd);
        mpfr_fms(t3, t1, zeta, t2, rnd);
        mpfr_mul_d(t4, zeta, k, rnd);
        mpfr_mul(t5, t4, N[0], rnd);
        mpfr_fma(t6, t3, N[1], t5, rnd);
        mpfr_set(N[0], N[1], rnd);
        mpfr_set(N[1], t6, rnd);
        mpfr_mul(t5, t4, D[0], rnd);
        mpfr_fma(t6, t3, D[1], t5, rnd);
        mpfr_set(D[0], D[1], rnd);
        mpfr_set(D[1], t6, rnd);
        mpfr_set(T[0], T[1], rnd);
        mpfr_div(T[1], N[1], D[1], rnd);
        k++;
    }
    isfinite_mpfr(T[1]) ? mpfr_set(ret, T[1], rnd) : mpfr_set(ret, T[0], rnd);
    mpfr_clear(epsilon);
    mpfr_clear(tol);
    mpfr_clear(absz);
    mpfr_clear(zeta);
    mpfr_clear(N[0]);
    mpfr_clear(N[1]);
    mpfr_clear(D[0]);
    mpfr_clear(D[1]);
    mpfr_clear(T[0]);
    mpfr_clear(T[1]);
    mpfr_clear(t1);
    mpfr_clear(t2);
    mpfr_clear(t3);
    mpfr_clear(t4);
    mpfr_clear(t5);
    mpfr_clear(t6);
    mpfr_clear(errtemps[0]);
    mpfr_clear(errtemps[1]);
    mpfr_clear(errtemps[2]);
    mpfr_clear(errtemps[3]);
    mpfr_clear(errtemps[4]);
    mpfr_clear(errtemps[5]);
    return;
}

double hf_weniger0f0(const double z, const int kmax) {
    if (fabs(z) < eps()) {return 1.0;}
    double zeta = 1.0/z;
    double N[2] = {zeta, (2.0*zeta-1.0)*zeta + 2.0*zeta};
    double D[2] = {zeta, (2.0*zeta-1.0)*zeta};
    double T[2] = {N[0]/D[0], N[1]/D[1]};
    double t[2];
    int k = 1;
    while (k < kmax && errcheck(T[0], T[1], 8.0*eps())) {
        t[0] = (4.0*k+2.0)*zeta*N[1] + N[0];
        t[1] = (4.0*k+2.0)*zeta*D[1] + D[0];
        hf_update2(&N, &D, &T, t);
        k++;
    }
    return isfinite(T[1]) ? T[1] : T[0];
}

double hf_weniger1f0(const double alpha, const double z, const int kmax) {
    double absalpha = fabs(alpha);
    if (fabs(z) < eps() || absalpha < absalpha*eps()) {return 1.0;}
    double zeta = 1.0/z;
    double N[2] = {zeta/alpha, (2.0*zeta-(alpha+1.0))*zeta/alpha + 2.0*zeta};
    double D[2] = {zeta/alpha, (2.0*zeta-(alpha+1.0))*zeta/alpha};
    double T[2] = {N[0]/D[0], N[1]/D[1]};
    if (fabs(alpha+1.0) < (absalpha+1.0)*eps()) {return T[1];}
    N[1] /= alpha+1.0;
    D[1] /= alpha+1.0;
    double t[2];
    int k = 1;
    while (k < kmax && errcheck(T[0], T[1], 8.0*eps())) {
        t[0] = (2.0*k+1.0)*(2.0*zeta-1.0)*N[1] - (k-alpha)*N[0];
        t[1] = (2.0*k+1.0)*(2.0*zeta-1.0)*D[1] - (k-alpha)*D[0];
        hf_update2(&N, &D, &T, t);
        if (fabs(alpha+k+1.0) < (absalpha+k+1.0)*eps()) {return T[1];}
        N[1] /= alpha+k+1.0;
        D[1] /= alpha+k+1.0;
        k++;
    }
    return isfinite(T[1]) ? T[1] : T[0];
}

double hf_weniger2f0(const double alpha, const double beta, const double z, const int kmax) {
    double absalpha = fabs(alpha);
    double absbeta = fabs(beta);
    if (fabs(z) < eps() || absalpha*absbeta < absalpha*absbeta*eps()) {return 1.0;}
    double zeta = 1.0/z;
    double a0 = (alpha+1.0)*(beta+1.0);
    double N[3];
    double D[3];
    double T[3];
    N[0] = zeta/(alpha*beta);
    D[0] = zeta/(alpha*beta);
    T[0] = N[0]/D[0];
    N[1] = (2.0*zeta-a0)*N[0] + 2.0*zeta;
    D[1] = (2.0*zeta-a0)*D[0];
    T[1] = N[1]/D[1];
    if (fabs(a0) < (absalpha+1.0)*(absbeta+1.0)*eps()) {return T[1];}
    N[1] /= a0;
    D[1] /= a0;
    double t[2];
    int k = 1;
    a0 = (alpha+2.0)*(beta+2.0);
    //a1 = T(2*(2-(2*α*β+α+β+1)))
    double t0 = 6.0*zeta-6.0+3.0*alpha*beta;
    double t1 = 6.0*zeta-2.0*(2.0*alpha*beta+alpha+beta-1.0);
    N[2] = t0*N[1] - t1*N[0] - 6.0*zeta;
    D[2] = t0*D[1] - t1*D[0];
    T[2] = N[2]/D[2];
    if (fabs(a0) < (absalpha+2.0)*(absbeta+2.0)*eps()) {return T[2];}
    N[2] /= a0;
    D[2] /= a0;
    k++;
    while (k < 3 || (k < kmax && errcheck(T[1], T[2], 8.0*eps()))) {
        a0 = (alpha+k+1.0)*(beta+k+1.0);
        t0 = (4.0*k+2.0)*zeta-(k*(alpha+beta+3.0*k)-(alpha+1.0)*(beta+1.0))*(2.0*k+1.0)/(2.0*k-1.0);
        t1 = (4.0*k+2.0)*zeta+k*(3.0*k-alpha-beta)-(alpha+1.0)*(beta+1.0);
        //a1 = T((2*k*k-(2*α*β+α+β+1)))*T(2k)/T(2k-1)
        double a2 = ((alpha+1.0-k)*(beta+1.0-k))*(2.0*k+1.0)/(2.0*k-1.0);
        t[0] = t0*N[2] - t1*N[1] - a2*N[0];
        t[1] = t0*D[2] - t1*D[1] - a2*D[0];
        hf_update3(&N, &D, &T, t);
        if (fabs(a0) < (absalpha+k+1.0)*(absbeta+k+1.0)*eps()) {return T[2];}
        N[2] /= a0;
        D[2] /= a0;
        k++;
    }
    return isfinite(T[2]) ? T[2] : isfinite(T[1]) ? T[1] : T[0];
}

double hf_weniger2f1(const double alpha, const double beta, const double gamma, const double z, const int kmax) {
    double absalpha = fabs(alpha);
    double absbeta = fabs(beta);
    if (fabs(z) < eps() || absalpha*absbeta < absalpha*absbeta*eps()) {return 1.0;}
    double zeta = 1.0/z;
    double a0 = (alpha+1.0)*(beta+1.0);
    double b0 = 2.0*(gamma+1.0);
    double N[3];
    double D[3];
    double T[3];
    N[0] = gamma*zeta/(alpha*beta);
    D[0] = gamma*zeta/(alpha*beta);
    T[0] = N[0]/D[0];
    N[1] = (b0*zeta-a0)*N[0] + b0*zeta;
    D[1] = (b0*zeta-a0)*D[0];
    T[1] = N[1]/D[1];
    if (fabs(a0) < (absalpha+1.0)*(absbeta+1.0)*eps()) {return T[1];}
    N[1] /= a0;
    D[1] /= a0;
    double t[2];
    int k = 1;
    a0 = (alpha+2.0)*(beta+2.0);
    //a1 = T(2*(2-(2*α*β+α+β+1)))
    b0 = 6.0*(gamma+2.0);
    double b1 = -6.0*gamma;
    double t0 = b0*zeta-6.0+3.0*alpha*beta;
    double t1 = b1*zeta+2.0*(2.0*alpha*beta+alpha+beta-1.0);
    N[2] = t0*N[1] + t1*N[0] + b1*zeta;
    D[2] = t0*D[1] + t1*D[0];
    T[2] = N[2]/D[2];
    if (fabs(a0) < (absalpha+2.0)*(absbeta+2.0)*eps()) {return T[2];}
    N[2] /= a0;
    D[2] /= a0;
    k++;
    while (k < 3 || (k < kmax && errcheck(T[1], T[2], 8.0*eps()))) {
        a0 = (alpha+k+1.0)*(beta+k+1.0);
        //a1 = T((2*k*k-(2*α*β+α+β+1)))*T(2k)/T(2k-1)
        double a2 = ((alpha+1.0-k)*(beta+1.0-k))*(2.0*k+1.0)/(2.0*k-1.0);
        b0 = (4.0*k+2.0)*(gamma+k+1.0);
        b1 = (4.0*k+2.0)*(k-gamma-1.0);
        t0 = b0*zeta-(k*(alpha+beta+3.0*k)-(alpha+1.0)*(beta+1.0))*(2.0*k+1.0)/(2.0*k-1.0);
        t1 = b1*zeta-k*(3.0*k-alpha-beta)+(alpha+1.0)*(beta+1.0);
        t[0] = t0*N[2] + t1*N[1] - a2*N[0];
        t[1] = t0*D[2] + t1*D[1] - a2*D[0];
        hf_update3(&N, &D, &T, t);
        if (fabs(a0) < (absalpha+k+1.0)*(absbeta+k+1.0)*eps()) {return T[2];}
        N[2] /= a0;
        D[2] /= a0;
        k++;
    }
    return isfinite(T[2]) ? T[2] : isfinite(T[1]) ? T[1] : T[0];
}
