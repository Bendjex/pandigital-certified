/*
 * motora.c - C3 certifier, motor A.
 *
 * Tables T(S) for |S| <= TBL (Theorem A) and window queries Q(S, W) (Lemma 8, Theorem B).
 * Values are nonnegative representatives |val E| (Lemma 1), powers have positive bases (Lemma 2),
 * the forms are those of Lemma 3.
 *
 * A class is either
 *   ordinary: key r (double) and code k; enclosure [max(0, pred(r - rho)), succ(r + rho)],
 *             rho = 2^(k-1) ulp(r) (k = 0: the exact value r); optionally the exact rational p/q;
 *   L-class:  a value x near 1 (x != 1) stored by an enclosure [llo, lhi] of lambda = ln x,
 *             which never contains 0; its enclosure for x is [exp_lo(llo), exp_hi(lhi)].
 * Derivations whose enclosure contains 0 (not the exact zero), contains 1 (not the exact one and
 * not representable as an L-class) or is relatively wider than 2^-20 are SINGULAR: they are written
 * out for motor B, whose exact answers (Lemma 7) are read back from a resolution file.
 *
 * Arithmetic: IEEE 754 binary64, round to nearest, no FMA contraction (-ffp-contract=off), and
 * CORE-MATH cr_pow (correctly rounded).  Every bound is pred()/succ() of a correctly rounded
 * result (Lemma 4).  exp and ln bounds come from cr_pow: e^z = 2^(z log2 e) with log2 e enclosed by
 * two doubles; a logarithm is guessed with a deterministic routine and the guess is accepted only
 * after it is validated with those exp bounds.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num(void) { return 0; }
static int omp_get_max_threads(void) { return 1; }
static double omp_get_wtime(void) { return 0; }
#endif

double cr_pow(double, double);

#define MAXD 9
#ifndef TBL
#define TBL 5
#endif
#define MAXDEPTH 8
typedef long long i64;
typedef __int128 i128;

static int ND, ORD, FULL;
static inline int popc(int m) { return __builtin_popcount(m); }
static inline int lowbit(int m) { return m & -m; }
static inline int is_interval(int m) { if (!m) return 0; int x = m / lowbit(m); return (x & (x + 1)) == 0; }
static inline int mask_ok(int m) { return !ORD || is_interval(m); }

/* ================================================================== doubles (Lemma 4) */
static inline double pred(double x) { return nextafter(x, -INFINITY); }
static inline double succ(double x) { return nextafter(x, INFINITY); }
static double RLO_SAFE, RHI_SAFE;                 /* pred(1e-300), succ(1e300): never discard a bound value */
static const double NEAR1 = 0x1p-10;              /* |ln x| below this: power results kept as L-classes */
static const double WIDE = 0x1p-20;               /* relative width above this: singular */
/* log2(e) = 1.44269504088896340735992468100189214...; RN(log2 e) = 0x1.71547652b82fep+0 lies below it
   (checked with mpmath in tests/t3_constants.py), so these two doubles enclose log2 e */
static const double L2E_LO = 0x1.71547652b82fep+0, L2E_HI = 0x1.71547652b82ffp+0;

static inline double ulp_of(double r) { int e; frexp(r, &e); return ldexp(1.0, e - 53); }

static inline void hull_of(double r, int code, double *lo, double *hi) {
    if (code == 0) { *lo = *hi = r; return; }
    double rho = ldexp(ulp_of(r), code - 1);
    double a = pred(r - rho), b = succ(r + rho);
    *lo = a < 0 ? 0 : a; *hi = b;
}
static inline int code_for(double r, double lo, double hi) {   /* lo <= r <= hi */
    double d1 = r - lo, d2 = hi - r;
    double d = succ(d1 > d2 ? d1 : d2);
    double q = succ(d / ulp_of(r));
    if (!(q > 1)) return 1;
    int e; double m = frexp(q, &e);
    int code = (m == 0.5) ? e : e + 1;
    return code > 254 ? 254 : code;
}

/* exp and ln bounds */
static inline double exp_hi(double z) {             /* >= e^z */
    if (z == -INFINITY) return 0;
    double w = succ(z >= 0 ? z * L2E_HI : z * L2E_LO);
    return succ(cr_pow(2.0, w));
}
static inline double exp_lo(double z) {             /* <= e^z */
    if (z == -INFINITY) return 0;
    double w = pred(z >= 0 ? z * L2E_LO : z * L2E_HI);
    double p = cr_pow(2.0, w);
    return p > 0 ? pred(p) : 0;
}
static double det_log(double x) {                   /* deterministic guess of ln x (IEEE ops only) */
    int e; double m = frexp(x, &e);
    if (m < 0.70710678118654752) { m *= 2; e -= 1; }
    double u = (m - 1) / (m + 1), u2 = u * u, s = 0, t = u;
    for (int k = 1; k < 41; k += 2) { s += t / k; t *= u2; }
    return 2 * s + e * 0.6931471805599453;
}
static double ln_lo(double x) {                     /* <= ln x  (x > 0) */
    if (x == INFINITY) return DBL_MAX;
    double g = det_log(x), step = fabs(g) * 0x1p-50 + 0x1p-1070;
    for (int it = 0; it < 2000; it++) { if (exp_hi(g) <= x) return g; g -= step; step *= 2; }
    return -INFINITY;
}
static double ln_hi(double x) {                     /* >= ln x  (x > 0) */
    if (x == INFINITY) return INFINITY;
    double g = det_log(x), step = fabs(g) * 0x1p-50 + 0x1p-1070;
    for (int it = 0; it < 2000; it++) { if (exp_lo(g) >= x) return g; g += step; step *= 2; }
    return INFINITY;
}
/* product / quotient of intervals, outward */
static inline void imul(double a1, double a2, double b1, double b2, double *lo, double *hi) {
    double p[4] = {a1 * b1, a1 * b2, a2 * b1, a2 * b2};
    double mn = p[0], mx = p[0];
    for (int i = 1; i < 4; i++) { if (p[i] < mn) mn = p[i]; if (p[i] > mx) mx = p[i]; }
    *lo = pred(mn); *hi = succ(mx);
}
static inline void idiv(double a1, double a2, double b1, double b2, double *lo, double *hi) {  /* 0 not in b */
    double p[4] = {a1 / b1, a1 / b2, a2 / b1, a2 / b2};
    double mn = p[0], mx = p[0];
    for (int i = 1; i < 4; i++) { if (p[i] < mn) mn = p[i]; if (p[i] > mx) mx = p[i]; }
    *lo = pred(mn); *hi = succ(mx);
}

/* ================================================================== exact rationals (Lemma 6) */
#define RLIM ((i64)1 << 53)
typedef struct { i64 p, q; } Rat;
static inline i128 gcd128(i128 a, i128 b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { i128 t = a % b; a = b; b = t; } return a; }
static inline int rat_set(i128 p, i128 q, Rat *o) {
    if (q == 0) return 0;
    if (q < 0) { p = -p; q = -q; }
    i128 g = gcd128(p, q);
    if (g > 1) { p /= g; q /= g; }
    if (p > RLIM || -p > RLIM || q > RLIM) return 0;
    o->p = (i64)p; o->q = (i64)q; return 1;
}
static inline int rat_add(Rat a, Rat b, Rat *o) { return rat_set((i128)a.p * b.q + (i128)b.p * a.q, (i128)a.q * b.q, o); }
static inline int rat_sub(Rat a, Rat b, Rat *o) { return rat_set((i128)a.p * b.q - (i128)b.p * a.q, (i128)a.q * b.q, o); }
static inline int rat_mul(Rat a, Rat b, Rat *o) { return rat_set((i128)a.p * b.p, (i128)a.q * b.q, o); }
static inline int rat_div(Rat a, Rat b, Rat *o) { if (b.p == 0) return 0; return rat_set((i128)a.p * b.q, (i128)a.q * b.p, o); }
static int ipow_lim(i64 b, i64 k, i64 *o) {
    i128 r = 1;
    for (i64 i = 0; i < k; i++) { r *= b; if (r > RLIM) return 0; if (b <= 1) break; }
    *o = (i64)r; return 1;
}
static int iroot(i64 x, int d, i64 *o) {
    if (x < 2) { *o = x; return 1; }
    double g = pow((double)x, 1.0 / d);           /* a guess only: checked exactly below */
    for (i64 c = (i64)g - 1; c <= (i64)g + 1; c++) {
        if (c < 0) continue;
        i128 r = 1; int ok = 1;
        for (int i = 0; i < d; i++) { r *= c; if (r > (i128)x) { ok = 0; break; } }
        if (ok && r == (i128)x) { *o = c; return 1; }
    }
    return 0;
}
static int rat_pow(Rat a, Rat e, Rat *o) {          /* a > 0 */
    if (a.p <= 0) return 0;
    if (e.q > 64 || e.p > 4096 || -e.p > 4096) return 0;
    i64 rp, rq;
    if (!iroot(a.p, (int)e.q, &rp) || !iroot(a.q, (int)e.q, &rq)) return 0;
    i64 k = e.p < 0 ? -e.p : e.p, np, nq;
    if (!ipow_lim(rp, k, &np) || !ipow_lim(rq, k, &nq)) return 0;
    return e.p >= 0 ? rat_set(np, nq, o) : rat_set(nq, np, o);
}
static inline double rat_double(Rat a) { return (double)a.p / (double)a.q; }
static inline int pow2(i64 q) { return q > 0 && (q & (q - 1)) == 0; }

