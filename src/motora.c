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
/* pred(x) = nextafter(x, -inf), succ(x) = nextafter(x, +inf), written with the bit pattern (same results for
   every input, including zeros, subnormals, infinities and NaN; checked against nextafter in tests) */
static inline double pred(double x) {
    if (x != x || x == -INFINITY) return x;
    if (x == 0) return -0x1p-1074;
    uint64_t u; memcpy(&u, &x, 8); u = x > 0 ? u - 1 : u + 1; memcpy(&x, &u, 8); return x;
}
static inline double succ(double x) {
    if (x != x || x == INFINITY) return x;
    if (x == 0) return 0x1p-1074;
    uint64_t u; memcpy(&u, &x, 8); u = x > 0 ? u + 1 : u - 1; memcpy(&x, &u, 8); return x;
}
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


/* ================================================================== anchored values
   An anchored class is a value  x = q e^mu  with q an exact rational (|p|, q <= 2^53) and
   0 < |mu| <= MU_MAX, mu known by an enclosure [m1, m2] that excludes 0.  Values that sit within a
   double rounding of a simple rational (1 + 10^-80, 2*9^(10^-60), .2 + .1^81, ...) are kept this way,
   so that they stay distinguishable from the rational and from each other. */
#define MU_MAX 0x1p-10
#define TINY_REL 0x1p-12

static inline void rat_hull(Rat q, double *lo, double *hi) {
    double r = rat_double(q);
    if (q.p == 0 || pow2(q.q)) { *lo = *hi = r; } else { *lo = pred(r); *hi = succ(r); }
}
static inline int rat_cmp(Rat a, Rat b) { i128 l = (i128)a.p * b.q, r = (i128)b.p * a.q; return (l > r) - (l < r); }
/* e^m - 1 for m in [m1, m2], |m| <= 2^-10:   m <= e^m - 1 <= m + m^2 */
static inline void expm1_hull(double m1, double m2, double *e1, double *e2) { *e1 = m1; *e2 = succ(m2 + succ(m2 * m2)); }
/* ln(1+u) for u in [u1, u2], |u| <= 2^-9:   u - u^2 <= ln(1+u) <= u */
static inline void log1p_hull(double u1, double u2, double *l1, double *l2) { *l1 = pred(u1 - succ(u1 * u1)); *l2 = u2; }
static inline void isub(double a1, double a2, double b1, double b2, double *lo, double *hi) { *lo = pred(a1 - b2); *hi = succ(a2 - b1); }
static inline void iadd(double a1, double a2, double b1, double b2, double *lo, double *hi) { *lo = pred(a1 + b1); *hi = succ(a2 + b2); }

/* ================================================================== tables */
typedef struct {
    i64 n;                       /* ordinary and rational classes, sorted by key */
    double *v; unsigned char *code;
    i64 *rp, *rq;                /* operand tables: rational (rq == 0: none) */
    double *lo, *hi, *lnlo, *lnhi;
    int maxcode;
    i64 na;                      /* anchored classes, sorted by (anchor, m1, m2) */
    i64 *ap, *aq; double *am1, *am2, *axlo, *axhi;
    double amaxw;
} Tab;
static Tab T[1 << MAXD];
static Tab Lit[1 << MAXD];

typedef struct { double r; unsigned char code, rat, anc; i64 p, q; double m1, m2; } Der;
static int cmp_der(const void *A, const void *B) {
    const Der *a = A, *b = B;
    if (a->anc != b->anc) return a->anc < b->anc ? -1 : 1;
    if (a->r != b->r) return a->r < b->r ? -1 : 1;
    if (a->anc) {
        if (a->p != b->p) return a->p < b->p ? -1 : 1;
        if (a->q != b->q) return a->q < b->q ? -1 : 1;
        if (a->m1 != b->m1) return a->m1 < b->m1 ? -1 : 1;
        if (a->m2 != b->m2) return a->m2 < b->m2 ? -1 : 1;
        return 0;
    }
    if (a->rat != b->rat) return a->rat > b->rat ? -1 : 1;
    if (a->rat) {
        if (a->p != b->p) return a->p < b->p ? -1 : 1;
        if (a->q != b->q) return a->q < b->q ? -1 : 1;
    }
    return 0;
}

enum { K_ORD, K_RAT, K_ANC };
typedef struct { int kind, zero, one, code; double v, lo, hi, lnlo, lnhi, m1, m2; Rat r; } Opd;

static inline void anc_ln(Rat q, double m1, double m2, double *l1, double *l2) {
    if (q.p == q.q) { *l1 = m1; *l2 = m2; return; }
    double ql, qh; rat_hull(q, &ql, &qh);
    *l1 = pred(ln_lo(ql) + m1); *l2 = succ(ln_hi(qh) + m2);
}
static inline void anc_x(Rat q, double m1, double m2, double *x1, double *x2) {
    double ql, qh; rat_hull(q, &ql, &qh);
    *x1 = pred(ql * exp_lo(m1)); *x2 = succ(qh * exp_hi(m2));
}
static inline Opd opd_ord(int m, i64 i) {
    Opd o; Tab *t = &T[m];
    memset(&o, 0, sizeof o);
    o.v = t->v[i]; o.code = t->code[i]; o.lo = t->lo[i]; o.hi = t->hi[i];
    o.lnlo = t->lnlo[i]; o.lnhi = t->lnhi[i];
    if (t->rq[i]) { o.kind = K_RAT; o.r.p = t->rp[i]; o.r.q = t->rq[i]; o.zero = o.r.p == 0; o.one = o.r.p == 1 && o.r.q == 1; }
    else o.kind = K_ORD;
    return o;
}
static inline Opd opd_anc(int m, i64 i) {
    Opd o; Tab *t = &T[m];
    memset(&o, 0, sizeof o);
    o.kind = K_ANC; o.r.p = t->ap[i]; o.r.q = t->aq[i]; o.m1 = t->am1[i]; o.m2 = t->am2[i];
    o.lo = t->axlo[i]; o.hi = t->axhi[i];
    anc_ln(o.r, o.m1, o.m2, &o.lnlo, &o.lnhi);
    o.v = o.lo + (o.hi - o.lo) * 0.5; if (o.v < o.lo) o.v = o.lo; if (o.v > o.hi) o.v = o.hi;
    o.code = 255;
    return o;
}

enum { F_ADD, F_SUB, F_MUL, F_DXY, F_DYX, F_PXY, F_PXYN, F_PYX, F_PYXN, NFORM };
static const char *FNAME[NFORM] = {"x+y", "|x-y|", "x*y", "x/y", "y/x", "x^y", "x^-y", "y^x", "y^-x"};
enum { D_SKIP, D_OK, D_SING_Z, D_SING_O, D_SING_W, D_OUT };

static inline int finish_real(double r, double lo, double hi, Der *d) {
    if (hi < RLO_SAFE) return D_OUT;
    if (lo > RHI_SAFE) return D_OUT;
    if (!(lo > 0)) return D_SING_Z;
    if (!isfinite(r) || !(r > 0)) return D_SING_W;
    if (r < lo) r = lo;
    if (r > hi) r = hi;
    int code = code_for(r, lo, hi);
    double slo, shi;
    hull_of(r, code, &slo, &shi);
    if (!(slo > 0)) return D_SING_Z;
    if (slo <= 1 && shi >= 1) return D_SING_O;
    if (shi - slo > ldexp(r, -20)) return D_SING_W;
    memset(d, 0, sizeof *d);
    d->r = r; d->code = (unsigned char)code;
    return D_OK;
}
static inline int finish_rat(Rat q, Der *d) {
    if (q.p < 0) q.p = -q.p;
    memset(d, 0, sizeof *d);
    d->rat = 1; d->p = q.p; d->q = q.q;
    d->r = rat_double(q);
    d->code = (q.p == 0 || pow2(q.q)) ? 0 : 1;
    return D_OK;
}
/* q e^mu, mu in [n1, n2] */
static inline int make_anc(Rat Q, double n1, double n2, Der *d) {
    if (Q.p < 0) Q.p = -Q.p;
    if (!(n1 > 0) && !(n2 < 0)) return D_SING_O;             /* possibly exactly Q */
    if (n2 < MU_MAX && n1 > -MU_MAX) {
        memset(d, 0, sizeof *d);
        d->anc = 1; d->p = Q.p; d->q = Q.q; d->m1 = n1; d->m2 = n2; d->r = rat_double(Q); d->code = 255;
        return D_OK;
    }
    double x1, x2; anc_x(Q, n1, n2, &x1, &x2);
    return finish_real(x1 + (x2 - x1) * 0.5, x1, x2, d);
}
/* a value with ln-enclosure [l1, l2]: anchored at 1 if near 1, else ordinary */
static inline int finish_from_ln(double l1, double l2, Der *d) {
    if (l2 < MU_MAX && l1 > -MU_MAX) { Rat one = {1, 1}; return make_anc(one, l1, l2, d); }
    double lo = exp_lo(l1), hi = exp_hi(l2);
    double r = cr_pow(2.0, (l1 + (l2 - l1) * 0.5) * L2E_LO);
    return finish_real(r, lo, hi, d);
}
/* |Q + delta| with Q exact (any sign) and delta in [d1, d2] */
static int anc_affine(Rat Q, double d1, double d2, Der *d) {
    if (Q.p == 0) {
        if (d1 > 0) return finish_real(d1 + (d2 - d1) * 0.5, d1, d2, d);
        if (d2 < 0) return finish_real(-d2 + (d2 - d1) * 0.5, -d2, -d1, d);
        return D_SING_Z;
    }
    if (Q.p < 0) { Q.p = -Q.p; double t = d1; d1 = -d2; d2 = -t; }
    double ql, qh; rat_hull(Q, &ql, &qh);
    if (d2 <= ql * TINY_REL && d1 >= -ql * TINY_REL) {
        double u1, u2, n1, n2;
        idiv(d1, d2, ql, qh, &u1, &u2);
        log1p_hull(u1, u2, &n1, &n2);
        return make_anc(Q, n1, n2, d);
    }
    double lo = pred(ql + d1), hi = succ(qh + d2);          /* Q + delta in [lo, hi]; the value is |Q + delta| */
    if (lo > 0) return finish_real(lo + (hi - lo) * 0.5, lo, hi, d);
    if (hi < 0) return finish_real(-hi + (hi - lo) * 0.5, -hi, -lo, d);
    return D_SING_Z;
}
/* additive decomposition  x = Q + delta : rational (Q, 0), anchored (q, q (e^mu - 1)), other (0, x) */
static inline void decomp(const Opd *x, Rat *Q, double *d1, double *d2) {
    if (x->kind == K_RAT) { *Q = x->r; *d1 = *d2 = 0; return; }
    if (x->kind == K_ANC) {
        double e1, e2, ql, qh; expm1_hull(x->m1, x->m2, &e1, &e2); rat_hull(x->r, &ql, &qh);
        *Q = x->r; imul(ql, qh, e1, e2, d1, d2); return;
    }
    Q->p = 0; Q->q = 1; *d1 = x->lo; *d2 = x->hi;
}
static inline int copy_opd(const Opd *x, Der *d) {
    if (x->kind == K_ANC) return make_anc(x->r, x->m1, x->m2, d);
    if (x->kind == K_RAT) return finish_rat(x->r, d);
    memset(d, 0, sizeof *d); d->r = x->v; d->code = (unsigned char)x->code;
    return D_OK;
}
static inline void pow_hull(double xlo, double xhi, double elo, double ehi, double *lo, double *hi) {
    double c1, c2;   /* monotone on the box: x not straddling 1, e not straddling 0 (Lemma 5) */
    if (xlo > 1) {
        if (elo > 0) { c1 = cr_pow(xlo, elo); c2 = cr_pow(xhi, ehi); }
        else         { c1 = cr_pow(xhi, elo); c2 = cr_pow(xlo, ehi); }
    } else {
        if (elo > 0) { c1 = cr_pow(xlo, ehi); c2 = cr_pow(xhi, elo); }
        else         { c1 = cr_pow(xhi, ehi); c2 = cr_pow(xlo, elo); }
    }
    *lo = pred(c1); *hi = succ(c2);
}

static int derive(int form, const Opd *x, const Opd *y, Der *d) {
    double lo, hi, r;
    Rat q;
    int anc_x = x->kind != K_ORD, anc_y = y->kind != K_ORD;
    switch (form) {
    case F_ADD: case F_SUB: {
        if (x->zero) return copy_opd(y, d);
        if (y->zero) return copy_opd(x, d);
        if (x->kind == K_RAT && y->kind == K_RAT) {
            if (form == F_ADD ? rat_add(x->r, y->r, &q) : rat_sub(x->r, y->r, &q)) return finish_rat(q, d);
        } else if (anc_x || anc_y) {
            Rat Qx, Qy, Q; double a1, a2, b1, b2, d1, d2;
            decomp(x, &Qx, &a1, &a2); decomp(y, &Qy, &b1, &b2);
            if (form == F_ADD ? rat_add(Qx, Qy, &Q) : rat_sub(Qx, Qy, &Q)) {
                if (form == F_ADD) iadd(a1, a2, b1, b2, &d1, &d2); else isub(a1, a2, b1, b2, &d1, &d2);
                return anc_affine(Q, d1, d2, d);
            }
        }
        if (form == F_ADD) { r = x->v + y->v; lo = pred(x->lo + y->lo); hi = succ(x->hi + y->hi); return finish_real(r, lo, hi, d); }
        r = fabs(x->v - y->v);
        lo = pred(x->lo - y->hi); hi = succ(x->hi - y->lo);
        if (lo > 0) return finish_real(r, lo, hi, d);
        if (hi < 0) return finish_real(r, -hi, -lo, d);
        return D_SING_Z;
    }
    case F_MUL:
        if (x->zero || y->zero) { Rat z = {0, 1}; return finish_rat(z, d); }
        if (x->one) return copy_opd(y, d);
        if (y->one) return copy_opd(x, d);
        if (x->kind == K_RAT && y->kind == K_RAT && rat_mul(x->r, y->r, &q)) return finish_rat(q, d);
        if (anc_x && anc_y && rat_mul(x->r, y->r, &q)) {
            double m1 = 0, m2 = 0;
            if (x->kind == K_ANC) { m1 = x->m1; m2 = x->m2; }
            if (y->kind == K_ANC) { if (x->kind == K_ANC) iadd(m1, m2, y->m1, y->m2, &m1, &m2); else { m1 = y->m1; m2 = y->m2; } }
            return make_anc(q, m1, m2, d);
        }
        r = x->v * y->v; lo = pred(x->lo * y->lo); hi = succ(x->hi * y->hi);
        return finish_real(r, lo, hi, d);
    case F_DXY: case F_DYX: {
        const Opd *a = form == F_DXY ? x : y, *b = form == F_DXY ? y : x;
        if (b->zero) return D_SKIP;
        if (a->zero) { Rat z = {0, 1}; return finish_rat(z, d); }
        if (b->one) return copy_opd(a, d);
        if (a->kind == K_RAT && b->kind == K_RAT && rat_div(a->r, b->r, &q)) return finish_rat(q, d);
        if (a->kind != K_ORD && b->kind != K_ORD && rat_div(a->r, b->r, &q)) {
            double m1 = 0, m2 = 0;
            if (a->kind == K_ANC) { m1 = a->m1; m2 = a->m2; }
            if (b->kind == K_ANC) { double t1, t2; isub(m1, m2, b->m1, b->m2, &t1, &t2); m1 = t1; m2 = t2; }
            return make_anc(q, m1, m2, d);
        }
        if (!(b->lo > 0)) return D_SING_Z;
        r = a->v / b->v; lo = pred(a->lo / b->hi); hi = succ(a->hi / b->lo);
        return finish_real(r, lo, hi, d);
    }
    default: {
        const Opd *a = (form == F_PXY || form == F_PXYN) ? x : y;     /* base */
        const Opd *b = (form == F_PXY || form == F_PXYN) ? y : x;     /* exponent, sign below */
        int negx = (form == F_PXYN || form == F_PYXN);
        if (b->zero) {
            if (a->zero || negx) return D_SKIP;          /* 0^0 undefined; a^-0 = a^0 counted once */
            Rat one = {1, 1}; return finish_rat(one, d);
        }
        if (a->zero) {
            if (negx) return D_SKIP;
            Rat z = {0, 1}; return finish_rat(z, d);
        }
        if (a->one) { Rat one = {1, 1}; return finish_rat(one, d); }
        if (b->one && !negx) return copy_opd(a, d);
        Rat e = b->r; if (negx) e.p = -e.p;
        if (a->kind == K_RAT && b->kind == K_RAT && rat_pow(a->r, e, &q)) return finish_rat(q, d);
        if (!(b->lo > 0)) return D_SING_Z;
        double elo = negx ? -b->hi : b->lo, ehi = negx ? -b->lo : b->hi;
        if (a->kind == K_ANC) {                          /* (q e^mu)^e */
            if (b->kind == K_RAT && rat_pow(a->r, e, &q)) {
                double m1, m2; imul(elo, ehi, a->m1, a->m2, &m1, &m2);
                return make_anc(q, m1, m2, d);
            }
            double l1, l2; imul(a->lnlo, a->lnhi, elo, ehi, &l1, &l2);
            return finish_from_ln(l1, l2, d);
        }
        if (a->kind == K_RAT && b->kind == K_ANC) {      /* r^(q e^mu) = r^q e^(q (e^mu - 1) ln r) */
            Rat qq = b->r; if (negx) qq.p = -qq.p;
            if (rat_pow(a->r, qq, &q)) {
                double e1, e2, k1, k2, n1, n2, ql, qh;
                expm1_hull(b->m1, b->m2, &e1, &e2);
                rat_hull(qq, &ql, &qh);
                imul(ql, qh, e1, e2, &k1, &k2);
                imul(k1, k2, a->lnlo, a->lnhi, &n1, &n2);
                return make_anc(q, n1, n2, d);
            }
        }
        if (!(a->lo > 0)) return D_SING_Z;
        if (a->lo <= 1 && a->hi >= 1) return D_SING_O;
        r = cr_pow(a->v, negx ? -b->v : b->v);
        pow_hull(a->lo, a->hi, elo, ehi, &lo, &hi);
        if (hi < RLO_SAFE || lo > RHI_SAFE) return D_OUT;
        if ((lo <= 1 && hi >= 1) || (lo > 1 - MU_MAX && hi < 1 + MU_MAX)) {
            double l1, l2;                               /* near 1: ln x = e ln a */
            imul(a->lnlo, a->lnhi, elo, ehi, &l1, &l2);
            return finish_from_ln(l1, l2, d);
        }
        return finish_real(r, lo, hi, d);
    }
    }
}

/* ================================================================== literals */
static void lit_rec(int *dg, int k, int i, Der **buf, i64 *n, i64 *cap) {
    if (i == k) {
        i64 v = 0;
        for (int j = 0; j < k; j++) v = v * 10 + dg[j];
        i64 p10 = 1;
        for (int j = 0; j <= k; j++) {
            if (*n >= *cap) { *cap = *cap * 2 + 64; *buf = realloc(*buf, *cap * sizeof(Der)); }
            Rat q = {0, 1};
            rat_set(v, p10, &q);
            finish_rat(q, &(*buf)[(*n)++]);
            p10 *= 10;
        }
        return;
    }
    for (int j = i; j < k; j++) {
        int t = dg[i]; dg[i] = dg[j]; dg[j] = t;
        lit_rec(dg, k, i + 1, buf, n, cap);
        t = dg[i]; dg[i] = dg[j]; dg[j] = t;
    }
}

static Tab merge(Der *buf, i64 n, int keep) {
    qsort(buf, n, sizeof(Der), cmp_der);
    i64 m = 0;
    for (i64 i = 0; i < n; i++) {
        if (m > 0) {
            Der *p = &buf[m - 1];
            if (!p->anc && !buf[i].anc && p->r == buf[i].r && p->rat == buf[i].rat &&
                (!p->rat || (p->p == buf[i].p && p->q == buf[i].q))) {
                if (buf[i].code > p->code) p->code = buf[i].code;
                continue;
            }
            if (p->anc && buf[i].anc && cmp_der(p, &buf[i]) == 0) continue;
        }
        buf[m++] = buf[i];
    }
    i64 na = 0;
    for (i64 i = 0; i < m; i++) na += buf[i].anc;
    i64 no = m - na;                                   /* ordinary first (cmp_der) */
    Tab t; memset(&t, 0, sizeof t);
    t.n = no; t.na = na;
    t.v = malloc((no + 1) * sizeof(double)); t.code = malloc(no + 1);
    t.maxcode = 0;
    for (i64 i = 0; i < no; i++) {
        t.v[i] = buf[i].r; t.code[i] = buf[i].code;
        if (buf[i].code > t.maxcode) t.maxcode = buf[i].code;
    }
    t.ap = malloc((na + 1) * sizeof(i64)); t.aq = malloc((na + 1) * sizeof(i64));
    t.am1 = malloc((na + 1) * sizeof(double)); t.am2 = malloc((na + 1) * sizeof(double));
    t.axlo = malloc((na + 1) * sizeof(double)); t.axhi = malloc((na + 1) * sizeof(double));
    t.amaxw = 0;
    for (i64 i = 0; i < na; i++) {
        Der *b = &buf[no + i];
        t.ap[i] = b->p; t.aq[i] = b->q; t.am1[i] = b->m1; t.am2[i] = b->m2;
        Rat a = {b->p, b->q};
        anc_x(a, b->m1, b->m2, &t.axlo[i], &t.axhi[i]);
        double w = succ(b->m2 - b->m1);
        if (w > t.amaxw) t.amaxw = w;
    }
    if (keep) {
        t.rp = malloc((no + 1) * sizeof(i64)); t.rq = malloc((no + 1) * sizeof(i64));
        t.lo = malloc((no + 1) * sizeof(double)); t.hi = malloc((no + 1) * sizeof(double));
        t.lnlo = malloc((no + 1) * sizeof(double)); t.lnhi = malloc((no + 1) * sizeof(double));
        for (i64 i = 0; i < no; i++) {
            t.rp[i] = buf[i].rat ? buf[i].p : 0; t.rq[i] = buf[i].rat ? buf[i].q : 0;
            hull_of(t.v[i], t.code[i], &t.lo[i], &t.hi[i]);
            if (t.lo[i] > 0) { t.lnlo[i] = ln_lo(t.lo[i]); t.lnhi[i] = ln_hi(t.hi[i]); }
            else { t.lnlo[i] = -INFINITY; t.lnhi[i] = -INFINITY; }
            if (buf[i].rat && buf[i].p > 0) {               /* exact rational near 1: ln(p/q) = log1p((p-q)/q), Lemma 8 */
                Rat u = {buf[i].p - buf[i].q, buf[i].q};
                double u1, u2; rat_hull(u, &u1, &u2);
                if (u1 > -0x1p-9 && u2 < 0x1p-9 && u.p != 0) {
                    double l1, l2; log1p_hull(u1, u2, &l1, &l2);
                    if (l1 > t.lnlo[i]) t.lnlo[i] = l1;
                    if (l2 < t.lnhi[i]) t.lnhi[i] = l2;
                }
            }
        }
    }
    return t;
}

static Tab make_literals(int mask, int keep) {
    int dg[MAXD], k = 0;
    for (int b = 0; b < ND; b++) if (mask >> b & 1) dg[k++] = b + 1;
    i64 n = 0, cap = 256; Der *buf = malloc(cap * sizeof(Der));
    lit_rec(dg, k, ORD ? k : 0, &buf, &n, &cap);
    Tab t = merge(buf, n, keep);
    free(buf);
    return t;
}

/* ================================================================== singular derivations and resolutions */
typedef struct { int m1, m2, form, kind; Opd x, y; } Sing;
typedef struct { Sing *s; i64 n, cap; } SingList;
static SingList SL[1 << MAXD];

static void opd_key(const Opd *o, char *buf) {
    if (o->kind == K_ANC) {
        if (o->r.p == 1 && o->r.q == 1) sprintf(buf, "l:%a:%a", o->m1, o->m2);
        else sprintf(buf, "a:%lld/%lld:%a:%a", (long long)o->r.p, (long long)o->r.q, o->m1, o->m2);
    } else sprintf(buf, "o:%a:%d", o->v, o->code);
}
static void sing_key(int S, const Sing *s, char *buf) {
    char a[160], b[160];
    opd_key(&s->x, a); opd_key(&s->y, b);
    sprintf(buf, "%d %d %d %s %d %s", S, s->form, s->m1, a, s->m2, b);
}
typedef struct { char *key; char *items; } Res;
static Res *RES; static i64 NRES;
static int cmp_res(const void *a, const void *b) { return strcmp(((const Res *)a)->key, ((const Res *)b)->key); }
/* portable getline: the whole line (with its newline) in a growing buffer; returns its length, 0 at EOF */
static long read_line(char **buf, size_t *cap, FILE *f) {
    long n = 0;
    if (!*buf || *cap < 256) { *cap = 4096; *buf = realloc(*buf, *cap); }
    for (;;) {
        if (!fgets(*buf + n, (int)(*cap - n), f)) return n;
        n += (long)strlen(*buf + n);
        if (n > 0 && (*buf)[n - 1] == '\n') return n;
        if ((size_t)n + 1 >= *cap) { *cap *= 2; *buf = realloc(*buf, *cap); }
    }
}
static void load_res(const char *fn) {
    NRES = 0; RES = NULL;
    if (!fn) return;
    FILE *f = fopen(fn, "r");
    if (!f) return;
    i64 cap = 0; char *line = NULL; size_t lcap = 0; long len;
    while ((len = read_line(&line, &lcap, f)) > 0) {
        if (line[len - 1] == '\n') line[len - 1] = 0;
        char *bar = strstr(line, " |");
        if (!bar) continue;
        *bar = 0;
        if (NRES >= cap) { cap = cap * 2 + 1024; RES = realloc(RES, cap * sizeof(Res)); }
        RES[NRES].key = strdup(line); RES[NRES].items = strdup(bar[2] ? bar + 3 : ""); NRES++;
    }
    free(line); fclose(f);
    qsort(RES, NRES, sizeof(Res), cmp_res);
}
static const char *find_res(const char *key) {
    i64 lo = 0, hi = NRES - 1;
    while (lo <= hi) {
        i64 mid = (lo + hi) / 2; int c = strcmp(RES[mid].key, key);
        if (c == 0) return RES[mid].items;
        if (c < 0) lo = mid + 1; else hi = mid - 1;
    }
    return NULL;
}
static void push_der(Der **buf, i64 *n, i64 *cap, const Der *d) {
    if (*n >= *cap) { *cap = *cap * 2 + 64; *buf = realloc(*buf, *cap * sizeof(Der)); }
    (*buf)[(*n)++] = *d;
}
/* items: q:p/q   x:lo:hi   l:m1:m2 (anchor 1)   a:p/q:m1:m2 */
static void parse_items(const char *s, Der **buf, i64 *n, i64 *cap) {
    char tok[256];
    while (*s) {
        while (*s == ' ') s++;
        if (!*s) break;
        int L = 0; while (s[L] && s[L] != ' ' && L < 255) { tok[L] = s[L]; L++; }
        tok[L] = 0; s += L;
        Der d; int res = D_SKIP;
        if (tok[0] == 'q') {
            long long p, q; Rat r;
            if (sscanf(tok + 2, "%lld/%lld", &p, &q) == 2 && rat_set(p, q, &r)) res = finish_rat(r, &d);
        } else if (tok[0] == 'x') {
            char *e; double lo = strtod(tok + 2, &e), hi = strtod(e + 1, NULL);
            res = finish_real(lo + (hi - lo) * 0.5, lo, hi, &d);
        } else if (tok[0] == 'l') {
            char *e; double a = strtod(tok + 2, &e), b = strtod(e + 1, NULL);
            Rat one = {1, 1}; res = make_anc(one, a, b, &d);
        } else if (tok[0] == 'a') {
            long long p, q; char *e; Rat r;
            if (sscanf(tok + 2, "%lld/%lld", &p, &q) == 2 && rat_set(p, q, &r)) {
                const char *c = strchr(tok + 2, ':');
                double a = strtod(c + 1, &e), b = strtod(e + 1, NULL);
                res = make_anc(r, a, b, &d);
            }
        }
        if (res == D_OK) push_der(buf, n, cap, &d);
        else if (res != D_OUT) { fprintf(stderr, "resolution item not usable: %s (%d)\n", tok, res); exit(3); }
    }
}

static FILE *PEND;
static int FREE_TOP = 0;                 /* stage mode: free the largest tables as soon as they are built */
static long long N_UNRES;

static void build_table(int S, int keep) {
    Tab lit = Lit[S];
    i64 cap = lit.n + 1024, n = 0;
    for (int A = (S - 1) & S; A; A = (A - 1) & S) {
        int B = S ^ A;
        if (ORD) { if (!is_interval(A) || !is_interval(B) || !(A & lowbit(S))) continue; }
        else if (A > B) continue;
        cap += (i64)NFORM * (T[A].n + T[A].na) * (T[B].n + T[B].na);
    }
    Der *buf = malloc(cap * sizeof(Der));
    for (i64 i = 0; i < lit.n; i++) {
        Der d; memset(&d, 0, sizeof d);
        d.r = lit.v[i]; d.code = lit.code[i]; d.rat = lit.rq[i] != 0; d.p = lit.rp[i]; d.q = lit.rq[i];
        buf[n++] = d;
    }
    SingList *sl = &SL[S]; sl->n = 0;
    static const int ORDF[6] = {F_ADD, F_SUB, F_MUL, F_DXY, F_PXY, F_PXYN};
    for (int A = (S - 1) & S; A; A = (A - 1) & S) {
        int B = S ^ A;
        if (ORD) { if (!is_interval(A) || !is_interval(B) || !(A & lowbit(S))) continue; }
        else if (A > B) continue;
        int nf = ORD ? 6 : NFORM;
        i64 na = T[A].n + T[A].na, nb = T[B].n + T[B].na;
        for (i64 i = 0; i < na; i++) {
            Opd x = i < T[A].n ? opd_ord(A, i) : opd_anc(A, i - T[A].n);
            for (i64 j = 0; j < nb; j++) {
                Opd y = j < T[B].n ? opd_ord(B, j) : opd_anc(B, j - T[B].n);
                for (int f = 0; f < nf; f++) {
                    int form = ORD ? ORDF[f] : f;
                    Der d;
                    int res = derive(form, &x, &y, &d);
                    if (res == D_OK) { if (n >= cap) { cap = cap * 2; buf = realloc(buf, cap * sizeof(Der)); } buf[n++] = d; continue; }
                    if (res >= D_SING_Z && res <= D_SING_W) {
                        if (sl->n >= sl->cap) { sl->cap = sl->cap * 2 + 64; sl->s = realloc(sl->s, sl->cap * sizeof(Sing)); }
                        Sing s = {A, B, form, res, x, y};
                        sl->s[sl->n++] = s;
                    }
                }
            }
        }
    }
    char key[512];
    for (i64 k = 0; k < sl->n; k++) {
        sing_key(S, &sl->s[k], key);
        const char *it = find_res(key);
        if (it) parse_items(it, &buf, &n, &cap);
        else {
            #pragma omp critical
            {
                N_UNRES++;
                if (PEND) fprintf(PEND, "%s | %d %s\n", key, sl->s[k].kind, FNAME[sl->s[k].form]);
            }
        }
    }
    T[S] = merge(buf, n, keep);
    free(buf);
}

static void build_tables(int upto, int verbose) {
    for (int m = 1; m <= FULL; m++) if (mask_ok(m)) Lit[m] = make_literals(m, popc(m) <= TBL);
    double t0 = omp_get_wtime();
    for (int k = 1; k <= upto && k <= ND; k++) {
        int cnt = 0, masks[512];
        for (int m = 1; m <= FULL; m++) if (popc(m) == k && mask_ok(m)) masks[cnt++] = m;
        long long before = N_UNRES;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cnt; i++) {
            if (k == 1) { T[masks[i]] = Lit[masks[i]]; continue; }
            build_table(masks[i], k <= 4);          /* sizes <= 4 are operands of queries (|S1| <= n/2) */
            if (FREE_TOP && k == upto) {            /* stage mode: top tables are not needed afterwards */
                Tab *t = &T[masks[i]];
                free(t->v); free(t->code); free(t->ap); free(t->aq); free(t->am1); free(t->am2); free(t->axlo); free(t->axhi);
                if (t->rp) { free(t->rp); free(t->rq); free(t->lo); free(t->hi); free(t->lnlo); free(t->lnhi); }
                t->v = NULL; t->code = NULL; t->ap = t->aq = NULL; t->am1 = t->am2 = t->axlo = t->axhi = NULL;
                t->rp = t->rq = NULL; t->lo = t->hi = t->lnlo = t->lnhi = NULL;
            }
        }
        i64 tot = 0, tota = 0, ns = 0;
        for (int i = 0; i < cnt; i++) { tot += T[masks[i]].n; tota += T[masks[i]].na; ns += SL[masks[i]].n; }
        if (verbose)
            fprintf(stderr, "tables size %d: %d sets, %lld classes + %lld anchored, %lld singular derivations "
                    "(%lld unresolved), %.1fs\n", k, cnt, tot, tota, ns, N_UNRES - before, omp_get_wtime() - t0);
        if (N_UNRES > before && k < upto) { fprintf(stderr, "unresolved singular derivations at size %d\n", k); exit(4); }
    }
}

/* ================================================================== window queries (Lemma 8, Theorem B) */
enum { Q_ADD, Q_SUB, Q_MUL, Q_SL, Q_LS, Q_PS, Q_PL, NQF };
static const char *QNAME[NQF] = {"add", "sub", "mul", "sl", "ls", "ps", "pl"};
typedef struct { int mask, fam, side, kind; double k1, k2, k3; i64 ap, aq; } Step;
typedef struct { Step st[MAXDEPTH]; int depth; } Path;
/* a window for l: always an x-enclosure [lo, hi]; optionally an enclosure of ln l [Llo, Lhi]; optionally
   an anchored form l in A e^[n1, n2] with A an exact rational (precise even below double resolution) */
typedef struct { double lo, hi, Llo, Lhi; int hasA; Rat A; double n1, n2; int lx; } Win;
/* lx (efficiency only): bit 0 means lo = exp_lo(Llo) unclamped, so ln_lo(lo) <= Llo and the log window's
   lower end is exactly Llo; bit 1 means hi = exp_hi(Lhi) unclamped, so ln_hi(hi) >= Lhi and the upper end is Lhi.
   win_logs then returns the same values without computing the logarithms. */

typedef struct { char **l; i64 n, cap; } Lines;
static Lines OUT[512];
#define NWPAD 16                    /* one 128-byte line per thread: no false sharing between threads */
static long long NWIN[512 * NWPAD];

#define SPILL_AT 1000000
static FILE *CUR_OUT; static const char *CUR_PREFIX = ""; static i64 N_SPILLED = 0;
/* MOTORA_DEDUP=1 (efficiency only): a line already emitted in the same work item (same first step, same
   thread) is not emitted again.  Repeated lines come from different windows of one step (|s - l| gives two,
   l^(+-s) two) and so from the same item; the set of distinct lines written is unchanged. */
static int DEDUP = 0;
typedef struct { uint64_t *h; char **s; i64 cap, n; char *arena; i64 acap, an; i64 dropped; char pad[64]; } DSet;
static DSet DS[512];
static uint64_t fnv(const char *s) { uint64_t h = 1469598103934665603ULL; for (; *s; s++) h = (h ^ (unsigned char)*s) * 1099511628211ULL; return h; }
static void dset_reset(int tid) { DSet *d = &DS[tid]; if (d->n) { memset(d->h, 0, d->cap * sizeof(uint64_t)); d->n = 0; } d->an = 0; }
static int dset_seen(int tid, const char *line) {      /* 1 if line was already added since the last reset */
    DSet *d = &DS[tid];
    if (2 * (d->n + 1) > d->cap) {                      /* grow (rehash) */
        i64 oc = d->cap, nc = oc ? 2 * oc : 1024; uint64_t *oh = d->h; char **os = d->s;
        d->h = calloc(nc, sizeof(uint64_t)); d->s = calloc(nc, sizeof(char *)); d->cap = nc;
        for (i64 i = 0; i < oc; i++) if (oh[i]) { i64 j = oh[i] & (nc - 1); while (d->h[j]) j = (j + 1) & (nc - 1); d->h[j] = oh[i]; d->s[j] = os[i]; }
        free(oh); free(os);
    }
    uint64_t h = fnv(line) | 1;
    i64 j = h & (d->cap - 1);
    while (d->h[j]) { if (d->h[j] == h && !strcmp(d->s[j], line)) return 1; j = (j + 1) & (d->cap - 1); }
    i64 len = strlen(line) + 1;
    if (d->an + len > d->acap) {                        /* arena full: start a new one (old strings stay valid) */
        d->acap = d->acap * 2 + (1 << 20) + len; d->arena = malloc(d->acap); d->an = 0;   /* (old arenas leak per item: bounded) */
    }
    char *c = d->arena + d->an; memcpy(c, line, len); d->an += len;
    d->h[j] = h; d->s[j] = c; d->n++;
    return 0;
}
static void emit(const Path *p, const char *leaf);
static void emit_(const Path *p, const char *leaf);
static void emit(const Path *p, const char *leaf) {
    if (!DEDUP) { emit_(p, leaf); return; }
    char buf[4096]; int k = 0;
    for (int i = 0; i < p->depth; i++) {
        const Step *s = &p->st[i];
        if (s->kind == 0) k += sprintf(buf + k, "%d:%s:%d:o:%a:%d ", s->mask, QNAME[s->fam], s->side, s->k1, (int)s->k2);
        else k += sprintf(buf + k, "%d:%s:%d:a:%lld/%lld:%a:%a ", s->mask, QNAME[s->fam], s->side,
                          (long long)s->ap, (long long)s->aq, s->k1, s->k2);
    }
    sprintf(buf + k, "%s", leaf);
    int tid = omp_get_thread_num();
    if (dset_seen(tid, buf)) { DS[tid].dropped++; return; }
    emit_(p, leaf);
}
static void emit_(const Path *p, const char *leaf) {
    char buf[4096]; int k = 0;
    for (int i = 0; i < p->depth; i++) {
        const Step *s = &p->st[i];
        if (s->kind == 0) k += sprintf(buf + k, "%d:%s:%d:o:%a:%d ", s->mask, QNAME[s->fam], s->side, s->k1, (int)s->k2);
        else k += sprintf(buf + k, "%d:%s:%d:a:%lld/%lld:%a:%a ", s->mask, QNAME[s->fam], s->side,
                          (long long)s->ap, (long long)s->aq, s->k1, s->k2);
    }
    sprintf(buf + k, "%s", leaf);
    Lines *L = &OUT[omp_get_thread_num()];
    if (L->n >= L->cap) { L->cap = L->cap * 2 + 64; L->l = realloc(L->l, L->cap * sizeof(char *)); }
    L->l[L->n++] = strdup(buf);
    /* very large windows: once a thread holds many lines, spill them (unsorted) to the output file,
       so that memory stays bounded; small outputs stay sorted and deterministic */
    if (L->n >= SPILL_AT) {
        #pragma omp critical(spill)
        {
            for (i64 i = 0; i < L->n; i++) { fprintf(CUR_OUT, "%sC %s\n", CUR_PREFIX, L->l[i]); free(L->l[i]); }
            N_SPILLED += L->n; L->n = 0;
        }
    }
}

static inline int rat_less_eq_key(i64 p1, i64 q1, i64 p2, i64 q2) { Rat a = {p1, q1}, b = {p2, q2}; return rat_cmp(a, b); }

/* the log window of a leaf query, computed only when a class needs it (same values as win_logs) */
typedef struct { const Win *w; double La, Lb; int have; } LZ;
static void win_logs(const Win *W, double *La, double *Lb);
static inline void lz_get(LZ *z) { if (!z->have) { win_logs(z->w, &z->La, &z->Lb); z->have = 1; } }
static void leaf_lookup(int S, const Win *w, LZ *z, const Path *p, int is_lit) {
    Tab *t = is_lit ? &Lit[S] : &T[S];
    char leaf[256];
    double a = w->lo, b = w->hi;
    if (t->n) {
        double delta = ldexp(1.0, t->maxcode - 52);
        double klo = pred(a * (1 - 2 * delta)); if (klo < 0) klo = 0;
        double khi = succ(b * (1 + 2 * delta));
        i64 lo = 0, hi = t->n;
        while (lo < hi) { i64 mid = (lo + hi) >> 1; if (t->v[mid] < klo) lo = mid + 1; else hi = mid; }
        for (i64 i = lo; i < t->n && t->v[i] <= khi; i++) {
            double h1, h2; hull_of(t->v[i], t->code[i], &h1, &h2);
            if (h2 < a || h1 > b) continue;
            if (t->code[i] == 0) {                       /* an exact double: decide with the precise forms */
                double v = t->v[i];
                if (w->hasA) {
                    double al, ah; rat_hull(w->A, &al, &ah);
                    if (al == ah && al == v) { if (!(w->n1 <= 0 && 0 <= w->n2)) continue; }
                } else if (v == 1) { lz_get(z); if (!(z->La <= 0 && 0 <= z->Lb)) continue; }
            } else if (h1 > 0.5 && h2 < 2) {             /* near 1: the log window decides */
                lz_get(z); if (ln_hi(h2) < z->La || ln_lo(h1) > z->Lb) continue;
            }
            sprintf(leaf, "%d:%s:o:%a:%d", S, is_lit ? "lit" : "leaf", t->v[i], t->code[i]);
            emit(p, leaf);
        }
    }
    if (t->na) {
        /* anchored classes: x-enclosure within [A e^-MU, A e^MU]; scan anchors in the window's range */
        double klo = pred(a * 0.998), khi = succ(b * 1.002);
        i64 lo = 0, hi = t->na;
        while (lo < hi) { i64 mid = (lo + hi) >> 1; if (rat_double((Rat){t->ap[mid], t->aq[mid]}) < klo) lo = mid + 1; else hi = mid; }
        for (i64 i = lo; i < t->na; i++) {
            Rat q = {t->ap[i], t->aq[i]};
            if (rat_double(q) > khi) break;
            if (t->axhi[i] < a || t->axlo[i] > b) continue;
            if (w->hasA && q.p == w->A.p && q.q == w->A.q) {
                if (t->am2[i] < w->n1 || t->am1[i] > w->n2) continue;
            } else {
                double l1, l2; anc_ln(q, t->am1[i], t->am2[i], &l1, &l2);
                lz_get(z); if (l2 < z->La || l1 > z->Lb) continue;
            }
            sprintf(leaf, "%d:leaf:a:%lld/%lld:%a:%a", S, (long long)q.p, (long long)q.q, t->am1[i], t->am2[i]);
            emit(p, leaf);
        }
    }
}

static int families(int side, int *f) {
    int n = 0;
    f[n++] = Q_ADD; f[n++] = Q_SUB; f[n++] = Q_MUL;
    if (side != 0) { f[n++] = Q_SL; f[n++] = Q_PS; }
    if (side != 1) { f[n++] = Q_LS; f[n++] = Q_PL; }
    return n;
}

/* l = 0 is admissible only where the family is defined at l = 0 and can give a value in [a, b];
   every other l is >= 1e-300 (in range), so a window may start there */
static inline Win mkwin(double l, double h, int allow0) {
    Win w; memset(&w, 0, sizeof w);
    if (l != l) l = 0;                               /* NaN bounds (never expected): widen, never drop */
    if (h != h) h = INFINITY;
    if (l < 0) l = 0;
    if (!allow0 && l < RLO_SAFE) l = RLO_SAFE;
    if (h > RHI_SAFE) h = RHI_SAFE;
    w.lo = l; w.hi = h; w.Llo = NAN; w.Lhi = NAN;
    return w;
}
static inline void add_win(Win *ws, int *n, Win w) { if (w.lo <= w.hi) ws[(*n)++] = w; }
static inline void push(Win *ws, int *n, double l, double h, int allow0) { add_win(ws, n, mkwin(l, h, allow0)); }
static inline void push_log(Win *ws, int *n, double u1, double u2, int allow0) {
    double e1 = exp_lo(u1), e2 = u2 == INFINITY ? INFINITY : exp_hi(u2);
    Win w = mkwin(e1, e2, allow0);
    w.Llo = u1; w.Lhi = u2;
    w.lx = (u1 == u1 && w.lo == e1 ? 1 : 0) | (u2 == u2 && u2 != INFINITY && w.hi == e2 ? 2 : 0);
    if (u1 > -MU_MAX && u2 < MU_MAX) { w.hasA = 1; w.A.p = 1; w.A.q = 1; w.n1 = u1; w.n2 = u2; }
    add_win(ws, n, w);
}
/* l in A e^[n1, n2] */
static inline void push_anc(Win *ws, int *n, Rat A, double n1, double n2) {
    if (A.p <= 0) return;
    double x1, x2; anc_x(A, n1, n2, &x1, &x2);
    Win w = mkwin(x1, x2, 0);
    double l1, l2; anc_ln(A, n1, n2, &l1, &l2);
    w.Llo = l1; w.Lhi = l2;
    if (n2 < MU_MAX && n1 > -MU_MAX) { w.hasA = 1; w.A = A; w.n1 = n1; w.n2 = n2; }
    add_win(ws, n, w);
}
/* l = Q + delta, delta in [d1, d2], l >= 0 */
static inline void push_affine(Win *ws, int *n, Rat Q, double d1, double d2) {
    if (Q.p == 0) { push(ws, n, d1, d2, 1); return; }
    double ql, qh; rat_hull(Q, &ql, &qh);
    if (Q.p > 0 && d2 <= ql * TINY_REL && d1 >= -ql * TINY_REL) {
        double u1, u2, n1, n2;
        idiv(d1, d2, ql, qh, &u1, &u2);
        log1p_hull(u1, u2, &n1, &n2);
        push_anc(ws, n, Q, n1, n2);
        return;
    }
    push(ws, n, pred(ql + d1), succ(qh + d2), 1);
}
/* the target window as Q + delta (anchored window: A + A(e^n - 1)) */
static inline int win_decomp(const Win *w, Rat *Q, double *d1, double *d2) {
    if (!w->hasA) return 0;
    double e1, e2, ql, qh;
    expm1_hull(w->n1, w->n2, &e1, &e2); rat_hull(w->A, &ql, &qh);
    *Q = w->A; imul(ql, qh, e1, e2, d1, d2);
    return 1;
}

/* preimage windows for l (Lemma 8); *cst >= 0 when the family is constant and its value is in W */
static int windows(int fam, const Opd *s, const Win *W, double La, double Lb, Win *w, double *cst) {
    int n = 0; *cst = -1;
    double a = W->lo, b = W->hi, slo = s->lo, shi = s->hi;
    Rat QW = {0, 1}; double wd1 = 0, wd2 = 0;
    int wa = win_decomp(W, &QW, &wd1, &wd2);
    int sa = s->kind != K_ORD;                       /* s rational or anchored */
    Rat Qs = {0, 1}; double sd1 = 0, sd2 = 0;
    if (sa) decomp(s, &Qs, &sd1, &sd2);
    switch (fam) {
    case Q_ADD:                                      /* l = x - s */
        if (!wa) { QW.p = 0; QW.q = 1; wd1 = a; wd2 = b; }   /* x = 0 + [a, b] */
        if (sa) { Rat Q; double d1, d2; if (rat_sub(QW, Qs, &Q)) { isub(wd1, wd2, sd1, sd2, &d1, &d2); push_affine(w, &n, Q, d1, d2); break; } }
        push(w, &n, pred(a - shi), succ(b - slo), 1); break;
    case Q_SUB: {                                    /* l = s + x  or  l = s - x */
        Rat Q1, Q2;
        if (!wa) { QW.p = 0; QW.q = 1; wd1 = a; wd2 = b; }
        if (sa && rat_add(Qs, QW, &Q1) && rat_sub(Qs, QW, &Q2)) {
            double d1, d2;
            iadd(sd1, sd2, wd1, wd2, &d1, &d2); push_affine(w, &n, Q1, d1, d2);
            isub(sd1, sd2, wd1, wd2, &d1, &d2); push_affine(w, &n, Q2, d1, d2);
            break;
        }
        push(w, &n, pred(slo + a), succ(shi + b), 1);
        push(w, &n, pred(slo - b), succ(shi - a), 1); break; }
    case Q_MUL:                                      /* l = x / s */
        if (s->zero) { if (a == 0) *cst = 0; break; }
        if (wa && sa) { Rat Q; if (rat_div(QW, s->r, &Q)) { double m1 = W->n1, m2 = W->n2; if (s->kind == K_ANC) isub(W->n1, W->n2, s->m1, s->m2, &m1, &m2); push_anc(w, &n, Q, m1, m2); break; } }
        if (La > -INFINITY) push_log(w, &n, pred(La - s->lnhi), succ(Lb - s->lnlo), 0);
        else push(w, &n, pred(a / shi), succ(b / slo), a == 0);
        break;
    case Q_SL:                                       /* l = s / x */
        if (s->zero) { if (a == 0) *cst = 0; break; }
        if (b == 0) break;
        if (wa && sa) { Rat Q; if (rat_div(s->r, QW, &Q)) { double m1, m2; if (s->kind == K_ANC) isub(s->m1, s->m2, W->n1, W->n2, &m1, &m2); else { m1 = -W->n2; m2 = -W->n1; } push_anc(w, &n, Q, m1, m2); break; } }
        if (La > -INFINITY) push_log(w, &n, pred(s->lnlo - Lb), succ(s->lnhi - La), 0);
        else push(w, &n, pred(slo / b), a > 0 ? succ(shi / a) : INFINITY, 0);
        break;
    case Q_LS:                                       /* l = x s */
        if (s->zero) break;
        if (wa && sa) { Rat Q; if (rat_mul(QW, s->r, &Q)) { double m1 = W->n1, m2 = W->n2; if (s->kind == K_ANC) iadd(W->n1, W->n2, s->m1, s->m2, &m1, &m2); push_anc(w, &n, Q, m1, m2); break; } }
        if (La > -INFINITY) push_log(w, &n, pred(La + s->lnlo), succ(Lb + s->lnhi), 0);
        else push(w, &n, pred(a * slo), succ(b * shi), a == 0);
        break;
    case Q_PS: {                                     /* l = |t|, s^t = x */
        if (s->zero) { if (a == 0) *cst = 0; break; }
        if (s->one) { if (a <= 1 && 1 <= b) *cst = 1; break; }
        if (b == 0) break;
        double l1 = s->lnlo, l2 = s->lnhi, tlo, thi;
        if (La == -INFINITY) {
            double q1 = Lb / l1, q2 = Lb / l2;
            if (l1 > 0) { tlo = -INFINITY; thi = succ(q1 > q2 ? q1 : q2); }
            else { thi = INFINITY; tlo = pred(q1 < q2 ? q1 : q2); }
        } else idiv(La, Lb, l1, l2, &tlo, &thi);
        int allow0 = a <= 1 && 1 <= b;
        if (tlo >= 0) push(w, &n, tlo, thi, allow0);
        else if (thi <= 0) push(w, &n, -thi, -tlo, allow0);
        else push(w, &n, 0, -tlo > thi ? -tlo : thi, allow0);
        break; }
    case Q_PL: {                                     /* l^(+-s) = x */
        if (s->zero) { if (a <= 1 && 1 <= b) *cst = 1; break; }
        if (b == 0) { push(w, &n, 0, 0, 1); break; }     /* l^s = 0 only for l = 0 (and +s) */
        for (int sg = 1; sg >= -1; sg -= 2) {
            if (wa && s->kind == K_RAT) {            /* l = A^(1/e) e^(n/e), exact when A^(1/e) is rational */
                Rat inv = {sg * s->r.q, s->r.p}, Q;
                if (inv.q < 0) { inv.p = -inv.p; inv.q = -inv.q; }
                if (rat_pow(QW, inv, &Q)) {
                    double el, eh, m1, m2; rat_hull(inv, &el, &eh);
                    imul(W->n1, W->n2, el, eh, &m1, &m2);
                    push_anc(w, &n, Q, m1, m2);
                    continue;
                }
            }
            double u1, u2;
            double A1 = sg > 0 ? La : -Lb, A2 = sg > 0 ? Lb : -La;     /* ln l = (+-ln x) / s */
            if (A1 == -INFINITY || A1 == INFINITY || A2 == -INFINITY || A2 == INFINITY) {
                if (sg > 0) { double q1 = Lb / slo, q2 = Lb / shi; u1 = -INFINITY; u2 = succ(q1 > q2 ? q1 : q2); }
                else { double q1 = -Lb / slo, q2 = -Lb / shi; u1 = pred(q1 < q2 ? q1 : q2); u2 = INFINITY; }
            } else idiv(A1, A2, slo, shi, &u1, &u2);
            push_log(w, &n, u1, u2, sg > 0 && a == 0);
        }
        break; }
    }
    return n;
}

static void query(int S, const Win *W, Path *p);

static void query_split(int S1, int S2, int side, const Win *W, double La, double Lb, Path *p, i64 only_i) {
    int fam[8], nf = families(side, fam);
    Tab *t = &T[S1];
    i64 tot = t->n + t->na;
    for (i64 i = (only_i >= 0 ? only_i : 0); i < (only_i >= 0 ? only_i + 1 : tot); i++) {
        Opd s = i < t->n ? opd_ord(S1, i) : opd_anc(S1, i - t->n);
        Step st; memset(&st, 0, sizeof st); st.mask = S1; st.side = side;
        if (s.kind == K_ANC) { st.kind = 1; st.k1 = s.m1; st.k2 = s.m2; st.ap = s.r.p; st.aq = s.r.q; }
        else { st.kind = 0; st.k1 = s.v; st.k2 = s.code; }
        for (int k = 0; k < nf; k++) {
            Win w[4]; double cst;
            int nw = windows(fam[k], &s, W, La, Lb, w, &cst);
            st.fam = fam[k];
            p->st[p->depth++] = st;
            if (cst >= 0) { char leaf[64]; sprintf(leaf, "%d:const:%g", S2, cst); emit(p, leaf); }
            for (int j = 0; j < nw; j++) { NWIN[omp_get_thread_num() * NWPAD]++; query(S2, &w[j], p); }
            p->depth--;
        }
    }
}

static int query_splits(int S, int *s1, int *s2, int *side) {
    int n = 0, pm = popc(S);
    if (ORD) {
        for (int A = (S - 1) & S; A; A = (A - 1) & S) {
            int B = S ^ A;
            if (!is_interval(A) || !is_interval(B) || !(A & lowbit(S))) continue;
            if (popc(A) <= popc(B)) { s1[n] = A; s2[n] = B; side[n] = 1; }
            else { s1[n] = B; s2[n] = A; side[n] = 0; }
            n++;
        }
    } else {
        for (int A = (S - 1) & S; A; A = (A - 1) & S) {
            int B = S ^ A, pa = popc(A);
            if (2 * pa > pm || (2 * pa == pm && A > B)) continue;
            s1[n] = A; s2[n] = B; side[n] = -1; n++;
        }
    }
    return n;
}

static void win_logs(const Win *W, double *La, double *Lb) {
    *La = (W->lx & 1) ? W->Llo : W->lo > 0 ? ln_lo(W->lo) : -INFINITY;
    *Lb = (W->lx & 2) ? W->Lhi : W->hi > 0 ? ln_hi(W->hi) : -INFINITY;
    if (W->Llo == W->Llo && W->Llo > *La) *La = W->Llo;
    if (W->Lhi == W->Lhi && W->Lhi < *Lb) *Lb = W->Lhi;
}

static void query(int S, const Win *W, Path *p) {
    if (p->depth >= MAXDEPTH - 1) { fprintf(stderr, "depth overflow\n"); exit(5); }
    if (popc(S) <= TBL && !(W->Llo == W->Llo) && !(W->Lhi == W->Lhi)) {
        /* no log interval: La = ln_lo(lo) <= ln_hi(hi) = Lb (or La = -inf), so the emptiness test below cannot
           fire, and the logarithms are computed only if a class needs them */
        LZ z = {W, 0, 0, 0}; leaf_lookup(S, W, &z, p, 0); return;
    }
    double La, Lb; win_logs(W, &La, &Lb);
    if (La > Lb) return;
    LZ z = {W, La, Lb, 1};
    if (popc(S) <= TBL) { leaf_lookup(S, W, &z, p, 0); return; }
    leaf_lookup(S, W, &z, p, 1);
    int s1[256], s2[256], sd[256];
    int ns = query_splits(S, s1, s2, sd);
    for (int k = 0; k < ns; k++) query_split(s1[k], s2[k], sd[k], W, La, Lb, p, -1);
}

static int cmp_str(const void *a, const void *b) { return strcmp(*(char *const *)a, *(char *const *)b); }

static void query_top(double a, double b, FILE *out, const char *prefix) {
    double t0 = omp_get_wtime();
    CUR_OUT = out; CUR_PREFIX = prefix; N_SPILLED = 0;
    Win W = mkwin(a, b, 1);
    double La, Lb; win_logs(&W, &La, &Lb);
    if (popc(FULL) <= TBL) {
        Path p; p.depth = 0;
        LZ z = {&W, La, Lb, 1}; leaf_lookup(FULL, &W, &z, &p, 0);
    } else {
        Path p0; p0.depth = 0;
        dset_reset(0);
        LZ z = {&W, La, Lb, 1}; leaf_lookup(FULL, &W, &z, &p0, 1);
        int s1[256], s2[256], sd[256];
        int ns = query_splits(FULL, s1, s2, sd);
        i64 nitems = 0;
        for (int k = 0; k < ns; k++) nitems += T[s1[k]].n + T[s1[k]].na;
        i64 *ik = malloc(nitems * sizeof(i64)), *ii = malloc(nitems * sizeof(i64)), q = 0;
        for (int k = 0; k < ns; k++) for (i64 i = 0; i < T[s1[k]].n + T[s1[k]].na; i++) { ik[q] = k; ii[q] = i; q++; }
        /* deterministic shuffle of the work items, so that the progress lines are representative
           (the output is sorted at the end, so the order does not change the result) */
        i64 *perm = malloc(nitems * sizeof(i64));
        for (i64 j = 0; j < nitems; j++) perm[j] = j;
        unsigned long long st = 88172645463325252ULL;
        for (i64 j = nitems - 1; j > 0; j--) {
            st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            i64 r = (i64)(st % (unsigned long long)(j + 1)), tmp = perm[j]; perm[j] = perm[r]; perm[r] = tmp;
        }
        i64 nrun = nitems;                     /* MOTORA_MAXITEMS (tests only): the first K items of the shuffled order */
        if (getenv("MOTORA_MAXITEMS") && atoll(getenv("MOTORA_MAXITEMS")) < nitems) nrun = atoll(getenv("MOTORA_MAXITEMS"));
        i64 done = 0, step = nrun / 100 > 0 ? nrun / 100 : 1;
        #pragma omp parallel for schedule(dynamic, 16)
        for (i64 j = 0; j < nrun; j++) {
            i64 it = perm[j];
            Path p; p.depth = 0;
            if (DEDUP) dset_reset(omp_get_thread_num());
            int k = (int)ik[it];
            query_split(s1[k], s2[k], sd[k], &W, La, Lb, &p, ii[it]);
            i64 d;
            #pragma omp atomic capture
            d = ++done;
            if (d % step == 0) {
                double el = omp_get_wtime() - t0;
                #pragma omp critical
                { fprintf(stderr, "%sprogress %lld/%lld (%.0f%%)  %.2f h elapsed, about %.2f h left\n", prefix, (long long)d,
                          (long long)nrun, 100.0 * d / nrun, el / 3600, el / 3600 * ((double)nrun / d - 1)); fflush(stderr); }
            }
        }
        free(perm);
        free(ik); free(ii);
    }
    i64 tot = 0; long long nw = 0;
    for (int t = 0; t < 512; t++) { tot += OUT[t].n; nw += NWIN[t * NWPAD]; }
    i64 spilled = N_SPILLED; N_SPILLED = 0;
    char **all = malloc((tot + 1) * sizeof(char *)); i64 k = 0;
    for (int t = 0; t < 512; t++) for (i64 i = 0; i < OUT[t].n; i++) all[k++] = OUT[t].l[i];
    qsort(all, tot, sizeof(char *), cmp_str);
    i64 nprint = 0, ndrop = 0;
    for (int t = 0; t < 512; t++) { ndrop += DS[t].dropped; DS[t].dropped = 0; }
    for (i64 i = 0; i < tot; i++) {
        if (DEDUP && i > 0 && !strcmp(all[i], all[i - 1])) { ndrop++; free(all[i - 1]); all[i - 1] = NULL; continue; }   /* sorted: equal lines adjacent */
        if (i > 0 && all[i - 1]) free(all[i - 1]);
        fprintf(out, "%sC %s\n", prefix, all[i]); nprint++;
    }
    if (tot) free(all[tot - 1]);
    free(all);
    for (int t = 0; t < 512; t++) { OUT[t].n = 0; NWIN[t * NWPAD] = 0; }
    fprintf(stderr, "%scandidates: %lld, windows queried: %lld, %.1fs\n", prefix, (long long)(DEDUP ? nprint + spilled : tot + spilled), nw, omp_get_wtime() - t0);
    if (DEDUP) fprintf(stderr, "%srepeated lines not written: %lld\n", prefix, (long long)ndrop);
    fflush(out); fflush(stderr);
}

#ifndef NO_MAIN
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: motora stage n ord k res pend | multi n ord res windows out\n"); return 1; }
    ND = atoi(argv[2]); ORD = atoi(argv[3]); FULL = (1 << ND) - 1;
    setvbuf(stderr, NULL, _IONBF, 0);
    DEDUP = getenv("MOTORA_DEDUP") && atoi(getenv("MOTORA_DEDUP"));
#if defined(__FMA__) && defined(__GNUC__)
    if (!__builtin_cpu_supports("fma")) { fprintf(stderr, "this binary needs a CPU with FMA (built with -mfma)\n"); return 6; }
#endif
    RLO_SAFE = pred(1e-300); RHI_SAFE = succ(1e300);
    if (!strcmp(argv[1], "stage")) {
        int k = atoi(argv[4]);
        load_res(argv[5]);
        PEND = fopen(argv[6], "w");
        FREE_TOP = 1;
        build_tables(k, 1);
        fclose(PEND);
        fprintf(stderr, "pending singular derivations: %lld\n", N_UNRES);
        return 0;
    }
    if (!strcmp(argv[1], "multi")) {
        load_res(argv[4]);
        build_tables(TBL, 1);
        if (N_UNRES) { fprintf(stderr, "ERROR: %lld unresolved singular derivations\n", N_UNRES); return 4; }
        FILE *in = fopen(argv[5], "r"), *out = fopen(argv[6], "w");
        char sa[64], sb[64]; int w = 0;
        while (fscanf(in, "%63s %63s", sa, sb) == 2) {
            char pre[32]; sprintf(pre, "W%d ", w++);
            query_top(strtod(sa, NULL), strtod(sb, NULL), out, pre);
        }
        fclose(in); fclose(out);
        return 0;
    }
    return 1;
}
#endif
