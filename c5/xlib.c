/*
 * xlib.c - C5 helper library (motor X): reuses motor A unchanged (motora.c with NO_MAIN).
 *
 *  - x_init: builds motor A's in-range tables T(S) for |S| <= upto (with its resolution file);
 *  - x_count / x_dump: every class of T(S) as an enclosure [llo, lhi] of ln(value) (outward, from
 *    motor A's enclosures and its rigorous ln bounds), with flags: 1 exact zero, 2 exact one;
 *  - v_exp_lo/hi, v_ln_lo/hi: rigorous vectorised bounds of e^x and ln x (motor A's exp_lo, exp_hi,
 *    ln_lo, ln_hi, built on CORE-MATH's correctly rounded pow; Lemma 4 of C4);
 *  - v_pred / v_succ: neighbouring doubles (outward rounding of + - * / results).
 */
#define NO_MAIN
#include "motora.c"     /* compile with -I <certifier src> */

/* motor A's build_tables, keeping the operand arrays (exact rationals, enclosures) at every size */
static void build_tables_keep(int upto) {
    for (int m = 1; m <= FULL; m++) if (mask_ok(m)) Lit[m] = make_literals(m, 1);
    for (int k = 1; k <= upto && k <= ND; k++) {
        int cnt = 0, masks[512];
        for (int m = 1; m <= FULL; m++) if (popc(m) == k && mask_ok(m)) masks[cnt++] = m;
        long long before = N_UNRES;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cnt; i++) {
            if (k == 1) { T[masks[i]] = Lit[masks[i]]; continue; }
            build_table(masks[i], 1);
        }
        i64 tot = 0;
        for (int i = 0; i < cnt; i++) tot += T[masks[i]].n + T[masks[i]].na;
        fprintf(stderr, "tables size %d: %d sets, %lld classes (kept), %lld unresolved\n", k, cnt, (long long)tot, N_UNRES - before);
        if (N_UNRES > before) { fprintf(stderr, "unresolved singular derivations at size %d\n", k); exit(4); }
    }
}

void x_init(int nd, int ord, const char *res, int upto) {
    ND = nd; ORD = ord; FULL = (1 << ND) - 1;
    RLO_SAFE = pred(1e-300); RHI_SAFE = succ(1e300);
    load_res(res);
    build_tables_keep(upto);
    if (N_UNRES) { fprintf(stderr, "ERROR: %lld unresolved singular derivations\n", N_UNRES); exit(4); }
}

long x_count(int mask) { return (long)(T[mask].n + T[mask].na); }

static void lnhull(double lo, double hi, double *a, double *b) {
    *a = lo > 0 ? ln_lo(lo) : -INFINITY;
    *b = ln_hi(hi);
    /* near 1: ln x = log1p(x - 1) with x - 1 exact (Sterbenz); log1p(u) in [u - u^2, u] for |u| < 2^-9 */
    if (lo > 1 - 0x1p-9 && lo < 1 + 0x1p-9) { double u = lo - 1; double l = pred(u - succ(u * u)); if (l > *a) *a = l; }
    if (hi > 1 - 0x1p-9 && hi < 1 + 0x1p-9) { double u = hi - 1; if (u < *b) *b = u; }
}

void x_dump(int mask, double *llo, double *lhi, unsigned char *flags) {
    Tab *t = &T[mask];
    for (i64 i = 0; i < t->n; i++) {
        unsigned char f = 0;
        double lo, hi;
        if (t->rq && t->rq[i]) {                                  /* exact rational p/q */
            Rat q = {t->rp[i], t->rq[i]};
            if (q.p == 0) { llo[i] = lhi[i] = -INFINITY; flags[i] = 1; continue; }
            if (q.p == q.q) { llo[i] = lhi[i] = 0; flags[i] = 2; continue; }
            rat_hull(q, &lo, &hi);
        } else {
            if (t->code[i] == 0 && t->v[i] == 0) { llo[i] = lhi[i] = -INFINITY; flags[i] = 1; continue; }   /* exact 0 */
            if (t->code[i] == 0 && t->v[i] == 1) { llo[i] = lhi[i] = 0; flags[i] = 2; continue; }           /* exact 1 */
            hull_of(t->v[i], t->code[i], &lo, &hi);
        }
        lnhull(lo, hi, &llo[i], &lhi[i]);
        if (t->lnlo && t->lnhi && t->lnlo[i] <= t->lnhi[i] && isfinite(t->lnlo[i]) && isfinite(t->lnhi[i])) {
            /* motor A's own ln enclosure (used near 1, where it is much tighter) */
            if (t->lnlo[i] > llo[i]) llo[i] = t->lnlo[i];
            if (t->lnhi[i] < lhi[i]) lhi[i] = t->lnhi[i];
        }
        flags[i] = f;
    }
    for (i64 j = 0; j < t->na; j++) {
        i64 i = t->n + j;
        Rat q = {t->ap[j], t->aq[j]};
        double l1, l2; anc_ln(q, t->am1[j], t->am2[j], &l1, &l2);
        llo[i] = l1; lhi[i] = l2; flags[i] = 0;
    }
}

/* exact rationals of the classes (q = 0: not an exact rational); tables of size <= 4 only */
void x_rats(int mask, long long *p, long long *q) {
    Tab *t = &T[mask];
    for (i64 i = 0; i < t->n; i++) { p[i] = t->rq ? t->rp[i] : 0; q[i] = t->rq ? t->rq[i] : 0; }
    for (i64 j = 0; j < t->na; j++) { p[t->n + j] = 0; q[t->n + j] = 0; }
}

/* elementwise pure functions: the OpenMP split of the loop (large arrays only) does not change any result */
#define VPAR _Pragma("omp parallel for schedule(static) if (n > 8192)")
void v_exp_lo(const double *x, double *y, long n) { VPAR for (long i = 0; i < n; i++) y[i] = exp_lo(x[i]); }
void v_exp_hi(const double *x, double *y, long n) { VPAR for (long i = 0; i < n; i++) y[i] = exp_hi(x[i]); }
void v_ln_lo(const double *x, double *y, long n) {
    VPAR for (long i = 0; i < n; i++) y[i] = x[i] > 0 ? ln_lo(x[i]) : (x[i] == 0 ? -INFINITY : NAN);
}
void v_ln_hi(const double *x, double *y, long n) {
    VPAR for (long i = 0; i < n; i++) y[i] = x[i] > 0 ? ln_hi(x[i]) : (x[i] == 0 ? -INFINITY : NAN);
}
void v_pred(const double *x, double *y, long n) { for (long i = 0; i < n; i++) y[i] = pred(x[i]); }
void v_succ(const double *x, double *y, long n) { for (long i = 0; i < n; i++) y[i] = succ(x[i]); }
