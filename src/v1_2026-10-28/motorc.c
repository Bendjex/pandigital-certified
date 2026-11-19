/*
 * Motor C: rigorous numerical pre-filter for motor A's candidate paths (C3, Lemma 15).
 *
 * For every candidate path it computes, with Arb balls (and exact rationals where the value is rational),
 * a SUPERSET of the values |val E| of the expressions E covered by the path that may lie in the window W0,
 * following exactly the structure of motor B's verification (verify.py): class values of S1 by complete
 * enumeration, preimage windows for l, recursion, leaf queries.  A path is DROPPED only when that superset is
 * empty, i.e. when every value it could produce is certainly outside W0, certainly out of range, or certainly
 * undefined.  Every other path (including any path where a ball could not decide something: "unknown") is
 * KEPT and goes to motor B, which decides it exactly.  So motor C never changes the result; it only removes
 * paths that motor B would also find empty.
 *
 * usage: motorc n ord window_file candidates_file kept_file [precision_bits]
 *        environment: MOTORC_JOBS (worker processes, tables shared copy-on-write), MOTORC_TPREC (table precision)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>

static slong PREC = 256;           /* working precision in bits (argument 6) */
#define MAXD 9
#define QBITS 4000                 /* rationals larger than this are kept as balls only */

static int ND, ORD;
static arb_t LO_R, HI_R, BIG, ZERO_A;  /* range bounds 1e-300, 1e300 */
static long long N_UNKNOWN = 0;
static int UPG = 0;                       /* table values are upgraded to the working precision when used */
#include <time.h>
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double T_FULL = 0, T_LEAF = 0, T_STEP = 0;

/* ------------------------------------------------------------------ values */
typedef struct { arb_t b; fmpq_t q; int isq; int unk; } Val;
typedef struct { Val *v; long n, cap; double *key; long *ord; long m; long *wide; long nwide; long *zl; long nzl; double zrad;
                 long *u; long nu; int dready; double *dlo, *dhi, *llo, *lhi;
                 short *pa; int *pi, *pj; signed char *pf; long pcap; } VSet;   /* p*: provenance of complete-table values */   /* u: distinct values; d*: outward double bounds */

static void val_init(Val *x) { arb_init(x->b); fmpq_init(x->q); x->isq = 0; x->unk = 0; }
static void val_clear(Val *x) { arb_clear(x->b); fmpq_clear(x->q); }
static void val_set(Val *d, const Val *s) { arb_set(d->b, s->b); fmpq_set(d->q, s->q); d->isq = s->isq; d->unk = s->unk; }
static void val_set_q(Val *x, const fmpq_t q) {
    fmpq_set(x->q, q); x->isq = 1; x->unk = 0; arb_set_fmpq(x->b, q, PREC);
}
static int val_is_zero(const Val *x) { return x->isq && fmpq_is_zero(x->q); }
static int val_is_one(const Val *x) { return x->isq && fmpq_is_one(x->q); }
static int q_ok(const fmpq_t q) {
    return fmpz_bits(fmpq_numref(q)) <= QBITS && fmpz_bits(fmpq_denref(q)) <= QBITS;
}
/* in range: 0 or 1e-300 <= |x| <= 1e300.  Returns 0 only if CERTAINLY out of range. */
static int in_range(const Val *x) {
    if (x->unk) return 1;
    if (val_is_zero(x)) return 1;
    arb_t a; arb_init(a); arb_abs(a, x->b);
    int out = 0;
    if (arb_is_finite(a)) {
        if (arb_gt(a, HI_R)) out = 1;
        else if (x->isq || !arb_contains_zero(a)) { if (arb_lt(a, LO_R)) out = 1; }
    }
    arb_clear(a);
    return !out;
}
static void vset_push(VSet *s, const Val *x) {
    if (s->n >= s->cap) {
        long nc = s->cap * 2 + 256;
        s->v = realloc(s->v, nc * sizeof(Val));
        for (long i = s->cap; i < nc; i++) val_init(&s->v[i]);
        s->cap = nc;
    }
    val_set(&s->v[s->n++], x);
}

/* ------------------------------------------------------------------ operations (values are >= 0) */
/* result flags: returns 0 = undefined, 1 = ok (r set) */
static int op_add(Val *r, const Val *a, const Val *b) {
    if (a->isq && b->isq) { fmpq_t t; fmpq_init(t); fmpq_add(t, a->q, b->q); if (q_ok(t)) { val_set_q(r, t); fmpq_clear(t); return 1; } fmpq_clear(t); }
    arb_add(r->b, a->b, b->b, PREC); r->isq = 0; r->unk = a->unk || b->unk; return 1;
}
static int op_absdiff(Val *r, const Val *a, const Val *b) {
    if (a->isq && b->isq) { fmpq_t t; fmpq_init(t); fmpq_sub(t, a->q, b->q); fmpq_abs(t, t); if (q_ok(t)) { val_set_q(r, t); fmpq_clear(t); return 1; } fmpq_clear(t); }
    arb_sub(r->b, a->b, b->b, PREC); arb_abs(r->b, r->b); r->isq = 0;
    r->unk = a->unk || b->unk;      /* a ball containing 0 is still a valid enclosure (possibly an exact 0) */
    return 1;
}
static int op_mul(Val *r, const Val *a, const Val *b) {
    if (val_is_zero(a) || val_is_zero(b)) { fmpq_t z; fmpq_init(z); val_set_q(r, z); fmpq_clear(z); return 1; }
    if (a->isq && b->isq) { fmpq_t t; fmpq_init(t); fmpq_mul(t, a->q, b->q); if (q_ok(t)) { val_set_q(r, t); fmpq_clear(t); return 1; } fmpq_clear(t); }
    arb_mul(r->b, a->b, b->b, PREC); r->isq = 0; r->unk = a->unk || b->unk; return 1;
}
static int op_div(Val *r, const Val *a, const Val *b) {
    if (val_is_zero(b)) return 0;
    if (val_is_zero(a)) { fmpq_t z; fmpq_init(z); val_set_q(r, z); fmpq_clear(z); return 1; }
    if (a->isq && b->isq) { fmpq_t t; fmpq_init(t); fmpq_div(t, a->q, b->q); if (q_ok(t)) { val_set_q(r, t); fmpq_clear(t); return 1; } fmpq_clear(t); }
    if (b->unk || arb_contains_zero(b->b)) {
        /* divisor in [0, rb] (exact 0 = undefined): the quotient is >= a/rb; out of range if that exceeds 1e300 */
        if (!b->unk && !a->unk && arb_is_positive(a->b)) {
            arb_t t, ub; arb_init(t); arb_init(ub);
            arf_t u; arf_init(u); arb_get_ubound_arf(u, b->b, PREC); arb_set_arf(ub, u); arf_clear(u);
            arb_div(t, a->b, ub, PREC);
            int out = arb_gt(t, HI_R);
            arb_clear(t); arb_clear(ub);
            if (out) { r->unk = -1; r->isq = 0; return 1; }        /* certainly out of range (or undefined) */
        }
        arb_zero_pm_inf(r->b); r->isq = 0; r->unk = 1; N_UNKNOWN++; return 1;
    }
    arb_div(r->b, a->b, b->b, PREC); r->isq = 0; r->unk = a->unk; return 1;
}
/* exact rational power (p/q)^(e) for rational e, when the result is rational and not huge */
static int q_pow(fmpq_t r, const fmpq_t x, const fmpq_t e) {
    if (fmpz_cmp_si(fmpq_denref(e), 64) > 0 || fmpz_bits(fmpq_numref(e)) > 13) return 0;
    slong d = fmpz_get_si(fmpq_denref(e)), k = fmpz_get_si(fmpq_numref(e));
    fmpz_t rn, rd, t; fmpz_init(rn); fmpz_init(rd); fmpz_init(t);
    int ok = 1;
    if (d > 1) {
        if (fmpz_sgn(fmpq_numref(x)) < 0) ok = 0;
        if (ok) { fmpz_root(rn, fmpq_numref(x), d); fmpz_pow_ui(t, rn, d); if (!fmpz_equal(t, fmpq_numref(x))) ok = 0; }
        if (ok) { fmpz_root(rd, fmpq_denref(x), d); fmpz_pow_ui(t, rd, d); if (!fmpz_equal(t, fmpq_denref(x))) ok = 0; }
    } else { fmpz_set(rn, fmpq_numref(x)); fmpz_set(rd, fmpq_denref(x)); }
    slong ak = k < 0 ? -k : k;
    if (ok && (fmpz_bits(rn) + fmpz_bits(rd)) * (ulong)ak > QBITS) ok = 0;
    if (ok) {
        fmpz_pow_ui(rn, rn, ak); fmpz_pow_ui(rd, rd, ak);
        if (k < 0) { if (fmpz_is_zero(rn)) ok = 0; else fmpq_set_fmpz_frac(r, rd, rn); }
        else fmpq_set_fmpz_frac(r, rn, rd);
    }
    fmpz_clear(rn); fmpz_clear(rd); fmpz_clear(t);
    return ok;
}
/* x^(sg*y), x >= 0 */
static int op_pow(Val *r, const Val *x, const Val *y, int sg) {
    if (val_is_zero(y)) { if (val_is_zero(x)) return 0; fmpq_t o; fmpq_init(o); fmpq_one(o); val_set_q(r, o); fmpq_clear(o); return 1; }
    if (val_is_zero(x)) { if (sg < 0) return 0; fmpq_t z; fmpq_init(z); val_set_q(r, z); fmpq_clear(z); return 1; }
    if (val_is_one(x)) { fmpq_t o; fmpq_init(o); fmpq_one(o); val_set_q(r, o); fmpq_clear(o); return 1; }
    if (x->isq && y->isq) {
        fmpq_t e, t; fmpq_init(e); fmpq_init(t); fmpq_set(e, y->q); if (sg < 0) fmpq_neg(e, e);
        int ok = q_pow(t, x->q, e);
        if (ok && q_ok(t)) { val_set_q(r, t); fmpq_clear(e); fmpq_clear(t); return 1; }
        fmpq_clear(e); fmpq_clear(t);
    }
    r->isq = 0;
    if (!x->unk && !y->unk && !arb_is_positive(x->b) && arb_is_finite(y->b)) {
        /* base in [0, rx] (possibly an exact zero): bound x^e from the corner rx^e */
        arb_t e, rx, t1, t2; arb_init(e); arb_init(rx); arb_init(t1); arb_init(t2);
        arb_set(e, y->b); if (sg < 0) arb_neg(e, e);
        arf_t u; arf_init(u); arb_get_ubound_arf(u, x->b, PREC); arb_set_arf(rx, u);
        int done = 0;
        if (arf_sgn(u) > 0 && arb_is_positive(e)) {          /* x^e in [0, max over e of rx^e] */
            arf_t el, eh; arf_init(el); arf_init(eh);
            arb_get_lbound_arf(el, e, PREC); arb_get_ubound_arf(eh, e, PREC);
            arb_t ea; arb_init(ea);
            arb_set_arf(ea, el); arb_pow(t1, rx, ea, PREC);
            arb_set_arf(ea, eh); arb_pow(t2, rx, ea, PREC);
            arb_union(t1, t1, t2, PREC);
            arf_t m; arf_init(m); arb_get_ubound_arf(m, t1, PREC);
            arb_zero(r->b); arb_union(r->b, r->b, t1, PREC); arb_set_arf(t2, m); arb_union(r->b, r->b, t2, PREC);
            arb_zero(t2); arb_union(r->b, r->b, t2, PREC);          /* [0, max] */
            r->unk = arb_is_finite(r->b) ? 0 : 1; if (r->unk) N_UNKNOWN++;
            arf_clear(m); arb_clear(ea); arf_clear(el); arf_clear(eh); done = 1;
        } else if (arf_sgn(u) > 0 && arb_is_negative(e)) {   /* x^e >= min over e of rx^e: out of range if that > 1e300 */
            arf_t el, eh; arf_init(el); arf_init(eh);
            arb_get_lbound_arf(el, e, PREC); arb_get_ubound_arf(eh, e, PREC);
            arb_t ea; arb_init(ea);
            arb_set_arf(ea, el); arb_pow(t1, rx, ea, PREC);
            arb_set_arf(ea, eh); arb_pow(t2, rx, ea, PREC);
            if (arb_gt(t1, HI_R) && arb_gt(t2, HI_R)) { r->unk = -1; done = 1; }
            arb_clear(ea); arf_clear(el); arf_clear(eh);
        }
        arf_clear(u); arb_clear(e); arb_clear(rx); arb_clear(t1); arb_clear(t2);
        if (done) return 1;
    }
    if (x->unk || y->unk || !arb_is_positive(x->b)) { arb_zero_pm_inf(r->b); r->unk = 1; N_UNKNOWN++; return 1; }
    arb_t l, e; arb_init(l); arb_init(e);
    arb_log(l, x->b, PREC);
    arb_mul(e, l, y->b, PREC); if (sg < 0) arb_neg(e, e);
    /* decide the range from y ln x before exponentiating (avoids useless huge balls) */
    arb_t lr; arb_init(lr); arb_set_ui(lr, 10); arb_log(lr, lr, PREC); arb_mul_ui(lr, lr, 300, PREC);
    if (arb_gt(e, lr)) { arb_set(r->b, BIG); arb_mul_ui(r->b, r->b, 10, PREC); r->unk = 0; }        /* certainly > 1e300 */
    else { arb_neg(lr, lr); if (arb_lt(e, lr)) { arb_zero(r->b); arb_add_error_2exp_si(r->b, -2000); r->unk = 0; r->isq = 0; arb_clear(l); arb_clear(e); arb_clear(lr);
               /* certainly < 1e-300 and nonzero: out of range */ r->unk = -1; return 1; }
           arb_exp(r->b, e, PREC); r->unk = 0; }
    arb_clear(l); arb_clear(e); arb_clear(lr);
    return 1;
}
/* in_range variant aware of the "certainly tiny" marker (unk = -1) */
static int keep_val(Val *r) { if (r->unk == -1) return 0; return in_range(r); }

/* ------------------------------------------------------------------ digit sets */
static int popc(int m) { return __builtin_popcount(m); }
static int lowbit(int m) { return m & -m; }
static int is_interval(int m) { if (!m) return 0; int x = m / lowbit(m); return (x & (x + 1)) == 0; }

static VSet FULL[1 << MAXD];
static int FULL_DONE[1 << MAXD];

static void literals(int S, VSet *out) {
    int dg[MAXD], k = 0;
    for (int b = 0; b < ND; b++) if (S >> b & 1) dg[k++] = b + 1;
    /* permutations (any order) or identity (in order) */
    int idx[MAXD]; for (int i = 0; i < k; i++) idx[i] = i;
    int c[MAXD]; memset(c, 0, sizeof c);
    fmpz_t v, p10; fmpz_init(v); fmpz_init(p10);
    fmpq_t q; fmpq_init(q);
    Val x; val_init(&x);
    int i = 0, first = 1;
    for (;;) {
        if (first || !ORD) {
            fmpz_zero(v);
            for (int j = 0; j < k; j++) { fmpz_mul_ui(v, v, 10); fmpz_add_ui(v, v, dg[idx[j]]); }
            fmpz_one(p10);
            for (int j = 0; j <= k; j++) { fmpq_set_fmpz_frac(q, v, p10); val_set_q(&x, q); vset_push(out, &x); fmpz_mul_ui(p10, p10, 10); }
        }
        first = 0;
        if (ORD) break;
        /* Heap's algorithm, iterative */
        while (i < k && c[i] >= i) { c[i] = 0; i++; }
        if (i >= k) break;
        if (i % 2 == 0) { int t = idx[0]; idx[0] = idx[i]; idx[i] = t; } else { int t = idx[c[i]]; idx[c[i]] = idx[i]; idx[i] = t; }
        c[i]++; i = 0;
    }
    fmpz_clear(v); fmpz_clear(p10); fmpq_clear(q); val_clear(&x);
}

/* provenance of a complete-table value: A (0 for a literal), indices i in FULL[A], j in FULL[S ^ A], form code c */
static void prov_set(VSet *o, long k, int A, long i, long j, int c) {
    if (k >= o->pcap) {
        long nc = o->pcap * 2 + 1024;
        o->pa = realloc(o->pa, nc * sizeof(short)); o->pi = realloc(o->pi, nc * sizeof(int));
        o->pj = realloc(o->pj, nc * sizeof(int)); o->pf = realloc(o->pf, nc * sizeof(signed char)); o->pcap = nc;
    }
    o->pa[k] = (short)A; o->pi[k] = (int)i; o->pj[k] = (int)j; o->pf[k] = (signed char)c;
}
/* form c of (s, l): 0 s+l, 1 |s-l|, 2 s*l, 3 s/l, 4 s^l, 5 s^-l, 6 l/s, 7 l^s, 8 l^-s; 0 if undefined */
static int form_op(Val *r, const Val *s, const Val *l, int c) {
    switch (c) {
    case 0: return op_add(r, s, l);
    case 1: return op_absdiff(r, s, l);
    case 2: return op_mul(r, s, l);
    case 3: return op_div(r, s, l);
    case 4: return op_pow(r, s, l, 1);
    case 5: return !val_is_zero(l) && op_pow(r, s, l, -1);
    case 6: return op_div(r, l, s);
    case 7: return op_pow(r, l, s, 1);
    case 8: return !val_is_zero(s) && op_pow(r, l, s, -1);
    }
    return 0;
}
/* forms of Lemma 3: side 1 = s is the left block (in order), side -1 = any order */
static void forms(const Val *s, const Val *l, int side, VSet *out, int A, long i, long j) {
    Val r; val_init(&r);
    int nc = side != 1 ? 9 : 6;
    for (int c = 0; c < nc; c++)
        if (form_op(&r, s, l, c) && keep_val(&r)) { prov_set(out, out->n, A, i, j, c); vset_push(out, &r); }
    val_clear(&r);
}

static double *CMP_KEY; static Val *CMP_V;
static int cmp_ord(const void *a, const void *b) {
    double x = CMP_KEY[*(const long *)a], y = CMP_KEY[*(const long *)b];
    return (x > y) - (x < y);
}
static int cmp_val(const void *A, const void *B) {
    const Val *x = &CMP_V[*(const long *)A], *y = &CMP_V[*(const long *)B];
    double a = CMP_KEY[*(const long *)A], b = CMP_KEY[*(const long *)B];
    if (a != b) return (a > b) - (a < b);
    return 0;
}
static void index_set(VSet *s) {
    double *mid = malloc((s->n + 1) * sizeof(double));
    s->ord = malloc((s->n + 1) * sizeof(long)); s->wide = malloc((s->n + 1) * sizeof(long)); s->zl = malloc((s->n + 1) * sizeof(long));
    s->m = 0; s->nwide = 0; s->nzl = 0; s->zrad = 0;
    long *all = malloc((s->n + 1) * sizeof(long));
    for (long i = 0; i < s->n; i++) { mid[i] = arf_get_d(arb_midref(s->v[i].b), ARF_RND_NEAR); all[i] = i; }
    /* remove exact duplicates (same rational, or identical ball) */
    CMP_KEY = mid; CMP_V = s->v; qsort(all, s->n, sizeof(long), cmp_val);
    long nu = 0;
    for (long i = 0; i < s->n; i++) {
        long k = all[i]; int dup = 0;
        for (long j = nu - 1; j >= 0 && mid[all[j]] == mid[k]; j--) {
            Val *x = &s->v[all[j]], *y = &s->v[k];
            /* with upgrades (UPG) only equal rationals are merged: two different exact values may share a ball */
            if (x->unk == y->unk && x->isq == y->isq && (x->isq ? fmpq_equal(x->q, y->q) : (!UPG && arb_equal(x->b, y->b)))) { dup = 1; break; }
        }
        if (!dup) all[nu++] = k;
    }
    for (long t = 0; t < nu; t++) {
        long i = all[t]; Val *x = &s->v[i]; double md = mid[i];
        if (x->unk || !arb_is_finite(x->b)) { s->wide[s->nwide++] = i; continue; }
        if (!x->isq && arb_contains_zero(x->b)) {       /* possibly an exact zero: kept apart, skipped when lo > zrad */
            double r = mag_get_d(arb_radref(x->b)) + fabs(md); if (r > s->zrad) s->zrad = r;
            s->zl[s->nzl++] = i; continue;
        }
        double rad = x->isq ? 0 : mag_get_d(arb_radref(x->b));
        if (!(rad <= fabs(md) * 1e-12)) { s->wide[s->nwide++] = i; continue; }
        s->ord[s->m++] = i;                 /* already in key order */
    }
    s->key = malloc((s->m + 1) * sizeof(double));
    for (long i = 0; i < s->m; i++) s->key[i] = mid[s->ord[i]];
    s->u = all; s->nu = nu;                 /* distinct values (all[0..nu) after the duplicate removal) */
    free(mid);
}

static void full_(int S);
static int FULL_DEPTH = 0;
static slong TPREC = 0;
static long long N_TABLE_UNK = 0;                    /* precision of the complete tables (0: the working precision) */
static void full(int S) {
    if (FULL_DONE[S]) return;
    double t0 = now(); slong P = PREC; if (TPREC) PREC = TPREC;
    long long U = N_UNKNOWN;               /* unknown table values are flagged in the values themselves */
    FULL_DEPTH++; full_(S); FULL_DEPTH--;
    PREC = P; if (!FULL_DEPTH) { T_FULL += now() - t0; N_TABLE_UNK += N_UNKNOWN - U; N_UNKNOWN = U; }
}
static void full_(int S) {
    if (FULL_DONE[S]) return;
    VSet *out = &FULL[S];
    literals(S, out);
    for (long k = 0; k < out->n; k++) prov_set(out, k, 0, 0, 0, -1);
    int n = popc(S);
    if (n > 1) {
        for (int A = (S - 1) & S; A; A = (A - 1) & S) {
            int B = S ^ A;
            int side;
            if (ORD) { if (!is_interval(A) || !is_interval(B) || !(A & lowbit(S))) continue; side = 1; }
            else { int pa = popc(A); if (2 * pa > n || (2 * pa == n && A > B)) continue; side = -1; }
            full(A); full(B);
            VSet *X = &FULL[A], *Y = &FULL[B];
            for (long i = 0; i < X->n; i++) for (long j = 0; j < Y->n; j++) forms(&X->v[i], &Y->v[j], side, out, A, i, j);
        }
    }
    index_set(out);
    FULL_DONE[S] = 1;
}

/* ------------------------------------------------------------------ upgrade of table values
   The complete tables may be built at a lower precision TPREC (memory).  A table value that is used is
   recomputed at the working precision from its provenance (the same exact operations on the same exact
   operands): a ball at any precision encloses the same exact value (C4, Lemma 15). */
#define UP_SLOTS (1 << 21)
typedef struct { unsigned long key; Val v; int used; } UPEnt;
static UPEnt *UPC;
static long long UP_HIT = 0, UP_MISS = 0;
static void upgrade_into(int S, long idx, Val *r) {
    VSet *X = &FULL[S];
    Val *x = &X->v[idx];
    if (x->isq) { val_set_q(r, x->q); return; }
    if (!UPC) { UPC = calloc(UP_SLOTS, sizeof(UPEnt)); }
    unsigned long key = ((unsigned long)S << 40) | (unsigned long)idx;
    unsigned long h = (key * 0x9E3779B97F4A7C15UL) >> 43;
    UPEnt *e = &UPC[h & (UP_SLOTS - 1)];
    if (e->used && e->key == key) { UP_HIT++; val_set(r, &e->v); return; }
    UP_MISS++;
    int A = X->pa[idx];
    int ok = 0;
    if (A) {
        Val a, b; val_init(&a); val_init(&b);
        upgrade_into(A, X->pi[idx], &a); upgrade_into(S ^ A, X->pj[idx], &b);
        ok = form_op(r, &a, &b, X->pf[idx]);
        val_clear(&a); val_clear(&b);
    }
    if (!ok) val_set(r, x);                    /* cannot happen (same exact operation); the table ball is valid */
    if (!e->used) { val_init(&e->v); e->used = 1; }
    e->key = key; val_set(&e->v, r);
}

/* ------------------------------------------------------------------ windows */
static int may_lie(const Val *x, const arb_t lo, const arb_t hi) {
    if (x->unk) return 1;
    if (arb_lt(x->b, lo)) return 0;
    if (arb_gt(x->b, hi)) return 0;
    return 1;
}
typedef struct { arb_t lo, hi; } Win;
static void win_init(Win *w) { arb_init(w->lo); arb_init(w->hi); }
static void win_clear(Win *w) { arb_clear(w->lo); arb_clear(w->hi); }

/* clamp to [0, BIG]; returns 0 if the window is certainly empty */
static int mkwin(Win *w, const arb_t a, const arb_t b) {
    if (arb_is_positive(a)) arb_set(w->lo, a); else arb_zero(w->lo);
    if (arb_lt(b, BIG)) arb_set(w->hi, b); else arb_set(w->hi, BIG);
    if (!arb_is_finite(w->lo)) arb_zero(w->lo);
    if (!arb_is_finite(w->hi)) arb_set(w->hi, BIG);
    return !arb_gt(w->lo, w->hi);
}

/* preimage windows for l (as motor B); returns count; *cst = 1 if the family is constant here */
static int windows(const char *fam, const Val *s, const arb_t lo, const arb_t hi, Win *w, int *cst) {
    int n = 0; *cst = 0;
    arb_t a, b, t1, t2; arb_init(a); arb_init(b); arb_init(t1); arb_init(t2);
    if (!strcmp(fam, "add")) {
        arb_sub(a, lo, s->b, PREC); arb_sub(b, hi, s->b, PREC); if (mkwin(&w[n], a, b)) n++;
    } else if (!strcmp(fam, "sub")) {
        arb_add(a, s->b, lo, PREC); arb_add(b, s->b, hi, PREC); if (mkwin(&w[n], a, b)) n++;
        arb_sub(a, s->b, hi, PREC); arb_sub(b, s->b, lo, PREC); if (mkwin(&w[n], a, b)) n++;
    } else if (!strcmp(fam, "mul")) {
        if (val_is_zero(s)) *cst = 1;
        else { arb_div(a, lo, s->b, PREC); arb_div(b, hi, s->b, PREC); if (mkwin(&w[n], a, b)) n++; }
    } else if (!strcmp(fam, "sl")) {
        if (val_is_zero(s)) *cst = 1;
        else { arb_div(a, s->b, hi, PREC); if (arb_is_positive(lo)) arb_div(b, s->b, lo, PREC); else arb_set(b, BIG); if (mkwin(&w[n], a, b)) n++; }
    } else if (!strcmp(fam, "ls")) {
        if (!val_is_zero(s)) { arb_mul(a, lo, s->b, PREC); arb_mul(b, hi, s->b, PREC); if (mkwin(&w[n], a, b)) n++; }
    } else if (!strcmp(fam, "ps")) {
        if (val_is_zero(s) || val_is_one(s)) *cst = 1;
        else {
            arb_t lam, A, Bh; arb_init(lam); arb_init(A); arb_init(Bh);
            if (s->unk || !arb_is_positive(s->b)) { arb_zero(a); arb_set(b, BIG); if (mkwin(&w[n], a, b)) n++; }
            else {
                arb_log(lam, s->b, PREC);
                if (arb_is_positive(lo)) arb_log(A, lo, PREC); else arb_neg(A, BIG);
                arb_log(Bh, hi, PREC);
                arb_div(t1, A, lam, PREC); arb_div(t2, Bh, lam, PREC);
                if (!arb_is_finite(t1) || !arb_is_finite(t2)) { arb_zero(a); arb_set(b, BIG); if (mkwin(&w[n], a, b)) n++; }
                else {
                    arf_t lo1, lo2, hi1, hi2; arf_init(lo1); arf_init(lo2); arf_init(hi1); arf_init(hi2);
                    arb_get_lbound_arf(lo1, t1, PREC); arb_get_lbound_arf(lo2, t2, PREC);
                    arb_get_ubound_arf(hi1, t1, PREC); arb_get_ubound_arf(hi2, t2, PREC);
                    arf_t tlo, thi; arf_init(tlo); arf_init(thi);
                    arf_min(tlo, lo1, lo2); arf_max(thi, hi1, hi2);
                    if (arf_sgn(tlo) > 0) { arb_set_arf(a, tlo); arb_set_arf(b, thi); }
                    else if (arf_sgn(thi) < 0) { arf_neg(thi, thi); arf_neg(tlo, tlo); arb_set_arf(a, thi); arb_set_arf(b, tlo); }
                    else { arf_neg(tlo, tlo); arf_max(thi, thi, tlo); arb_zero(a); arb_set_arf(b, thi); }
                    if (mkwin(&w[n], a, b)) n++;
                    arf_clear(lo1); arf_clear(lo2); arf_clear(hi1); arf_clear(hi2); arf_clear(tlo); arf_clear(thi);
                }
            }
            arb_clear(lam); arb_clear(A); arb_clear(Bh);
        }
    } else if (!strcmp(fam, "pl")) {
        if (val_is_zero(s)) *cst = 1;
        else if (s->unk || !arb_is_positive(s->b)) { arb_zero(a); arb_set(b, BIG); if (mkwin(&w[n], a, b)) n++; }
        else {
            arb_t inv, ll, lh; arb_init(inv); arb_init(ll); arb_init(lh);
            arb_inv(inv, s->b, PREC);
            arb_log(lh, hi, PREC);
            int lpos = arb_is_positive(lo);
            if (lpos) arb_log(ll, lo, PREC);
            /* l^s = x: l = x^(1/s) */
            if (lpos) { arb_mul(a, ll, inv, PREC); arb_exp(a, a, PREC); } else arb_zero(a);
            arb_mul(b, lh, inv, PREC); arb_exp(b, b, PREC);
            if (mkwin(&w[n], a, b)) n++;
            /* l^-s = x: l = x^(-1/s) */
            arb_neg(inv, inv);
            arb_mul(a, lh, inv, PREC); arb_exp(a, a, PREC);
            if (lpos) { arb_mul(b, ll, inv, PREC); arb_exp(b, b, PREC); } else arb_set(b, BIG);
            if (mkwin(&w[n], a, b)) n++;
            arb_clear(inv); arb_clear(ll); arb_clear(lh);
        }
    }
    arb_clear(a); arb_clear(b); arb_clear(t1); arb_clear(t2);
    return n;
}

static void family_values(const char *fam, const Val *s, const Val *l, int side, VSet *out) {
    Val r; val_init(&r);
    if (!strcmp(fam, "add")) { if (op_add(&r, s, l) && keep_val(&r)) vset_push(out, &r); }
    else if (!strcmp(fam, "sub")) { if (op_absdiff(&r, s, l) && keep_val(&r)) vset_push(out, &r); }
    else if (!strcmp(fam, "mul")) { if (op_mul(&r, s, l) && keep_val(&r)) vset_push(out, &r); }
    else if (!strcmp(fam, "sl")) { if (op_div(&r, s, l) && keep_val(&r)) vset_push(out, &r); }
    else if (!strcmp(fam, "ls")) { if (op_div(&r, l, s) && keep_val(&r)) vset_push(out, &r); }
    else if (!strcmp(fam, "ps")) {
        if (op_pow(&r, s, l, 1) && keep_val(&r)) vset_push(out, &r);
        if (!val_is_zero(l) && op_pow(&r, s, l, -1) && keep_val(&r)) vset_push(out, &r);
    } else if (!strcmp(fam, "pl")) {
        if (op_pow(&r, l, s, 1) && keep_val(&r)) vset_push(out, &r);
        if (!val_is_zero(s) && op_pow(&r, l, s, -1) && keep_val(&r)) vset_push(out, &r);
    }
    (void)side;
    val_clear(&r);
}

/* ------------------------------------------------------------------ queries */
static void filter_full(int S, const arb_t lo, const arb_t hi, VSet *out) {
    full(S);
    VSet *X = &FULL[S];
    long m = X->m;
    arf_t t; arf_init(t);
    arb_get_lbound_arf(t, lo, PREC); double flo = arf_get_d(t, ARF_RND_DOWN);
    arb_get_ubound_arf(t, hi, PREC); double fhi = arf_get_d(t, ARF_RND_UP);
    arf_clear(t);
    flo -= fabs(flo) * 1e-9 + 1e-300; fhi += fabs(fhi) * 1e-9 + 1e-300;
    long a = 0, b = m;
    while (a < b) { long c = (a + b) / 2; if (X->key[c] < flo) a = c + 1; else b = c; }
    Val u; val_init(&u);
#define EMIT(IDX) { long ix_ = (IDX); Val *x = &X->v[ix_]; \
        if (may_lie(x, lo, hi)) { if (UPG) { upgrade_into(S, ix_, &u); if (keep_val(&u) && may_lie(&u, lo, hi)) vset_push(out, &u); } \
                                   else vset_push(out, x); } }
    for (long i = a; i < m && X->key[i] <= fhi; i++) EMIT(X->ord[i])
    for (long i = 0; i < X->nwide; i++) EMIT(X->wide[i])
    if (X->nzl) {
        arb_t z; arb_init(z); arb_set_d(z, X->zrad * (1 + 1e-9));
        if (!arb_gt(lo, z))                     /* some possibly-zero value may reach the window */
            for (long i = 0; i < X->nzl; i++) EMIT(X->zl[i])
        arb_clear(z);
    }
#undef EMIT
    val_clear(&u);
}

static const char *FAMS_ALL[7] = {"add", "sub", "mul", "sl", "ls", "ps", "pl"};
static const char *FAMS_L[5] = {"add", "sub", "mul", "sl", "ps"};
static const char *FAMS_R[5] = {"add", "sub", "mul", "ls", "pl"};

static void positive_literal(int S, Val *w) {
    fmpz_t v; fmpz_init(v); fmpq_t q; fmpq_init(q);
    for (int b = 0; b < ND; b++) if (S >> b & 1) { fmpz_mul_ui(v, v, 10); fmpz_add_ui(v, v, b + 1); }
    fmpq_set_fmpz(q, v); val_set_q(w, q);
    fmpz_clear(v); fmpq_clear(q);
}

/* ---- outward-rounded double interval arithmetic for the preimage windows of query() (C4, Lemma 15).
   Every double operation is correctly rounded to nearest (IEEE 754), so the exact result lies within
   half an ulp of the computed one and the next double below / above it is a rigorous bound. Logarithms
   and exponentials are computed with Arb and their bounds rounded outward (floor / ceil). */
static double DBIG;                        /* a double >= 1e300 (the clamp of mkwin) */
static double dn(double x) { return nextafter(x, -INFINITY); }
static double up(double x) { return nextafter(x, INFINITY); }
static void dbounds(const arb_t x, double *l, double *h) {
    arf_t t; arf_init(t);
    arb_get_lbound_arf(t, x, PREC); *l = arf_get_d(t, ARF_RND_FLOOR);
    arb_get_ubound_arf(t, x, PREC); *h = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
}
static void dready(VSet *X) {              /* bounds of the distinct values and of their logarithms */
    if (X->dready) return;
    X->dlo = malloc((X->nu + 1) * sizeof(double)); X->dhi = malloc((X->nu + 1) * sizeof(double));
    X->llo = malloc((X->nu + 1) * sizeof(double)); X->lhi = malloc((X->nu + 1) * sizeof(double));
    arb_t t; arb_init(t);
    for (long k = 0; k < X->nu; k++) {
        Val *x = &X->v[X->u[k]];
        X->dlo[k] = X->dhi[k] = X->llo[k] = X->lhi[k] = NAN;
        if (x->unk || !arb_is_finite(x->b)) continue;
        dbounds(x->b, &X->dlo[k], &X->dhi[k]);
        if (arb_is_positive(x->b)) { arb_log(t, x->b, PREC); dbounds(t, &X->llo[k], &X->lhi[k]); }
    }
    arb_clear(t);
    X->dready = 1;
}
static double dexp_lo(double x) { arb_t t; arb_init(t); arb_set_d(t, x); arb_exp(t, t, 64); double l, h; dbounds(t, &l, &h); arb_clear(t); return l; }
static double dexp_hi(double x) { arb_t t; arb_init(t); arb_set_d(t, x); arb_exp(t, t, 64); double l, h; dbounds(t, &l, &h); arb_clear(t); return h; }
/* product / quotient of intervals, outward */
static void imul(double a1, double a2, double b1, double b2, double *l, double *h) {
    double c[4] = {a1 * b1, a1 * b2, a2 * b1, a2 * b2}; *l = c[0]; *h = c[0];
    for (int i = 1; i < 4; i++) { if (c[i] < *l) *l = c[i]; if (c[i] > *h) *h = c[i]; }
    *l = dn(*l); *h = up(*h);
}
static void idiv(double a1, double a2, double b1, double b2, double *l, double *h) {   /* 0 not in [b1, b2] */
    double c[4] = {a1 / b1, a1 / b2, a2 / b1, a2 / b2}; *l = c[0]; *h = c[0];
    for (int i = 1; i < 4; i++) { if (c[i] < *l) *l = c[i]; if (c[i] > *h) *h = c[i]; }
    *l = dn(*l); *h = up(*h);
}
typedef struct { double a, b; } DWin;
static long long N_ARBWIN = 0;
static int dwin(DWin *w, double a, double b) {             /* clamp to [0, DBIG]; 0 if empty */
    if (!(a > 0)) a = 0;                                     /* also NaN */
    if (!(b < DBIG)) b = DBIG;
    w->a = a; w->b = b;
    return a <= b;
}
/* the windows of windows(), for the k-th distinct value of X, in outward double arithmetic;
   L = [Ll, Lh] bounds of lo, H = [Hl, Hh] bounds of hi, (A1, A2), (B1, B2): bounds of log lo, log hi */
static int dwindows(const char *fam, VSet *X, long k, double Ll, double Hh, int lpos,
                    double A1, double A2, double B1, double B2, DWin *w, int *cst) {
    int n = 0; *cst = 0;
    Val *s = &X->v[X->u[k]];
    double sl = X->dlo[k], sh = X->dhi[k];
    int sok = !s->unk && !isnan(sl) && !isnan(sh);
    if (!strcmp(fam, "add")) {
        if (!sok) { if (dwin(&w[n], 0, DBIG)) n++; }
        else if (dwin(&w[n], dn(Ll - sh), up(Hh - sl))) n++;
    } else if (!strcmp(fam, "sub")) {
        if (!sok) { if (dwin(&w[n], 0, DBIG)) n++; }
        else {
            if (dwin(&w[n], dn(sl + Ll), up(sh + Hh))) n++;
            if (dwin(&w[n], dn(sl - Hh), up(sh - Ll))) n++;
        }
    } else if (!strcmp(fam, "mul")) {
        if (val_is_zero(s)) *cst = 1;
        else if (!sok || !(sl > 0)) { if (dwin(&w[n], 0, DBIG)) n++; }
        else if (dwin(&w[n], dn(Ll / sh), up(Hh / sl))) n++;
    } else if (!strcmp(fam, "sl")) {
        if (val_is_zero(s)) *cst = 1;
        else if (!sok) { if (dwin(&w[n], 0, DBIG)) n++; }
        else if (dwin(&w[n], Hh > 0 ? dn(sl / Hh) : 0, Ll > 0 ? up(sh / Ll) : DBIG)) n++;
    } else if (!strcmp(fam, "ls")) {
        if (!val_is_zero(s)) {
            if (!sok) { if (dwin(&w[n], 0, DBIG)) n++; }
            else if (dwin(&w[n], dn(Ll * sl), up(Hh * sh))) n++;
        }
    } else if (!strcmp(fam, "ps")) {
        if (val_is_zero(s) || val_is_one(s)) *cst = 1;
        else {
            double l1 = X->llo[k], l2 = X->lhi[k];
            if (!sok || isnan(l1) || isnan(l2) || (l1 <= 0 && l2 >= 0)) { if (dwin(&w[n], 0, DBIG)) n++; }
            else {
                double t1l, t1h, t2l, t2h;
                idiv(A1, A2, l1, l2, &t1l, &t1h); idiv(B1, B2, l1, l2, &t2l, &t2h);
                double tlo = t1l < t2l ? t1l : t2l, thi = t1h > t2h ? t1h : t2h;
                if (isnan(tlo) || isnan(thi)) { if (dwin(&w[n], 0, DBIG)) n++; }
                else if (tlo > 0) { if (dwin(&w[n], tlo, thi)) n++; }
                else if (thi < 0) { if (dwin(&w[n], -thi, -tlo)) n++; }
                else { double m = thi > -tlo ? thi : -tlo; if (dwin(&w[n], 0, m)) n++; }
            }
        }
    } else if (!strcmp(fam, "pl")) {
        if (val_is_zero(s)) *cst = 1;
        else if (!sok || !(sl > 0)) { if (dwin(&w[n], 0, DBIG)) n++; }
        else {
            double i1 = dn(1.0 / sh), i2 = up(1.0 / sl), pl_, ph_;
            /* l^s = x: l = exp(log x / s) */
            double a = 0, b;
            if (lpos) { imul(A1, A2, i1, i2, &pl_, &ph_); a = dexp_lo(pl_); }
            imul(B1, B2, i1, i2, &pl_, &ph_); b = dexp_hi(ph_);
            if (dwin(&w[n], a, b)) n++;
            /* l^-s = x: l = exp(-log x / s) */
            imul(B1, B2, -i2, -i1, &pl_, &ph_); a = dexp_lo(pl_);
            b = DBIG;
            if (lpos) { imul(A1, A2, -i2, -i1, &pl_, &ph_); b = dexp_hi(ph_); }
            if (dwin(&w[n], a, b)) n++;
        }
    }
    return n;
}

static void query(int S, const arb_t lo, const arb_t hi, VSet *out) {
    if (popc(S) <= 4) { filter_full(S, lo, hi, out); return; }
    VSet L; memset(&L, 0, sizeof L); literals(S, &L);
    for (long i = 0; i < L.n; i++) if (may_lie(&L.v[i], lo, hi)) vset_push(out, &L.v[i]);
    for (long i = 0; i < L.cap; i++) val_clear(&L.v[i]); free(L.v);
    int n = popc(S);
    double Ll, Lh, Hl, Hh, A1 = -DBIG, A2 = -DBIG, B1, B2;     /* bounds of lo, hi, log lo, log hi */
    dbounds(lo, &Ll, &Lh); dbounds(hi, &Hl, &Hh);
    int lpos = arb_is_positive(lo);
    { arb_t t; arb_init(t);
      if (lpos) { arb_log(t, lo, PREC); dbounds(t, &A1, &A2); }
      arb_log(t, hi, PREC); dbounds(t, &B1, &B2);
      if (!arb_is_positive(hi)) { B1 = -DBIG; B2 = -DBIG; }
      arb_clear(t); }
    if (isnan(Ll) || Ll < 0) Ll = 0;
    if (isnan(Hh) || !(Hh < DBIG)) Hh = DBIG;
    for (int A = (S - 1) & S; A; A = (A - 1) & S) {
        int B = S ^ A, S1, S2, side;
        if (ORD) {
            if (!is_interval(A) || !is_interval(B) || !(A & lowbit(S))) continue;
            if (popc(A) <= popc(B)) { S1 = A; S2 = B; side = 1; } else { S1 = B; S2 = A; side = 0; }
        } else { int pa = popc(A); if (2 * pa > n || (2 * pa == n && A > B)) continue; S1 = A; S2 = B; side = -1; }
        full(S1);
        VSet *X = &FULL[S1];
        dready(X);
        const char **fams = side == 1 ? FAMS_L : side == 0 ? FAMS_R : FAMS_ALL;
        int nf = side == -1 ? 7 : 5;
        Val SU; val_init(&SU);
        for (long k = 0; k < X->nu; k++) {
            Val *s = &X->v[X->u[k]];
            long sidx = X->u[k]; int sup = 0;      /* s at the working precision, computed when first needed */
#define SV() (UPG ? (sup ? &SU : (upgrade_into(S1, sidx, &SU), sup = 1, &SU)) : s)
            for (int f = 0; f < nf; f++) {
                DWin dw[2]; int cst;
                int nw = dwindows(fams[f], X, k, Ll, Hh, lpos, A1, A2, B1, B2, dw, &cst);
                Win w[2]; win_init(&w[0]); win_init(&w[1]);
                int use_arb = 0;          /* double windows lose the information near the ends of the double range
                                             and when s is much larger than the window (cancellation): use Arb */
                for (int j = 0; j < nw; j++) if (dw[j].b < 1e-280 || dw[j].a > 1e280) use_arb = 1;
                if (f < 2 && X->dhi[k] > 1e6 * Hh) use_arb = 1;          /* add, sub */
                if (use_arb) { N_ARBWIN++; nw = windows(fams[f], s, lo, hi, w, &cst); }
                else for (int j = 0; j < nw; j++) { arb_set_d(w[j].lo, dw[j].a); arb_set_d(w[j].hi, dw[j].b); }
                if (cst) {
                    Val wl; val_init(&wl); positive_literal(S2, &wl);
                    VSet tmp; memset(&tmp, 0, sizeof tmp);
                    family_values(fams[f], SV(), &wl, side, &tmp);
                    for (long k = 0; k < tmp.n; k++) if (may_lie(&tmp.v[k], lo, hi)) vset_push(out, &tmp.v[k]);
                    for (long k = 0; k < tmp.cap; k++) val_clear(&tmp.v[k]); free(tmp.v); val_clear(&wl);
                }
                for (int j = 0; j < nw; j++) {
                    VSet sub; memset(&sub, 0, sizeof sub);
                    query(S2, w[j].lo, w[j].hi, &sub);
                    for (long k = 0; k < sub.n; k++) {
                        VSet tmp; memset(&tmp, 0, sizeof tmp);
                        family_values(fams[f], SV(), &sub.v[k], side, &tmp);
                        for (long q = 0; q < tmp.n; q++) if (may_lie(&tmp.v[q], lo, hi)) vset_push(out, &tmp.v[q]);
                        for (long q = 0; q < tmp.cap; q++) val_clear(&tmp.v[q]); free(tmp.v);
                    }
                    for (long k = 0; k < sub.cap; k++) val_clear(&sub.v[k]); free(sub.v);
                }
                win_clear(&w[0]); win_clear(&w[1]);
            }
#undef SV
        }
        val_clear(&SU);
    }
}

/* ------------------------------------------------------------------ class windows */
static double ulp_of(double r) { int e; frexp(r, &e); return ldexp(1.0, e - 53); }
static void class_window(char **tok, int ntok, arb_t lo, arb_t hi) {
    /* tok: kind, then fields */
    if (!strcmp(tok[0], "o")) {
        double v = strtod(tok[1], NULL); int code = atoi(tok[2]);
        if (code == 0) { arb_set_d(lo, v); arb_set_d(hi, v); return; }
        double rho = ldexp(ulp_of(v), code - 1);
        double a = nextafter(v - rho, -INFINITY), b = nextafter(v + rho, INFINITY);
        arb_set_d(lo, a < 0 ? 0 : a); arb_set_d(hi, b); return;
    }
    /* anchored: a:p/q:m1:m2   (l: is anchor 1) */
    fmpq_t Q; fmpq_init(Q);
    double m1, m2;
    if (!strcmp(tok[0], "l")) { fmpq_one(Q); m1 = strtod(tok[1], NULL); m2 = strtod(tok[2], NULL); }
    else { fmpq_set_str(Q, tok[1], 10); m1 = strtod(tok[2], NULL); m2 = strtod(tok[3], NULL); }
    arb_t q, e; arb_init(q); arb_init(e);
    arb_set_fmpq(q, Q, PREC);
    arb_set_d(e, m1); arb_exp(e, e, PREC); arb_mul(e, e, q, PREC);
    arf_t t; arf_init(t); arb_get_lbound_arf(t, e, PREC); arb_set_arf(lo, t);
    arb_set_d(e, m2); arb_exp(e, e, PREC); arb_mul(e, e, q, PREC);
    arb_get_ubound_arf(t, e, PREC); arb_set_arf(hi, t);
    arf_clear(t); arb_clear(q); arb_clear(e); fmpq_clear(Q);
    (void)ntok;
}

/* ------------------------------------------------------------------ path verification */
#define MAXTOK 16
typedef struct { int mask; char what[8]; char fam[8]; int side; char *spec[6]; int nspec; } Tok;

static int split(char *s, char sep, char **parts, int maxp) {
    int n = 0; parts[n++] = s;
    for (char *p = s; *p; p++) if (*p == sep && n < maxp) { *p = 0; parts[n++] = p + 1; }
    return n;
}

/* The values of S1 in a step's class window depend only on the step token (mask and class), not on the
   enclosing window: they are computed once and kept in a cache (results are identical to recomputing). */
static long long UNK0 = 0;
/* copy a set into exactly sized storage (cached sets) */
static void vset_compact(VSet *d, VSet *t) {
    memset(d, 0, sizeof *d);
    d->v = malloc((t->n + 1) * sizeof(Val)); d->cap = t->n; d->n = t->n;
    for (long i = 0; i < t->n; i++) { val_init(&d->v[i]); val_set(&d->v[i], &t->v[i]); }
    t->n = 0;
}
static VSet QTMP;
#define SQ_SLOTS 65536
typedef struct { char *key; VSet v; long long unk; } SQEnt;
static SQEnt SQ[SQ_SLOTS];
static long long SQ_HIT = 0, SQ_MISS = 0;
static VSet *step_query(Tok *st) {
    char key[512]; int o = snprintf(key, sizeof key, "%d", st->mask);
    for (int i = 0; i < st->nspec && o < (int)sizeof key - 2; i++) o += snprintf(key + o, sizeof key - o, ":%s", st->spec[i]);
    unsigned long h = 1469598103934665603UL; for (char *p = key; *p; p++) h = (h ^ (unsigned char)*p) * 1099511628211UL;
    SQEnt *e = &SQ[h % SQ_SLOTS];
    if (e->key && !strcmp(e->key, key)) { SQ_HIT++; N_UNKNOWN += e->unk; return &e->v; }
    SQ_MISS++;
    if (e->key) { free(e->key); for (long i = 0; i < e->v.cap; i++) val_clear(&e->v.v[i]); free(e->v.v); }
    e->key = strdup(key); memset(&e->v, 0, sizeof e->v);
    arb_t clo, chi; arb_init(clo); arb_init(chi);
    class_window(st->spec, st->nspec, clo, chi);
    long long u0 = N_UNKNOWN; double t0 = now(), f0 = T_FULL;
    QTMP.n = 0; query(st->mask, clo, chi, &QTMP); vset_compact(&e->v, &QTMP);
    T_STEP += now() - t0 - (T_FULL - f0);
    e->unk = N_UNKNOWN - u0;            /* unknown events are replayed on every hit (keep decision unchanged) */
    arb_clear(clo); arb_clear(chi);
    return &e->v;
}

#define LQ_SLOTS (1 << 20)
static SQEnt *LQ;
static long long LQ_DIRECT = 0, LQ_HIT = 0, LQ_MISS = 0, LQ_VALS = 0, LQ_FLUSH = 0;
static long long LQ_MAXVALS = 4000000;      /* about 0.4 GB of cached values */
static VSet *leaf_query(int digs, Tok *lf, const arb_t clo, const arb_t chi) {
    if (!LQ) LQ = calloc(LQ_SLOTS, sizeof(SQEnt));
    char key[512]; int o = snprintf(key, sizeof key, "%d:%s", digs, lf->what);
    for (int i = 0; i < lf->nspec && o < (int)sizeof key - 2; i++) o += snprintf(key + o, sizeof key - o, ":%s", lf->spec[i]);
    unsigned long h = 1469598103934665603UL; for (char *p = key; *p; p++) h = (h ^ (unsigned char)*p) * 1099511628211UL;
    SQEnt *e = &LQ[h % LQ_SLOTS];
    if (e->key && !strcmp(e->key, key)) { LQ_HIT++; N_UNKNOWN += e->unk; return &e->v; }
    LQ_MISS++;
    if (e->key) { LQ_VALS -= e->v.cap; free(e->key); for (long i = 0; i < e->v.cap; i++) val_clear(&e->v.v[i]); free(e->v.v); e->key = NULL; }
    if (LQ_VALS > LQ_MAXVALS) {                      /* memory bound: empty the whole cache */
        for (long i = 0; i < LQ_SLOTS; i++) if (LQ[i].key) { free(LQ[i].key); for (long j = 0; j < LQ[i].v.cap; j++) val_clear(&LQ[i].v.v[j]); free(LQ[i].v.v); LQ[i].key = NULL; }
        LQ_VALS = 0; LQ_FLUSH++;
    }
    e->key = strdup(key); memset(&e->v, 0, sizeof e->v);
    long long u0 = N_UNKNOWN; double t0 = now(), f0 = T_FULL;
    QTMP.n = 0; query(digs, clo, chi, &QTMP); vset_compact(&e->v, &QTMP);
    T_LEAF += now() - t0 - (T_FULL - f0);
    e->unk = N_UNKNOWN - u0; LQ_VALS += e->v.cap;
    return &e->v;
}

/* all values covered by a (sub)path within [lo, hi] (superset) */
static void collect_verify(int digs, Tok *steps, int nsteps, Tok *leaf, const arb_t lo, const arb_t hi, VSet *out, int depth) {
    if (nsteps == 0) {
        arb_t clo, chi, a, b; arb_init(clo); arb_init(chi); arb_init(a); arb_init(b);
        class_window(leaf->spec, leaf->nspec, clo, chi);
        if (!strcmp(leaf->what, "lit")) {
            VSet L; memset(&L, 0, sizeof L); literals(digs, &L);
            for (long k = 0; k < L.n; k++) if (may_lie(&L.v[k], clo, chi) && may_lie(&L.v[k], lo, hi)) vset_push(out, &L.v[k]);
            for (long k = 0; k < L.cap; k++) val_clear(&L.v[k]); free(L.v);
        } else {
            /* superset of the values of digs in [clo, chi] (cached per leaf token), then restricted to [lo, hi]:
               every value in the intersection lies in both windows, so it is kept (Lemma 15 unchanged) */
            if (arb_gt(clo, lo)) arb_set(a, clo); else arb_set(a, lo);
            if (arb_lt(chi, hi)) arb_set(b, chi); else arb_set(b, hi);
            if (!arb_gt(a, b)) {
                long long u0 = N_UNKNOWN;
                VSet *XL = leaf_query(digs, leaf, clo, chi);
                int has_unk = N_UNKNOWN != u0;
                for (long k = 0; k < XL->n && !has_unk; k++) if (XL->v[k].unk) has_unk = 1;
                if (has_unk) {                  /* undecided values: query the intersection directly (tighter) */
                    N_UNKNOWN = u0; LQ_DIRECT++;
                    query(digs, a, b, out);
                } else
                    for (long k = 0; k < XL->n; k++) if (may_lie(&XL->v[k], lo, hi)) vset_push(out, &XL->v[k]);
            }
        }
        arb_clear(clo); arb_clear(chi); arb_clear(a); arb_clear(b);
        return;
    }
    Tok *st = &steps[0];
    int S1 = st->mask, S2 = digs & ~S1;
    VSet *XP = step_query(st);
    for (long i = 0; i < XP->n; i++) {
        if (depth == 0 && (out->n > 0 || N_UNKNOWN != UNK0)) break;   /* path already kept: stop early */
        Val *s = &XP->v[i];
        if (nsteps == 1 && !strcmp(leaf->what, "const")) {
            Val wl; val_init(&wl); positive_literal(S2, &wl);
            VSet tmp; memset(&tmp, 0, sizeof tmp);
            family_values(st->fam, s, &wl, st->side, &tmp);
            for (long k = 0; k < tmp.n; k++) if (may_lie(&tmp.v[k], lo, hi)) vset_push(out, &tmp.v[k]);
            for (long k = 0; k < tmp.cap; k++) val_clear(&tmp.v[k]); free(tmp.v); val_clear(&wl);
            continue;
        }
        Win w[2]; win_init(&w[0]); win_init(&w[1]); int cst;
        int nw = windows(st->fam, s, lo, hi, w, &cst);
        for (int j = 0; j < nw; j++) {
            VSet sub; memset(&sub, 0, sizeof sub);
            collect_verify(S2, steps + 1, nsteps - 1, leaf, w[j].lo, w[j].hi, &sub, depth + 1);
            for (long k = 0; k < sub.n; k++) {
                VSet tmp; memset(&tmp, 0, sizeof tmp);
                family_values(st->fam, s, &sub.v[k], st->side, &tmp);
                for (long q = 0; q < tmp.n; q++) if (may_lie(&tmp.v[q], lo, hi)) vset_push(out, &tmp.v[q]);
                for (long q = 0; q < tmp.cap; q++) val_clear(&tmp.v[q]); free(tmp.v);
            }
            for (long k = 0; k < sub.cap; k++) val_clear(&sub.v[k]); free(sub.v);
        }
        win_clear(&w[0]); win_clear(&w[1]);
    }
}

static int parse_tok(char *s, Tok *t) {
    char *p[MAXTOK]; int n = split(s, ':', p, MAXTOK);
    t->mask = atoi(p[0]);
    strncpy(t->what, p[1], 7); t->what[7] = 0;
    if (!strcmp(p[1], "leaf") || !strcmp(p[1], "lit")) { t->nspec = n - 2; for (int i = 0; i < t->nspec; i++) t->spec[i] = p[2 + i]; return 1; }
    if (!strcmp(p[1], "const")) { t->nspec = 0; return 1; }
    strncpy(t->fam, p[1], 7); t->fam[7] = 0; strcpy(t->what, "step");
    t->side = atoi(p[2]);
    t->nspec = n - 3; for (int i = 0; i < t->nspec; i++) t->spec[i] = p[3 + i];
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: motorc n ord window_file candidates_file kept_file [precision_bits]\n"); return 1; }
    double T0 = now();
    ND = atoi(argv[1]); ORD = atoi(argv[2]);
    if (argc > 6) PREC = atol(argv[6]);
    arb_init(LO_R); arb_init(HI_R); arb_init(BIG); arb_init(ZERO_A);
    arb_set_str(LO_R, "1e-300", PREC); arb_set_str(HI_R, "1e300", PREC);
    arb_set(BIG, HI_R);
    DBIG = up(1e300);
    FILE *wf = fopen(argv[3], "r"); char sa[64], sb[64];
    if (fscanf(wf, "%63s %63s", sa, sb) != 2) return 2;
    fclose(wf);
    arb_t lo, hi; arb_init(lo); arb_init(hi);
    arb_set_d(lo, strtod(sa, NULL)); arb_set_d(hi, strtod(sb, NULL));
    int digs = (1 << ND) - 1;
    /* MOTORC_TPREC: precision of the complete tables; MOTORC_JOBS: worker processes sharing the tables */
    if (getenv("MOTORC_TPREC")) TPREC = atol(getenv("MOTORC_TPREC"));
    UPG = TPREC && TPREC < PREC;
    int J = getenv("MOTORC_JOBS") ? atoi(getenv("MOTORC_JOBS")) : 1, me = 0;
    if (J > 1) {
        for (int S = 1; S < (1 << ND); S++) if (popc(S) <= 4 && (!ORD || is_interval(S))) full(S);
        fprintf(stderr, "tables built in %.1f s\n", T_FULL);
        for (me = 0; me < J; me++) { pid_t p = fork(); if (p == 0) break; }
        if (me == J) {                                   /* parent: wait, then join the parts in order */
            int bad = 0, st;
            for (int k = 0; k < J; k++) { if (wait(&st) < 0 || !WIFEXITED(st) || WEXITSTATUS(st)) bad = 1; }
            if (bad) { fprintf(stderr, "a worker failed\n"); return 3; }
            FILE *of = fopen(argv[5], "w"); char pn[4096], *ln = NULL; size_t lc = 0;
            for (int k = 0; k < J; k++) {
                snprintf(pn, sizeof pn, "%s.part%d", argv[5], k); FILE *pf = fopen(pn, "r");
                while (getline(&ln, &lc, pf) > 0) fputs(ln, of);
                fclose(pf); remove(pn);
            }
            fclose(of);
            fprintf(stderr, "all %d workers done, total %.1f s\n", J, now() - T0);
            return 0;
        }
    }
    char partname[4096]; if (J > 1) snprintf(partname, sizeof partname, "%s.part%d", argv[5], me); else snprintf(partname, sizeof partname, "%s", argv[5]);
    FILE *cf = fopen(argv[4], "r"), *of = fopen(partname, "w");
    char *line = NULL; size_t cap = 0; long nl = 0, kept = 0, idx = -1;
    char buf[8192];
    while (getline(&line, &cap, cf) > 0) {
        size_t L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (!L) continue;
        if (++idx % J != me) continue;                  /* this worker's share */
        strncpy(buf, line, sizeof buf - 1); buf[sizeof buf - 1] = 0;
        char *tk[MAXTOK]; int nt = split(buf, ' ', tk, MAXTOK);
        int first = (!strcmp(tk[1], "C")) ? 2 : 1;      /* "W0 C step ... leaf" */
        Tok steps[MAXTOK], leaf;
        int ns = 0;
        for (int i = first; i < nt - 1; i++) parse_tok(tk[i], &steps[ns++]);
        parse_tok(tk[nt - 1], &leaf);
        long unk_before = N_UNKNOWN;
        VSet out; memset(&out, 0, sizeof out);
        UNK0 = N_UNKNOWN;
        collect_verify(digs, steps, ns, &leaf, lo, hi, &out, 0);
        int keep = out.n > 0;
        for (long i = 0; i < out.cap; i++) val_clear(&out.v[i]); free(out.v);
        if (N_UNKNOWN != unk_before) keep = 1;            /* anything undecided: motor B decides */
        nl++;
        if (keep) { fprintf(of, "%s\n", line); kept++; }
        if (nl % 10000 == 0) { fprintf(stderr, "[%d] %ld paths, %ld kept, %lld unknown\n", me, nl, kept, N_UNKNOWN); fflush(of); }
    }
    fclose(cf); fclose(of);
    fprintf(stderr, "[%d] time (s): tables %.1f, step queries %.1f, leaf queries %.1f, total %.1f\n", me, T_FULL, T_STEP, T_LEAF, now() - T0);
    fprintf(stderr, "[%d] done: %ld paths, %ld kept, %lld unknown events; step cache %lld/%lld, leaf cache %lld/%lld (hits/misses) %lld direct, upgrades %lld/%lld, unknown table values %lld, Arb windows %lld\n", me, nl, kept, N_UNKNOWN, SQ_HIT, SQ_MISS, LQ_HIT, LQ_MISS, LQ_DIRECT, UP_HIT, UP_MISS, N_TABLE_UNK, N_ARBWIN);
    return 0;
}
