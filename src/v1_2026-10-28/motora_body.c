
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
