
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
typedef struct { double lo, hi, Llo, Lhi; int hasA; Rat A; double n1, n2; } Win;

typedef struct { char **l; i64 n, cap; } Lines;
static Lines OUT[512];
static long long NWIN[512];

#define SPILL_AT 1000000
static FILE *CUR_OUT; static const char *CUR_PREFIX = ""; static i64 N_SPILLED = 0;
static void emit(const Path *p, const char *leaf) {
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

static void leaf_lookup(int S, const Win *w, double La, double Lb, const Path *p, int is_lit) {
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
                } else if (v == 1 && !(La <= 0 && 0 <= Lb)) continue;
            } else if (h1 > 0.5 && h2 < 2) {             /* near 1: the log window decides */
                if (ln_hi(h2) < La || ln_lo(h1) > Lb) continue;
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
                if (l2 < La || l1 > Lb) continue;
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
    Win w = mkwin(exp_lo(u1), u2 == INFINITY ? INFINITY : exp_hi(u2), allow0);
    w.Llo = u1; w.Lhi = u2;
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
            for (int j = 0; j < nw; j++) { NWIN[omp_get_thread_num()]++; query(S2, &w[j], p); }
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
    *La = W->lo > 0 ? ln_lo(W->lo) : -INFINITY;
    *Lb = W->hi > 0 ? ln_hi(W->hi) : -INFINITY;
    if (W->Llo == W->Llo && W->Llo > *La) *La = W->Llo;
    if (W->Lhi == W->Lhi && W->Lhi < *Lb) *Lb = W->Lhi;
}

static void query(int S, const Win *W, Path *p) {
    if (p->depth >= MAXDEPTH - 1) { fprintf(stderr, "depth overflow\n"); exit(5); }
    double La, Lb; win_logs(W, &La, &Lb);
    if (La > Lb) return;
    if (popc(S) <= TBL) { leaf_lookup(S, W, La, Lb, p, 0); return; }
    leaf_lookup(S, W, La, Lb, p, 1);
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
        leaf_lookup(FULL, &W, La, Lb, &p, 0);
    } else {
        Path p0; p0.depth = 0;
        leaf_lookup(FULL, &W, La, Lb, &p0, 1);
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
        i64 done = 0, step = nitems / 100 > 0 ? nitems / 100 : 1;
        #pragma omp parallel for schedule(dynamic, 16)
        for (i64 j = 0; j < nitems; j++) {
            i64 it = perm[j];
            Path p; p.depth = 0;
            int k = (int)ik[it];
            query_split(s1[k], s2[k], sd[k], &W, La, Lb, &p, ii[it]);
            i64 d;
            #pragma omp atomic capture
            d = ++done;
            if (d % step == 0) {
                double el = omp_get_wtime() - t0;
                #pragma omp critical
                { fprintf(stderr, "%sprogress %lld/%lld (%.0f%%)  %.2f h elapsed, about %.2f h left\n", prefix, (long long)d,
                          (long long)nitems, 100.0 * d / nitems, el / 3600, el / 3600 * ((double)nitems / d - 1)); fflush(stderr); }
            }
        }
        free(perm);
        free(ik); free(ii);
    }
    i64 tot = 0; long long nw = 0;
    for (int t = 0; t < 512; t++) { tot += OUT[t].n; nw += NWIN[t]; }
    i64 spilled = N_SPILLED; N_SPILLED = 0;
    char **all = malloc((tot + 1) * sizeof(char *)); i64 k = 0;
    for (int t = 0; t < 512; t++) for (i64 i = 0; i < OUT[t].n; i++) all[k++] = OUT[t].l[i];
    qsort(all, tot, sizeof(char *), cmp_str);
    for (i64 i = 0; i < tot; i++) { fprintf(out, "%sC %s\n", prefix, all[i]); free(all[i]); }
    free(all);
    for (int t = 0; t < 512; t++) { OUT[t].n = 0; NWIN[t] = 0; }
    fprintf(stderr, "%scandidates: %lld, windows queried: %lld, %.1fs\n", prefix, (long long)(tot + spilled), nw, omp_get_wtime() - t0);
    fflush(out); fflush(stderr);
}

#ifndef NO_MAIN
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: motora stage n ord k res pend | multi n ord res windows out\n"); return 1; }
    ND = atoi(argv[2]); ORD = atoi(argv[3]); FULL = (1 << ND) - 1;
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
