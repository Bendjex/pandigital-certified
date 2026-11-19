#define NO_MAIN
#include "motora_instr.c"
int main(int argc, char **argv) {
    ND = atoi(argv[1]); ORD = atoi(argv[2]); FULL = (1 << ND) - 1; int upto = atoi(argv[4]);
    RLO_SAFE = pred(1e-300); RHI_SAFE = succ(1e300);
    load_res(argv[3]);
    for (int m = 1; m <= FULL; m++) if (mask_ok(m)) Lit[m] = make_literals(m, 1);
    for (int k = 1; k <= upto && k <= ND; k++) {
        int cnt = 0, masks[512];
        for (int m = 1; m <= FULL; m++) if (popc(m) == k && mask_ok(m)) masks[cnt++] = m;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cnt; i++) { if (k == 1) { T[masks[i]] = Lit[masks[i]]; continue; } build_table(masks[i], 1); }
        i64 tot = 0; for (int i = 0; i < cnt; i++) tot += T[masks[i]].n + T[masks[i]].na;
        printf("n=%d ord=%d size %d: %d sets, %lld classes; finish_real calls %lld; with hi=+inf and lo<=RHI: %lld (returned D_OK: %lld); unresolved %lld\n",
               ND, ORD, k, cnt, (long long)tot, N_FR, N_HI_INF, N_HI_INF_OK, N_UNRES);
        fflush(stdout);
    }
    if (N_HI_INF_OK) printf("example: r=%g lo=%g\n", EX_R, EX_LO);
    return 0;
}
