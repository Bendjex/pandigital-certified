/* recuento de finish_real con cota superior +inf, usando build_tables de motor A (libera las tablas del tamaño máximo,
   como el modo "stage"): ./probe_inf2 n ord res tam_max */
#define NO_MAIN
#include "motora_instr.c"
int main(int argc, char **argv) {
    ND = atoi(argv[1]); ORD = atoi(argv[2]); FULL = (1 << ND) - 1; int upto = atoi(argv[4]);
    RLO_SAFE = pred(1e-300); RHI_SAFE = succ(1e300);
    load_res(argv[3]);
    FREE_TOP = 1;
    build_tables(upto, 1);
    printf("n=%d ord=%d tablas hasta %d: llamadas a finish_real %lld; con hi=+inf y lo<=succ(1e300): %lld (devolvieron D_OK: %lld); sin resolver %lld\n",
           ND, ORD, upto, N_FR, N_HI_INF, N_HI_INF_OK, N_UNRES);
    if (N_HI_INF_OK) printf("ejemplo: r=%g lo=%g\n", EX_R, EX_LO);
    return 0;
}
