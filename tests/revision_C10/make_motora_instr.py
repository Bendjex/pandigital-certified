"""genera motora_instr.c: motora.c con un contador de las llamadas a finish_real cuya cota superior es +inf
(y lo <= succ(1e300)).  uso: python3 make_motora_instr.py RUTA/src/motora.c > motora_instr.c
compilar: gcc -O2 -ffp-contract=off -fopenmp -I RUTA/src -o probe_inf probe_inf.c RUTA/src/coremath/pow.c -lm
ejecutar: ./probe_inf n ord res_n_ord.txt tam_max"""
import sys
src = open(sys.argv[1]).read()
old = """static inline int finish_real(double r, double lo, double hi, Der *d) {
    if (hi < RLO_SAFE) return D_OUT;
    if (lo > RHI_SAFE) return D_OUT;"""
new = """static long long N_HI_INF = 0, N_HI_INF_OK = 0, N_FR = 0; static double EX_R, EX_LO;
static inline int finish_real_(double r, double lo, double hi, Der *d);
static inline int finish_real(double r, double lo, double hi, Der *d) {
    int res = finish_real_(r, lo, hi, d);
    #pragma omp atomic
    N_FR++;
    if (hi == INFINITY && !(lo > RHI_SAFE)) {
        #pragma omp critical(cnt)
        { N_HI_INF++; if (res == D_OK) { N_HI_INF_OK++; EX_R = r; EX_LO = lo; } }
    }
    return res;
}
static inline int finish_real_(double r, double lo, double hi, Der *d) {
    if (hi < RLO_SAFE) return D_OUT;
    if (lo > RHI_SAFE) return D_OUT;"""
assert old in src
sys.stdout.write(src.replace(old, new))
