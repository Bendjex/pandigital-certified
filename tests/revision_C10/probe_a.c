#define NO_MAIN
#include "motora.c"
int main(void) {
    RLO_SAFE = pred(1e-300); RHI_SAFE = succ(1e300);
    /* code_for with a non-finite upper bound */
    double r = 1e300, lo = 9.99e299, hi = INFINITY;
    int c = code_for(r, lo, hi);
    double slo, shi; hull_of(r, c, &slo, &shi);
    printf("code_for(1e300, 9.99e299, inf) = %d  hull = [%a, %a]  covers hi? %d\n", c, slo, shi, shi >= hi);
    Der d; int res = finish_real(r, lo, hi, &d);
    printf("finish_real(r=1e300, lo=9.99e299, hi=inf) -> %d (D_OK=%d, D_SING_W=%d, D_OUT=%d)  code=%d\n", res, D_OK, D_SING_W, D_OUT, res == D_OK ? d.code : -1);
    r = 1e299; lo = 1e299; hi = INFINITY;
    res = finish_real(r, lo, hi, &d);
    printf("finish_real(r=1e299, lo=1e299, hi=inf) -> %d  code=%d\n", res, res == D_OK ? d.code : -1);
    /* DBL_MAX upper bound (finite) for comparison */
    res = finish_real(1e299, 1e299, DBL_MAX, &d);
    printf("finish_real(r=1e299, lo=1e299, hi=DBL_MAX) -> %d\n", res);
    printf("L2E_LO=%a L2E_HI=%a  pred(1e-300)=%a succ(1e300)=%a  1e300=%.25g 1e-300=%.25g\n", L2E_LO, L2E_HI, RLO_SAFE, RHI_SAFE, 1e300, 1e-300);
    return 0;
}
