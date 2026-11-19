/* roundeven() (C23) is missing from the MinGW runtime; this exact replacement is used only in the
   Windows build.  Round to nearest integer, ties to even, exact for all doubles. */
#include <math.h>
double roundeven(double x) {
    double t = trunc(x);
    if (fabs(x - t) == 0.5) return 2.0 * round(x / 2.0);   /* tie: x/2 is exact, round half away of x/2 */
    return round(x);
}
