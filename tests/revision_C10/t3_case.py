import sys, os
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np, mpmath
from mpmath import mp, iv
from fractions import Fraction as Fr
import ivl, motorx, lvl
sys.set_int_max_str_digits(0)
mp.prec = 2000
# [l^s | [l^-s | [|s-l| | [l^-s | [s^l | 56 | 28] | [l^-s | 95 | 8]] | 6] | 86] | 1/10]
# innermost: s^l with s=56, l=28: 56^28.   [l^-s | 95 | 8] = 8^-95.   [l^-s | X=56^28 | Y=8^-95] = Y^-X = (8^-95)^-(56^28) = 8^(95*56^28)
# [|s-l| | Z | 6] = |Z - 6|.   [l^-s | that | 86] = 86^-(|Z-6|).  [l^s | W | 1/10] = (1/10)^W,  W = 86^-(Z-6) tiny => value = 1 - W ln 10 ...
X = mp.mpf(56) ** 28
lnZ = 95 * X * mp.log(8)                 # ln Z, Z = 8^(95 * 56^28)  (huge)
# W = 86^-(Z - 6):  ln W = -(Z - 6) ln 86 ~ -Z ln 86 ;  ln|ln W| = ln Z + ln ln 86 + ln(1 - 6/Z)
lnlnW = lnZ + mp.log(mp.log(86))         # (6/Z is negligible: Z ~ e^(1e51))
# u = (1/10)^W = exp(-W ln 10): ln u = -W ln 10;  ln|ln u| = ln W + ln ln 10 = -e^(lnlnW) + ln ln 10 = -e^nu
# nu = ln(e^lnlnW - ln ln 10) = lnlnW + ln(1 - lnln10 * e^-lnlnW) = lnlnW (to all digits)
nu = lnlnW
print("nu =", mp.nstr(nu, 40))
print("ivl M3, N3 = 1.7576608412311433e+51, 1.75766084123342e+51")
print("M3 <= nu <= N3 :", mp.mpf(1.7576608412311433e+51) <= nu <= mp.mpf(1.75766084123342e+51))
print("rel width of [M3,N3]:", (1.75766084123342e+51 - 1.7576608412311433e+51) / 1.7576608412311433e+51)
