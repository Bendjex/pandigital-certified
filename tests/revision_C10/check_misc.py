from mpmath import mp, mpf, log10, floor, e, pi, log, power, mpc, im, re as mre
import sys
mp.dps = 120
# ---- e records: decimals k = floor(log10(N/e))
def k_of(logN10):   # log10 N
    return floor(logN10 - log10(e))
L5 = log10(5); L3 = log10(3); L2 = log10(2)
print("Sabey    N=3^(2^85):  k =", int(k_of(mpf(2) ** 85 * L3)))
print("Bamberger N=5^(3^84): k =", int(k_of(mpf(3) ** 84 * L5)))
print("Fastfaxr N=3^(2^153): k ~", mp.nstr(mpf(2) ** 153 * L3, 6))
k09 = k_of(mpf(3) ** 140 * L5); print("0-9 N=5^(3^140): k =", int(k09), " digits of N =", int(floor(mpf(3) ** 140 * L5)) + 1)
print("1-8 N=2^(2^21): -log10(err) =", mp.nstr(-log10(e / (2 * mpf(2) ** (2 ** 21))) if False else (2 ** 21) * L2 + log10(2) - log10(e), 14), " k =", int(k_of(2 ** 21 * L2)))
for name, N in (("2^72", mpf(2) ** 72), ("5^34", mpf(5) ** 34), ("5^81", mpf(5) ** 81)):
    mp.dps = 250
    err = e - (1 + 1 / N) ** N
    print(name, "error =", mp.nstr(err, 6), " -log10 =", mp.nstr(-log10(err), 12), " k =", int(floor(log10(N / e))), " e/2N =", mp.nstr(e / (2 * N), 6))
mp.dps = 120
print("McKernon ~3.1e30 decimals claim (not checkable here)")
# three-storey example
print("2^(3^49): decimals ~", mp.nstr(mpf(3) ** 49 * L2, 6), "; ln ln N =", mp.nstr(log(mpf(3) ** 49 * log(2)), 8))
for nm, lnN in (("5^(3^84)", mpf(3) ** 84 * log(5)), ("2^(3^84)", mpf(3) ** 84 * log(2)), ("3^(2^85)", mpf(2) ** 85 * log(3)), ("2^(2^70)", mpf(2) ** 70 * log(2)), ("2^72", 72 * log(2)), ("2^(2^21)", mpf(2) ** 21 * log(2)), ("5^(3^140)", mpf(3) ** 140 * log(5))):
    print("  ln ln N of", nm, "=", mp.nstr(log(lnN), 20))
print("SB(3) e^21.456 =", mp.nstr(mp.exp(21.456), 4), "; e^3.032 =", mp.nstr(mp.exp(3.032), 4), "; e^21.5 =", mp.nstr(mp.exp(21.5), 4))
# ---- section 1.1 formulas
c = pi
f1 = mpf(2) ** (mpf(5) ** mpf("0.4")) - mpf("0.6") - (mpf("0.3") ** 9 / 7) ** (mpf("0.8") ** mpf("0.1"))
f2 = 3 + mpf(5) ** (-(mpf(7) ** mpf("0.1"))) + mpf("0.49") ** 8 / mpf(2) ** 6
f3 = (2 / mpf("0.98") - mpf("0.3")) * (mpf("0.4") + mpf(5) ** (mpf(7) ** mpf("-0.6") - mpf("0.1")))
print("Sabey/Ziv pi:", mp.nstr(abs(f1 - c), 5), "| Ioffe:", mp.nstr(abs(f2 - c), 5), "| Vlasii:", mp.nstr(abs(f3 - c), 5))
# ---- section 6 example
print("2^(3^7)/9^(5^4) =", mp.nstr(mpf(2) ** (3 ** 7) / mpf(9) ** (5 ** 4), 5))
# ---- Aizenbud
mp.dps = 250
eps = mpf(4) ** (-72)
val = -mp.sqrt(mpc(-(mpf(8) ** 96))) * (mp.exp(mpc(0, 1) * pi * eps) - mpf(3) ** 0)
print("Aizenbud: re - pi =", mp.nstr(mre(val) - pi, 5), " im =", mp.nstr(im(val), 5))
# ---- delta
print("delta - centre = 0.2016e-30 ->", mp.nstr(mpf("4.669201609102990671853203820466201617258185577") - mpf("4.669201609102990671853203820466"), 5))
# ---- exponential quantiles
print("quantiles:", mp.nstr(-log(mpf("0.95")), 4), mp.nstr(-log(mpf("0.05")), 4))
# n = 1 errors
for nm, cc in (("pi", pi), ("gamma", mp.euler), ("delta", mpf("4.669201609102990671853203820466")), ("e", e)):
    print(" n=1", nm, mp.nstr(abs(1 - cc), 3), end=";")
print()
