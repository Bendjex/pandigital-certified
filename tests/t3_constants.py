"""Check of the two constants trusted by motor A (Lemma 7 of Appendix A) and of the digits of delta (Section 2.2).
lambda1 = 0x1.71547652b82fep+0 and lambda2 = succ(lambda1) must enclose log2(e); the 31 digits of delta used by the
certifier must agree with OEIS A006890 (its first 35 digits are copied below)."""
import math
import mpmath as mp
mp.mp.dps = 60
l1 = float.fromhex("0x1.71547652b82fep+0"); l2 = math.nextafter(l1, math.inf)
ok1 = mp.mpf(l1) < 1 / mp.log(2) < mp.mpf(l2)
DELTA_CERT = "4.669201609102990671853203820466"
# first 35 significant digits of OEIS A006890 (checked against oeis.org on 2 Nov 2026; Section 2.2 of the article)
A006890 = "4.6692016091029906718532038204662016"
ok2 = A006890.startswith(DELTA_CERT) and abs(mp.mpf(A006890) - mp.mpf(DELTA_CERT)) < mp.mpf("1e-30")
print("lambda1 < log2 e < lambda2:", ok1, "| delta digits agree with A006890 and lie in the ball of radius 1e-30:", ok2)
raise SystemExit(0 if ok1 and ok2 else 1)
