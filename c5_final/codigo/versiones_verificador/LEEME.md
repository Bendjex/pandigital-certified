# Versions of the checker (lvl.py + hitcheck.py) used in the records of n = 9 in any order

Fingerprint = first 16 hex digits of SHA-256(lvl.py + hitcheck.py), the two files concatenated in that order, with
CRLF turned into LF (`python huella.py <folder>`; the same as `_code()` in `c5run.py`, written in the field `code` of
every record since 13 November 2026).

| Folder | Fingerprint | Date | Records used in the closure (pieces / splits) |
|---|---|---|---|
| p0_original_c90dbdc006cb5be9 | c90dbdc006cb5be9 | 11-12 Nov | together with p3: 1303 / 40 (splits 0-16 of the laboratory; no `code` field) |
| p3_f70deff8e1e04e0e | f70deff8e1e04e0e | 13 Nov | (see above: p3 changes only printing; the records cannot tell them apart) |
| p4_5e2fbc16882dfa91 | 5e2fbc16882dfa91 | 13-14 Nov | 95 / 11 |
| p5_72a64c06e8a95b0f | 72a64c06e8a95b0f | 14-17 Nov | 1086 / 111 |
| p7_bc9d449da8d28463 | bc9d449da8d28463 | 17-18 Nov | 3459 / 597 |
| p7c_1e05db8f91fe15ec | 1e05db8f91fe15ec | 18 Nov (final) | 6 / 6 |

Patch 6 (fingerprint d2ead1fbf0e0e33b, sign error 8) is not here: no record of it is used.
hitcheck.py is the same file (025afa1d...) in p4, p5, p7 and p7c.

To repeat a record with its own checker, copy lvl.py and hitcheck.py of its folder into codigo/c5 (REPRODUCIR.md).

## Difference between p0 (original) and p3: printing only

The repr of a huge exact rational, and of a value that cannot be printed, no longer raises an exception. In p0 an
exception would have stopped the piece before any record was written, so every record of p0 is one p3 would have
written, and vice versa.

```diff
--- p0_original_c90dbdc006cb5be9/lvl.py	2026-11-18 21:52:41.388352230 +0000
+++ p3_f70deff8e1e04e0e/lvl.py	2026-11-18 21:52:41.392106050 +0000
@@ -75,7 +75,12 @@
         if self.far:
             return f"<far {self.far}>"
         if self.q is not None:
-            return f"<q {self.q if self.q.denominator < 10**20 else float(self.q)}>"
+            if self.q.denominator < 10**20:
+                return f"<q {self.q}>"
+            try:
+                return f"<q {float(self.q)}>"
+            except OverflowError:            # (13/11) beyond the doubles: only the order of magnitude, from the bit lengths
+                return f"<q ~2^{self.q.numerator.bit_length() - self.q.denominator.bit_length()}>"
         if self.anc is not None:
             aq, sg, D = self.anc
             return f"<{aq} (1 {'+' if sg > 0 else '-'} e^[{mpmath.nstr(D.a, 8)}, {mpmath.nstr(D.b, 8)}])>"
--- p0_original_c90dbdc006cb5be9/hitcheck.py	2026-11-18 21:52:41.389643547 +0000
+++ p3_f70deff8e1e04e0e/hitcheck.py	2026-11-18 21:52:41.393313885 +0000
@@ -152,24 +152,32 @@
     res, undecided = [], 0
     for i in range(len(H)):
         vals = C.hit(int(H.pid[i]))
-        bad = [(e, repr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
+        bad = [(e, srepr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
         prec = lvl.PREC
         for p in PRECS[1:]:                  # undecided: the same evaluation with more bits (lvl.set_prec)
             if not bad:
                 break
             lvl.set_prec(p); C.memo.clear(); prec = p
             vals = C.hit(int(H.pid[i]))
-            bad = [(e, repr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
+            bad = [(e, srepr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
         if lvl.PREC != PRECS[0]:
             lvl.set_prec(PRECS[0]); C.memo.clear()
         undecided += len(bad)
         res.append({"const": c, "enclosure_ln": [float(H.A[i]), float(H.B[i])], "values": len(vals), "prec_bits": prec,
-                    "undecided": bad[:5], "sample": [[e[:400], repr(v)] for (v, e) in vals[:3]]})
+                    "undecided": bad[:5], "sample": [[e[:400], srepr(v)] for (v, e) in vals[:3]]})
         print(f"[progress] {c} hit {i + 1}/{len(H)} ({100 * (i + 1) / len(H):.1f}%) values {len(vals)} "
               f"undecided {len(bad)}", (bad[0][0][:300] if bad else ""), flush=True)
     return undecided, res
 
 
+def srepr(v):
+    """(13/11) repr for the records: a value that cannot be printed must never stop the check"""
+    try:
+        return repr(v)
+    except Exception as ex:
+        return f"<unprintable value: {type(ex).__name__}>"
+
+
 def certainly_inside(v, lo, hi):
     """True if the magnitude v is certainly inside [lo, hi] (an exact rational, or an enclosure, within it): such a value
     stays inside at any precision; used only to stop raising the precision and to report it"""
@@ -206,14 +214,14 @@
         if lvl.PREC != PRECS[0]:
             lvl.set_prec(PRECS[0]); C.memo.clear()
         rec = {"const": c, "enclosure_ln": [float(H.A[i]), float(H.B[i])], "values": len(vals), "prec_bits": prec,
-               "sample": [[e[:400], repr(v)] for (v, e) in vals[:3]], "undecided": []}
+               "sample": [[e[:400], srepr(v)] for (v, e) in vals[:3]], "undecided": []}
         for k, b in bad.items():
             ins = sum(1 for (e, v) in b if certainly_inside(v, *vws[k]))
             tot[k][0] += len(b); tot[k][1] += ins
-            rec["undecided_" + k] = [(e, repr(v)) for (e, v) in b[:5]]
+            rec["undecided_" + k] = [(e, srepr(v)) for (e, v) in b[:5]]
             rec["inside_" + k] = ins
             if b:
-                rec["undecided"] = rec["undecided"] or [(e, repr(v)) for (e, v) in b[:5]]
+                rec["undecided"] = rec["undecided"] or [(e, srepr(v)) for (e, v) in b[:5]]
         res.append(rec)
         print(f"[progress] {c} hit {i + 1}/{len(H)} ({100 * (i + 1) / len(H):.1f}%) values {len(vals)} " +
               " ".join(f"undecided_{k} {len(b)}" for k, b in bad.items()), flush=True)
```
