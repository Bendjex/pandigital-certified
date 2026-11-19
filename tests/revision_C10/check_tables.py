"""Reviewer's check of Tables 1-3 and 5 against results/*.jsonl, with an independent parser/evaluator (mpmath)."""
import sys, json, re, os
from mpmath import mp, mpf
mp.dps = 80
MS = sys.argv[1]; RES = sys.argv[2]
DELTA = mpf("4.669201609102990671853203820466201617258185577475768632745651343004134330211314737138689744023948013817165984855189815134408627142027932522312442988890890859944935463236713411532481714219947455644365823793202009561058330575458617652222070385410646749494284981453391726200568755665952339875603825637225")
CONST = {"pi": mp.pi, "e": mp.e, "gamma": mp.euler, "delta": DELTA}

class P:
    """recursive descent: expr := term (('+'|'-') term)*; term := unary (('*'|'/') unary)*;
    unary := '-' unary | power; power := atom ('^' unary)?   (so a^-b^c = a^(-(b^c)), -a^b = -(a^b))"""
    def __init__(self, s):
        self.s = s.replace(" ", "").replace("−", "-").replace("·", "*").replace("×", "*"); self.i = 0; self.digits = []; self.nodes = []
    def peek(self): return self.s[self.i] if self.i < len(self.s) else ""
    def expr(self):
        v = self.term()
        while self.peek() in "+-" and self.peek():
            op = self.s[self.i]; self.i += 1
            w = self.term()
            v = v + w if op == "+" else v - w
            self.nodes.append(v)
        return v
    def term(self):
        v = self.unary()
        while True:
            c = self.peek()
            if c in ("*", "/"):
                self.i += 1
                w = self.unary()
                if c == "/" and w == 0: raise ZeroDivisionError
                v = v * w if c == "*" else v / w
                self.nodes.append(v)
            elif c == "(" or c.isdigit() or c == ".":      # implicit multiplication is not used in the tables
                raise ValueError("implicit multiplication at %d" % self.i)
            else:
                return v
    def unary(self):
        if self.peek() == "-":
            self.i += 1
            v = -self.unary(); self.nodes.append(v); return v
        return self.power()
    def power(self):
        b = self.atom()
        if self.peek() == "^":
            self.i += 1
            e = self.unary()
            if b < 0:
                if e != int(e): raise ValueError("negative base, non-integer exponent")
                v = mpf(b) ** int(e)
            elif b == 0:
                if e <= 0: raise ValueError("0^nonpositive")
                v = mpf(0)
            else:
                v = b ** e
            self.nodes.append(v)
            return v
        return b
    def atom(self):
        if self.peek() == "(":
            self.i += 1
            v = self.expr()
            assert self.peek() == ")", (self.s, self.i)
            self.i += 1
            return v
        m = re.match(r"\d*\.?\d+", self.s[self.i:])
        assert m, (self.s, self.i)
        t = m.group(0); self.i += len(t)
        self.digits += [c for c in t if c.isdigit()]
        v = mpf(t); self.nodes.append(v)
        return v

def evaluate(f):
    p = P(f); v = p.expr()
    assert p.i == len(p.s), (f, p.i)
    return v, p.digits, p.nodes

def check_formula(f, n, order, digits0=False):
    v, dg, nodes = evaluate(f)
    exp = [str(d) for d in (range(0, 10) if digits0 else range(1, n + 1))]
    ok_digits = sorted(dg) == sorted(exp)
    ok_order = (dg == exp) if order else True
    rng = all(x == 0 or mpf(10) ** -300 <= abs(x) <= mpf(10) ** 300 for x in nodes)
    return v, ok_digits, ok_order, rng

# ---------------------------------------------------------------- results files
cert = {}
for fn in ("certified.jsonl", "e_7_1.jsonl"):
    for l in open(os.path.join(RES, fn)):
        r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
c11 = {}
for l in open(os.path.join(RES, "c11.jsonl")):
    r = json.loads(l); c11[(r["const"], r["n"], r["ord"])] = r
print("certified cells:", len(cert), "| c11 cells:", len(c11))
def ball_mid(s):
    m = re.match(r"\[(.*) \+/- (.*)\]", s); return mpf(m.group(1)), mpf(m.group(2))
bad = 0
for k, r in sorted(cert.items()):
    c, n, o = k
    if r["ties"] or r["unresolved"] or len(r["optimal"]) != 1:
        print("CERT anomaly", k, r["ties"], r["unresolved"], len(r["optimal"])); bad += 1
    op = r["optimal"][0]
    mid, rad = ball_mid(op["error"])
    for ch in op["checked"]:
        v, okd, oko, rng = check_formula(ch["formula"], n, o)
        err = abs(v - CONST[c])
        if not (okd and oko and rng): print("CERT formula rule problem", k, ch["formula"], okd, oko, rng); bad += 1
        if abs(err - mid) > rad + mpf(10) ** -30 * 5: print("CERT error mismatch", k, ch["formula"], mp.nstr(err, 30), op["error"]); bad += 1
print("certified.jsonl: formulas re-evaluated, problems:", bad)

# ---------------------------------------------------------------- manuscript tables 2 and 3
ms = open(MS, encoding="utf-8").read()
def table(title):
    i = ms.index(title); j = ms.index("\n\n", ms.index("|---|", i))
    return [ln for ln in ms[i:j].split("\n") if ln.startswith("|")][2:]
def parse_err(s):
    s = s.strip().replace("·10", "e").translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-"))
    return mpf(s)
consts = ["pi", "gamma", "delta", "e"]
bad = 0; ncell = 0
for title, o in (("**Table 2.**", 1), ("**Table 3.**", 0)):
    for row in table(title):
        cells = [x.strip() for x in row.strip("|").split("|")]
        n = int(cells[0])
        for c, cell in zip(consts, cells[1:]):
            if "Section 7" in cell: continue
            m = re.match(r"(.*?) `(.*)`", cell)
            errs, f = m.group(1), m.group(2)
            ncell += 1
            r = cert.get((c, n, o))
            if r is None: print("TABLE cell without certificate", c, n, o); bad += 1; continue
            mid, rad = ball_mid(r["optimal"][0]["error"])
            terr = parse_err(errs)
            # 15 significant digits
            if abs(terr - mid) > abs(mid) * mpf(10) ** -14 * 6: print("TABLE error differs from certificate", c, n, o, errs, r["optimal"][0]["error"]); bad += 1
            if mp.nstr(mid, 15) != mp.nstr(terr, 15): print("TABLE 15-digit rounding", c, n, o, errs, "cert:", mp.nstr(mid, 15)); bad += 1
            v, okd, oko, rng = check_formula(f, n, o)
            err = abs(v - CONST[c])
            if not (okd and oko and rng): print("TABLE formula rule problem", c, n, o, f, okd, oko, rng); bad += 1
            if abs(err - mid) > rad + mpf(10) ** -40: print("TABLE formula does not give the certified error", c, n, o, f, mp.nstr(err, 20)); bad += 1
            fs = [ch["formula"] for ch in r["optimal"][0]["checked"]]
            if f not in fs: print("TABLE formula is not the certificate's writing", c, n, o, f, fs)
print("Tables 2-3:", ncell, "cells with formula checked; problems:", bad)

# ---------------------------------------------------------------- Table 5 (no unary minus)
rows = table("**Table 5.**")
cols = [("pi", 1), ("gamma", 1), ("delta", 1), ("e", 1), ("pi", 0), ("gamma", 0), ("delta", 0), ("e", 0)]
bad = 0; same = 0; tot = 0; same_n7 = []
for row in rows:
    cells = [x.strip() for x in row.strip("|").split("|")]
    n = int(cells[0])
    for (c, o), cell in zip(cols, cells[1:]):
        if cell == "—": 
            if (c, n, o) in c11: print("T5: dash but c11 has the cell", c, n, o); bad += 1
            continue
        tot += 1
        m = re.match(r"([^ ]+)(?: \((.*)\))?$", cell)
        e11 = parse_err(m.group(1)); ewith = parse_err(m.group(2)) if m.group(2) else None
        r = c11.get((c, n, o))
        if r is None: print("T5 cell missing in c11.jsonl", c, n, o); bad += 1; continue
        if r["ties"] or r["unresolved"]: print("T5/c11 anomaly", c, n, o, r["ties"], r["unresolved"]); bad += 1
        mid, rad = ball_mid(r["optimal_no_minus"][0]["error"])
        if mp.nstr(mid, 3) != mp.nstr(e11, 3): print("T5 value differs", c, n, o, cell, mp.nstr(mid, 5)); bad += 1
        rc = cert[(c, n, o)]; midc, _ = ball_mid(rc["optimal"][0]["error"])
        if ewith is None:
            same += 1
            if n >= 7: same_n7.append((c, n, o))
            if abs(mid - midc) > rad + mpf(10) ** -30: print("T5: no bracket but optimum differs", c, n, o, mp.nstr(mid, 8), mp.nstr(midc, 8)); bad += 1
        else:
            if mp.nstr(midc, 3) != mp.nstr(ewith, 3): print("T5 bracket differs from certificate", c, n, o, cell, mp.nstr(midc, 5)); bad += 1
            if not mid > midc: print("T5: bracket but not worse", c, n, o); bad += 1
        # formula
        op = r["optimal_no_minus"][0]
        f = op["formula_no_minus"]
        v, okd, oko, rng = check_formula(f, n, o)
        if "-" in re.sub(r"(?<=[0-9)])-", "", f): print("T5/c11 formula has a unary minus?", c, n, o, f); bad += 1
        if not (okd and oko and rng): print("c11 formula rule problem", c, n, o, f, okd, oko, rng); bad += 1
        if abs(abs(v - CONST[c]) - mid) > rad + mpf(10) ** -40: print("c11 formula does not give its error", c, n, o, f); bad += 1
print("Table 5:", tot, "cells; unchanged optimum:", same, "; unchanged with n >= 7:", same_n7, "; problems:", bad)
print("c11 cells total:", len(c11), "; cells where no_minus == with minus (from files):",
      sum(1 for k, r in c11.items() if abs(ball_mid(r["optimal_no_minus"][0]["error"])[0] - ball_mid(cert[k]["optimal"][0]["error"])[0]) < mpf(10) ** -40))

# ---------------------------------------------------------------- Table 1
bad = 0
names = {"π": "pi", "γ": "gamma", "δ": "delta", "e": "e"}
for row in table("**Table 1.**"):
    cells = [x.strip() for x in row.strip("|").split("|")]
    c = names[cells[0]]; var = cells[1]; n = int(cells[2]); f = cells[3]; errs = cells[4]; prev = cells[5]
    if c == "e" and (n >= 7): 
        print("T1 (e, towers) not evaluated here:", var, n, f, errs); continue
    o = 1 if var == "in order" else 0
    v, okd, oko, rng = check_formula(f, n, o)
    err = abs(v - CONST[c])
    mid, rad = ball_mid(cert[(c, n, o)]["optimal"][0]["error"])
    if not (okd and oko and rng): print("T1 formula rule problem", c, var, n, f, okd, oko, rng); bad += 1
    if abs(err - mid) > rad + mpf(10) ** -40: print("T1 formula is not the certified optimum", c, var, n, f, mp.nstr(err, 12), mp.nstr(mid, 12)); bad += 1
    if mp.nstr(err, 3) != mp.nstr(parse_err(errs), 3): print("T1 error column", c, var, n, errs, mp.nstr(err, 5)); bad += 1
    pm = re.match(r"([0-9.]+·10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", prev)
    if pm and not parse_err(pm.group(1)) > err: print("T1 previous not worse", c, var, n); bad += 1
print("Table 1 (pi, gamma, delta rows): problems:", bad)
