"""
C3: certification of one cell (constant c, n digits, mode) - the main theorem for that cell.

  1. W0 = [c - E, c + E] with E a guess of the optimum (the best error found by pisearch5, slightly
     enlarged).  Any E works: if W0 holds no expression, E is doubled and the run repeated.
  2. Motor A lists candidate paths covering every in-range expression with value in W0 (Theorem B).
  3. Motor B verifies each path exactly and returns the values in W0 with their expressions (Theorem C).
  4. The minimum of |v - c| is decided with Arb at increasing precision; values attaining it are
     compared exactly (same value = several writings; different values with equal error = a tie).
  5. The optimal expressions are re-evaluated by an independent evaluator (mpmath, string parsing)
     that also checks the digit rule.

usage: certify.py n ord const E_guess resfile [outjson]
"""
import sys, os, re, json, time, math, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath
from flint import arb, ctx
import motorb, verify, exact as ex
from exact import Unresolved

ctx.prec = 256
DELTA_STR = "4.669201609102990671853203820466"     # Feigenbaum delta, 30 decimals (radius 1e-30)


def constant(name, prec=256):
    old = ctx.prec
    ctx.prec = prec
    try:
        if name == "pi": return arb.pi()
        if name == "e": return arb(1).exp()
        if name == "gamma": return arb.const_euler()
        if name == "delta": return arb(DELTA_STR) + arb(0, 1e-30)
    finally:
        ctx.prec = old
    raise ValueError(name)


def mp_constant(name):
    if name == "pi": return mpmath.pi
    if name == "e": return mpmath.e
    if name == "gamma": return mpmath.euler
    if name == "delta": return mpmath.mpf(DELTA_STR)   # 30 decimals only


def independent_eval(expr, n, ordered, dps=80):
    """independent check of one formula: digit rule, validity and range of EVERY subexpression, value.
    The string is parsed with Python's ast (numbers become mpmath numbers at dps digits) and each node
    is evaluated separately; returns (digits_ok, valid_and_in_range, value, smallest relative margin to
    a range bound)."""
    import ast
    digs = [ch for ch in expr if ch.isdigit()]
    if ordered:
        dok = digs == [str(d) for d in range(1, n + 1)]
    else:
        dok = sorted(digs) == [str(d) for d in range(1, n + 1)]
    py = re.sub(r"(\d*\.\d+|\d+)", lambda m: f"N('{m.group(1)}')", expr).replace("^", "**")
    tree = ast.parse(py, mode="eval").body
    LO, HI = mpmath.mpf(10) ** -300, mpmath.mpf(10) ** 300
    state = {"ok": True, "margin": mpmath.inf}
    def chk(v):
        if v == 0:
            return v
        a = abs(v)
        if not (LO <= a <= HI):
            state["ok"] = False
        else:
            state["margin"] = min(state["margin"], abs(mpmath.log(a / LO)), abs(mpmath.log(HI / a)))
        return v
    def ev(t):
        if isinstance(t, ast.Call):                       # N('digits')
            return chk(mpmath.mpf(t.args[0].value))
        if isinstance(t, ast.UnaryOp) and isinstance(t.op, ast.USub):
            return -ev(t.operand)
        if isinstance(t, ast.BinOp):
            x, y = ev(t.left), ev(t.right)
            if isinstance(t.op, ast.Add): return chk(x + y)
            if isinstance(t.op, ast.Sub): return chk(x - y)
            if isinstance(t.op, ast.Mult): return chk(x * y)
            if isinstance(t.op, ast.Div):
                if y == 0: state["ok"] = False; return mpmath.mpf(0)
                return chk(x / y)
            if isinstance(t.op, ast.Pow):
                if x == 0:
                    if not y > 0: state["ok"] = False
                    return chk(mpmath.mpf(0))
                if x < 0:
                    if y != mpmath.floor(y): state["ok"] = False; return mpmath.mpf(0)
                    return chk(mpmath.power(x, y))
                return chk(mpmath.power(x, y))
        raise ValueError("unexpected node in " + expr)
    with mpmath.workdps(dps):
        v = ev(tree)
    return dok, state["ok"], v, state["margin"]


def window_doubles(c, E):
    lo, hi = c - E, c + E
    a = math.nextafter(float(lo.lower()), -math.inf)
    b = math.nextafter(float(hi.upper()), math.inf)
    return a, b


def certify(n, ordered, cname, E, resfile, workdir="."):
    c = constant(cname)
    # the guess comes from a double-precision search: widen it by more than its rounding error
    E = E * (1 + 1e-6) + 8 * math.ulp(float(c.mid()))
    t0 = time.time()
    V = verify.Verifier(n, ordered, fullmax=int(os.environ.get("FULLMAX", "3")))
    while True:
        a, b = window_doubles(c, arb(E))
        wf = os.path.join(workdir, f"w_{cname}_{n}_{ordered}.txt")
        cf = os.path.join(workdir, f"c_{cname}_{n}_{ordered}.txt")
        with open(wf, "w") as f:
            f.write(f"{a.hex()} {b.hex()}\n")
        if os.environ.get("REUSE_CAND") and os.path.exists(cf) and os.path.getsize(cf) > 0:
            pass                                   # motor A already ran for this exact window file
        else:
            r = subprocess.run([os.environ.get("MOTORA", "./motora"), "multi", str(n), str(ordered), resfile, wf, cf],
                               stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise RuntimeError("motor A failed: " + r.stderr[-2000:])
        tA = time.time() - t0
        if os.environ.get("ONLY_A") or os.path.exists(os.path.join(workdir, "cells", "ONLY_A")):
            print("motor A done:", cf, round(tA, 1), "s")
            sys.exit(0)
        cands = list(dict.fromkeys(l.split(" ", 1)[1].strip() for l in open(cf) if l.strip()))
        lo, hi = arb(a), arb(b)
        vals = {}
        for cand in cands:
            vals.update(V.verify_line(cand, lo, hi))
        # keep values whose error is <= E (the window may be slightly wider)
        errs = []
        for key, (v, e) in vals.items():
            errs.append((v, e))
        if errs:
            break
        E *= 2
    return finish(n, ordered, cname, E, len(cands), len(vals), errs, list(V.B.unresolved), t0, tA)


def finish(n, ordered, cname, E, ncand, nvals, errs, unresolved, t0, tA):
    """minimum over the values in the window (errs: list of (value, writing)), decided with Arb"""
    # decide the minimum: error enclosures at increasing precision (adaptive per value)
    def err(v):
        best = None
        for prec in ex.TPRECS:
            try:
                b = v.ball(prec)
            except Unresolved:
                continue
            old = ctx.prec
            ctx.prec = prec
            try:
                best = abs(b - constant(cname, prec))
            finally:
                ctx.prec = old
            if best.is_finite() and best.rad() <= best.mid() * arb(2) ** -50:
                break
        if best is None:
            raise Unresolved("no enclosure of an error")
        return best
    es = [(err(v), v, e) for (v, e) in errs]
    up = min(es, key=lambda t: float(t[0].upper()))[0].upper()
    # values whose error cannot be shown larger than the smallest upper bound
    close = [t for t in es if not (t[0].lower() > up)]
    mn = min(close, key=lambda t: float(t[0].mid()))
    # group the close values by exact equality; values that agree to 16384 bits but whose equality the
    # exact arithmetic cannot prove are grouped too, with the flag 'proven': False
    groups = []
    for t in close:
        for g in groups:
            d = ex.sub(t[1], g[0][1])
            try:
                same = ex.sign(d) == 0
                proven = True
            except Unresolved:
                b16 = d.ball(16384) if d.kind == "N" else None
                same = b16 is not None and b16.contains(0)
                proven = False
            if same:
                g.append(t)
                if not proven:
                    g[0] = g[0] + ("unproven",) if len(g[0]) == 3 else g[0]
                break
        else:
            groups.append([t])
    opt_groups = []
    for g in groups:
        exprs = sorted({t[2] for t in g}, key=len)
        okv = []
        for e in exprs[:20]:
            dok, rok, val, margin = independent_eval(e, n, ordered)
            with mpmath.workdps(80):
                er = abs(val - mp_constant(cname))
            okv.append({"formula": e, "digits_ok": dok, "valid_in_range": rok,
                        "range_margin_ln": mpmath.nstr(margin, 6), "error_mpmath": mpmath.nstr(er, 25)})
        opt_groups.append({"equality_proven": len(g[0]) == 3,
                           "value": ex.tight_ball(g[0][1]).mid().str(30, radius=False),
                           "error": g[0][0].str(25),
                           "writings": len(exprs), "checked": okv})
    res = {"n": n, "ord": ordered, "const": cname, "E_window": E, "candidates": ncand,
           "values_in_window": nvals, "optimal": opt_groups, "ties": len(groups) > 1,
           "unresolved": len(unresolved), "seconds": round(time.time() - t0, 1), "seconds_A": round(tA, 1)}
    return res


if __name__ == "__main__":
    n, ordered, cname, E, resfile = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
    res = certify(n, ordered, cname, E, resfile)
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 6:
        with open(sys.argv[6], "a") as f:
            f.write(json.dumps(res) + "\n")
