"""
C11: certification of one cell in the variant without unary minus.

Same pipeline as certify.py (C3), with the flagged motor B (nm.py):
  1. window W = [c - E, c + E];  motor A lists candidate paths covering every C1 expression in range
     with |value| in W (a superset of the minus-free ones);
  2. the flagged verifier returns every magnitude in W with its POS/NEG flags;
  3. only magnitudes with POS (reachable as +v without unary minus) count; if there is none, E is
     multiplied by 4 and the run is repeated (the with-minus optimum is a lower bound);
  4. minimum decided with Arb; minus-free formula re-evaluated independently (digits, validity,
     range, and absence of unary minus).
usage: c11_certify.py n ord const E_guess resfile [outjson]
"""
import sys, os, json, time, math, subprocess, ast
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SRC = os.environ.get("PANDI_SRC", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, SRC)
import nm, certify, exact as ex
from exact import Unresolved
from flint import arb, ctx
import mpmath


def has_unary_minus(expr):
    py = expr.replace("^", "**")
    import re
    py = re.sub(r"(\d*\.\d+|\d+)", lambda m: "N", py)
    tree = ast.parse(py, mode="eval")
    return any(isinstance(t, ast.UnaryOp) for t in ast.walk(tree))


def run(n, ordered, cname, E, resfile, workdir):
    c = certify.constant(cname)
    E = E * (1 + 1e-6) + 8 * math.ulp(float(c.mid()))
    t0 = time.time()
    V = nm.NVerifier(n, ordered, fullmax=int(os.environ.get("FULLMAX", "3")))
    rounds = []
    while True:
        a, b = certify.window_doubles(c, arb(E))
        wf = os.path.join(workdir, f"w_{cname}_{n}_{ordered}_{E:.3e}.txt")
        cf = os.path.join(workdir, f"c_{cname}_{n}_{ordered}_{E:.3e}.txt")
        with open(wf, "w") as f:
            f.write(f"{a.hex()} {b.hex()}\n")
        if not (os.path.exists(cf) and os.path.getsize(cf) > 0):
            r = subprocess.run([os.path.join(SRC, "motora"), "multi", str(n), str(ordered), resfile, wf, cf],
                               stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise RuntimeError("motor A failed: " + r.stderr[-2000:])
        cands = list(dict.fromkeys(l.split(" ", 1)[1].strip() for l in open(cf) if l.strip()))
        lo, hi = arb(a), arb(b)
        vals = {}
        for cand in cands:
            for k, ent in V.verify_line(cand, lo, hi).items():
                V.B.merge(vals, k, ent[0], ent[1], ent[2], ent[3], ent[4])
        pos = [ent for ent in vals.values() if ent[2] & nm.POS]
        rounds.append({"E": E, "candidates": len(cands), "values": len(vals), "pos_values": len(pos)})
        print("round", rounds[-1], flush=True)
        if pos:
            break
        E *= 4
    return finish(n, ordered, cname, pos, rounds, t0, len(V.B.unresolved), len(nm.NB.WARN))


def finish(n, ordered, cname, pos, rounds, t0, n_unres, n_warn):
    """minimum among the POS magnitudes, decided with Arb; independent re-evaluation of the formula"""
    def err(v):
        best = None
        for prec in ex.TPRECS:
            try:
                bb = v.ball(prec)
            except Unresolved:
                continue
            old = ctx.prec; ctx.prec = prec
            try:
                best = abs(bb - certify.constant(cname, prec))
            finally:
                ctx.prec = old
            if best.is_finite() and best.rad() <= best.mid() * arb(2) ** -50:
                break
        if best is None:
            raise Unresolved("no enclosure of an error")
        return best
    es = [(err(ent[0]), ent) for ent in pos]
    up = min(es, key=lambda t: float(t[0].upper()))[0].upper()
    close = [t for t in es if not (t[0].lower() > up)]
    groups = []
    for t in close:
        for g in groups:
            d = ex.sub(t[1][0], g[0][1][0])
            try:
                same = ex.sign(d) == 0
            except Unresolved:
                same = False
            if same:
                g.append(t); break
        else:
            groups.append([t])
    opt = []
    for g in groups:
        e_nm = g[0][1][3]
        dok, rok, val, margin = certify.independent_eval(e_nm, n, ordered)
        with mpmath.workdps(80):
            er = abs(val - certify.mp_constant(cname))
        opt.append({"value": ex.tight_ball(g[0][1][0]).mid().str(30, radius=False),
                    "error": g[0][0].str(25), "formula_no_minus": e_nm,
                    "has_unary_minus": has_unary_minus(e_nm), "digits_ok": dok, "valid_in_range": rok,
                    "error_mpmath": mpmath.nstr(er, 25), "sign_ok": bool(val > 0),
                    "with_minus_writing": g[0][1][1]})
    return {"n": n, "ord": ordered, "const": cname, "rounds": rounds, "optimal_no_minus": opt,
            "ties": len(groups) > 1, "unresolved": n_unres, "warn_parity": n_warn,
            "seconds": round(time.time() - t0, 1)}


if __name__ == "__main__":
    n, o, cname, E, resfile = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
    wd = os.environ.get("WORK", "run")
    os.makedirs(wd, exist_ok=True)
    res = run(n, o, cname, E, resfile, wd)
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 6:
        with open(sys.argv[6], "a") as f:
            f.write(json.dumps(res) + "\n")
