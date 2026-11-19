"""Audit of review point 4.4 over every C11 cell: rerun the flagged verification with the instrumented nm.py
(it only counts negative bases whose exponent is not an exact rational but may be an integer)."""
import sys, os, json, math, subprocess, time
H = "/home/user/audit44"
sys.path[:0] = [H + "/src", H + "/c11"]
import certify
from flint import arb
rows = [json.loads(l) for l in open("/home/user/articulo/paquete_v3/pandigital_cert/results/c11.jsonl") if l.strip()]
env = dict(os.environ, PYTHONHASHSEED="0", FULLMAX="3", PANDI_SRC=H + "/src", PYTHONPATH=H + "/src:" + H + "/c11",
           C11_AUDIT_DIR=H + "/out", C11_NOWIDEN="1", OMP_NUM_THREADS="1")
RES = {(9, 0): "/tmp/lab_c7_x/c7/res_9_0.txt"}
def res(n, o):
    return RES.get((n, o), f"/home/user/work/pandigital_cert/src/res_{n}_{o}.txt")
def log(m):
    with open(H + "/audit.log", "a") as f: f.write(time.strftime("%d/%m %H:%M:%S ") + m + "\n")
def verify(n, o, c, E, cands=None, tag=""):
    e = dict(env, WORK=f"{H}/work/{c}_{n}_{o}{tag}")
    if cands: e["C11_CANDS"] = cands
    out = f"{H}/out/res_{c}_{n}_{o}.jsonl"
    r = subprocess.run([sys.executable, "-B", H + "/c11/c11_par.py", str(n), str(o), c, repr(E), res(n, o), "4", out],
                       env=e, capture_output=True, text=True)
    last = [json.loads(l) for l in open(out)][-1] if os.path.exists(out) else None
    log(f"{c} n={n} ord={o} rc={r.returncode} rounds={last and last['rounds'][-1]} unresolved={last and last.get('unresolved')}"
        + ("" if r.returncode == 0 else " ERR " + r.stderr[-500:]))
def chain(n, o, c, E):
    d = f"{H}/work/chain_{c}_{n}_{o}"; os.makedirs(d, exist_ok=True)
    cc = certify.constant(c)
    ew = E * (1 + 1e-6) + 8 * math.ulp(float(cc.mid()))
    a, b = certify.window_doubles(cc, arb(ew))
    open(d + "/w.txt", "w").write(f"{a.hex()} {b.hex()}\n")
    e = dict(env, OMP_NUM_THREADS="4", MOTORA_DEDUP="1")
    subprocess.run([H + "/src/motora", "multi", str(n), str(o), res(n, o), d + "/w.txt", d + "/cand.txt"], env=e, capture_output=True)
    e = dict(env, MOTORC_FAST="1", MOTORC_TPREC="128")
    subprocess.run([H + "/src/motorc", str(n), str(o), d + "/w.txt", d + "/cand.txt", d + "/k1.txt", "1024"], env=e, capture_output=True)
    subprocess.run([H + "/src/motorc", str(n), str(o), d + "/w.txt", d + "/k1.txt", d + "/k2.txt", "2048"], env=e, capture_output=True)
    nk = sum(1 for _ in open(d + "/k2.txt"))
    log(f"{c} n={n} ord={o} chain: {sum(1 for _ in open(d + '/cand.txt'))} candidates, {nk} after motor C")
    verify(n, o, c, E, cands=d + "/k2.txt", tag="_chain")
order = sorted(rows, key=lambda r: (r["n"] == 9 and r["ord"] == 0, r["n"] == 8 and r["ord"] == 0, r["n"]))
log(f"start: {len(order)} cells")
for r in order:
    n, o, c = r["n"], r["ord"], r["const"]
    E = r["rounds"][-1]["E"]
    if os.path.exists(f"{H}/out/res_{c}_{n}_{o}.jsonl"):          # resumable: a cell already finished is skipped
        continue
    try:
        if n == 9 and o == 0:
            verify(n, o, c, E, cands=f"/home/user/recog2/work/c11/{c}/kept2.txt")
        elif n == 8 and o == 0:
            chain(n, o, c, E)
        else:
            verify(n, o, c, E)
    except Exception as ex_:
        log(f"{c} n={n} ord={o} EXCEPTION {ex_!r}")
    na = sum(sum(1 for _ in open(os.path.join(H + "/out", f))) for f in os.listdir(H + "/out") if f.startswith("audit_"))
    log(f"   audit events so far: {na}")
log("DONE")
