"""
C11 for the large cells: the same pipeline as c11_certify.py (motor A candidates, flagged motor B,
minimum decided with Arb), with the verification of the candidate paths split across processes.

Each worker runs its own NVerifier over chunks of candidate lines and returns, per magnitude, its key,
motor B's writing, the POS/NEG flags and the minus-free writings. Exact values are sent
without their cache of Arb balls (which is what made them unpicklable).

The first window is E_start (normally an upper bound given by a verified minus-free formula, so that
one pass suffices); without a POS magnitude the window is multiplied by 4, as in c11_certify.py.
usage: c11_par.py n ord const E_start resfile workers [outjson]
env: C11_CANDS (given candidate paths), C11_NOWIDEN, WORK, FULLMAX, C11_MAXTASKS (optional worker recycling)
"""
import sys, os, json, time, math, subprocess, pickle, base64
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SRC = os.environ.get("PANDI_SRC", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, SRC)
import nm, certify, exact as ex, exprparse
import c11_certify as cc
from flint import arb


def _x_getstate(self):                 # exact values travel between processes without their Arb-ball cache
    return (self.kind, self.terms, self.tree, self._key)


def _x_setstate(self, st):
    self.kind, self.terms, self.tree, self._key = st
    self._ball = {}


ex.X.__getstate__ = _x_getstate
ex.X.__setstate__ = _x_setstate

G = {}


def init(n, o, a, b, fullmax):
    G["V"] = nm.NVerifier(n, o, fullmax=fullmax)
    G["lo"], G["hi"] = arb(a), arb(b)


def work(chunk):
    V = G["V"]
    vals = {}
    u0, w0 = len(V.B.unresolved), len(nm.NB.WARN)
    t0 = time.time()
    for cand in chunk:
        for k, ent in V.verify_line(cand, G["lo"], G["hi"]).items():
            V.B.merge(vals, k, ent[0], ent[1], ent[2], ent[3], ent[4])
    out = [(k, ent[0], ent[1], ent[2], ent[3], ent[4]) for k, ent in vals.items()]
    return out, len(V.B.unresolved) - u0, len(nm.NB.WARN) - w0, len(chunk), time.time() - t0


def run(n, ordered, cname, E, resfile, workdir, workers, chunk=40):
    c = certify.constant(cname)
    E = E * (1 + 1e-6) + 8 * math.ulp(float(c.mid()))
    t0 = time.time()
    fullmax = int(os.environ.get("FULLMAX", "3"))
    B = nm.NB(ordered=bool(ordered), fullmax=fullmax)          # only for merge()
    rounds = []
    n_unres = n_warn = 0
    while True:
        a, b = certify.window_doubles(c, arb(E))
        wf = os.path.join(workdir, f"w_{cname}_{n}_{ordered}_{E:.3e}.txt")
        cf = os.path.join(workdir, f"c_{cname}_{n}_{ordered}_{E:.3e}.txt")
        with open(wf, "w") as f:
            f.write(f"{a.hex()} {b.hex()}\n")
        tA = time.time()
        if os.environ.get("C11_CANDS"):
            cf = os.environ["C11_CANDS"]                 # given candidate paths (e.g. motor C survivors)
        elif not (os.path.exists(cf) and os.path.getsize(cf) > 0):
            r = subprocess.run([os.path.join(SRC, "motora"), "multi", str(n), str(ordered), resfile, wf, cf],
                               stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise RuntimeError("motor A failed: " + r.stderr[-2000:])
        tA = time.time() - tA
        cands = list(dict.fromkeys(l.split(" ", 1)[1].strip() for l in open(cf) if l.strip()))
        print(f"window E={E:.4e}: {len(cands)} candidate paths (motor A {tA:.0f}s)", flush=True)
        chunks = [cands[i:i + chunk] for i in range(0, len(cands), chunk)]
        ck = os.path.join(workdir, f"chunks_{cname}_{n}_{ordered}_{E:.3e}.jsonl")   # checkpoint
        done = {}
        if os.path.exists(ck):
            for l in open(ck):
                try:
                    r = json.loads(l)
                except ValueError:
                    break
                done[r["ci"]] = r
        vals = {}
        nd, tB = 0, time.time()
        todo = [ci for ci in range(len(chunks)) if ci not in done]

        def apply(out):
            for k, v0, e, fl, epos, eneg in out:
                if k not in vals:
                    v = v0
                    if v.key() != k:
                        raise AssertionError("magnitude and key from the worker disagree: " + e)
                else:
                    v = vals[k][0]
                B.merge(vals, k, v, e, fl, epos, eneg)

        mt = int(os.environ.get("C11_MAXTASKS", "0")) or None   # recycle a worker after mt chunks (bounds its memory)
        pool = Pool(workers, initializer=init, initargs=(n, ordered, a, b, fullmax), maxtasksperchild=mt) if todo else None
        try:
            it = pool.imap(work, [chunks[ci] for ci in todo]) if pool else iter(())
            with open(ck, "a") as fck:
                for ci in range(len(chunks)):            # applied in candidate order, as c11_certify.py
                    if ci in done:
                        r = done[ci]; out = pickle.loads(base64.b64decode(r["out"])); du, dw = r["du"], r["dw"]
                    else:
                        out, du, dw, m, sec = next(it)
                        fck.write(json.dumps({"ci": ci, "out": base64.b64encode(pickle.dumps(out)).decode(),
                                              "du": du, "dw": dw}) + "\n"); fck.flush()
                    n_unres += du; n_warn += dw; nd += len(chunks[ci])
                    apply(out)
                    if ci % max(1, len(chunks) // 100) == 0 or ci == len(chunks) - 1:
                        el = time.time() - tB
                        print(f"  [progress] {cname}: {nd}/{len(cands)} paths ({100 * nd / max(1, len(cands)):.1f}%), "
                              f"{el:.0f}s, magnitudes {len(vals)}", flush=True)
        finally:
            if pool:
                pool.close(); pool.join()
        pos = [ent for ent in vals.values() if ent[2] & nm.POS]
        rounds.append({"E": E, "candidates": len(cands), "values": len(vals), "pos_values": len(pos),
                       "seconds_A": round(tA, 1), "seconds_B": round(time.time() - tB, 1), "workers": workers})
        print("round", rounds[-1], flush=True)
        if pos:
            break
        if os.environ.get("C11_NOWIDEN"):
            return {"n": n, "ord": ordered, "const": cname, "rounds": rounds, "optimal_no_minus": [],
                    "status": "no POS magnitude in this window", "unresolved": n_unres, "warn_parity": n_warn,
                    "seconds": round(time.time() - t0, 1)}
        E *= 4
    return cc.finish(n, ordered, cname, pos, rounds, t0, n_unres, n_warn)


if __name__ == "__main__":
    n, o, cname, E, resfile = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
    workers = int(sys.argv[6])
    wd = os.environ.get("WORK", f"run_{n}_{o}")
    os.makedirs(wd, exist_ok=True)
    res = run(n, o, cname, E, resfile, wd, workers)
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 7:
        with open(sys.argv[7], "a") as f:
            f.write(json.dumps(res) + "\n")
