"""
Motor B (C3 certification) with the verification of the candidate paths split across processes.
Same result as certify.py on the same candidate file:
  - each worker runs verify.Verifier.verify_line on chunks of candidate lines;
  - the parent applies the results IN THE ORIGINAL ORDER of the candidates (ordered imap), exactly like
    certify.py's `vals.update(...)`, so the stored writing of each value is the same;
  - the workers send the exact values themselves (pickled without their Arb-ball cache);
  - the minimum is decided by certify.finish (the code of certify.py).
Resumable: every finished chunk is appended to <out>.chunks.jsonl; a restart skips them.
usage: certify_par.py n ord const E_guess candfile workers out.jsonl [windowfile line]
  E_guess is widened as in certify.py; if a window file and line are given, the window must coincide.
"""
import sys, os, json, time, math, pickle, base64
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify, verify, exprparse, exact as ex
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
    G["V"] = verify.Verifier(n, o, fullmax=fullmax)
    G["lo"], G["hi"] = arb(a), arb(b)


def work(job):
    ci, chunk = job
    V = G["V"]
    u0 = len(V.B.unresolved)
    out = []
    for cand in chunk:
        out.append([(k, v, e) for k, (v, e) in V.verify_line(cand, G["lo"], G["hi"]).items()])
    return ci, out, len(V.B.unresolved) - u0


def main():
    n, o, cname, Eg, cf, workers, outf = (int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]),
                                          sys.argv[5], int(sys.argv[6]), sys.argv[7])
    fullmax = int(os.environ.get("FULLMAX", "3"))
    chunk = int(os.environ.get("CHUNK", "40"))
    c = certify.constant(cname)
    E = Eg * (1 + 1e-6) + 8 * math.ulp(float(c.mid()))
    a, b = certify.window_doubles(c, arb(E))
    if len(sys.argv) > 9:
        wl = open(sys.argv[8]).read().split("\n")[int(sys.argv[9])].split()
        if (a.hex(), b.hex()) != (wl[0], wl[1]) and (float.fromhex(wl[0]), float.fromhex(wl[1])) != (a, b):
            raise SystemExit(f"window mismatch: computed {a.hex()} {b.hex()}, file {wl}")
    t0 = time.time()
    cands = list(dict.fromkeys(l.split(" ", 1)[1].strip() for l in open(cf) if l.strip()))
    chunks = [(i, cands[i * chunk:(i + 1) * chunk]) for i in range((len(cands) + chunk - 1) // chunk)]
    ck = outf + f".{cname}.chunks.jsonl"
    done = {}
    if os.path.exists(ck):
        for l in open(ck):
            try:
                r = json.loads(l)
            except ValueError:
                break                                   # a truncated last line (interrupted write)
            done[r["ci"]] = r
    print(f"{cname} {n} {o}: {len(cands)} candidate paths, {len(chunks)} chunks, {len(done)} already done", flush=True)
    vals, keys = {}, {}
    unres = 0
    nd, tB = 0, time.time()
    todo = [j for j in chunks if j[0] not in done]

    def apply(ci, out, du):
        nonlocal unres, nd
        unres += du
        for items in out:
            for k, v, e in items:
                if k not in vals:
                    if v.key() != k:
                        raise AssertionError("value and key from the worker disagree: " + e)
                    vals[k] = (v, e)
                else:
                    vals[k] = (vals[k][0], e)           # same value; the later writing, as vals.update
        nd += len(chunks[ci][1])

    with open(ck, "a") as fck:
        pool = Pool(workers, initializer=init, initargs=(n, o, a, b, fullmax)) if todo else None
        try:
            res_iter = pool.imap(work, todo) if pool else iter(())
            # apply in chunk order: chunks already done come from the checkpoint file
            for ci in range(len(chunks)):
                if ci in done:
                    r = done[ci]
                    apply(ci, pickle.loads(base64.b64decode(r["out"])), r["du"])
                else:
                    cj, out, du = next(res_iter)
                    assert cj == ci
                    fck.write(json.dumps({"ci": ci, "out": base64.b64encode(pickle.dumps(out)).decode(), "du": du}) + "\n")
                    fck.flush()
                    apply(ci, out, du)
                if ci % max(1, len(chunks) // 100) == 0 or ci == len(chunks) - 1:
                    el = time.time() - tB
                    print(f"  [progress] {cname}: {nd}/{len(cands)} paths ({100 * nd / max(1, len(cands)):.1f}%), "
                          f"{el:.0f}s, values {len(vals)}", flush=True)
        finally:
            if pool:
                pool.close(); pool.join()
    errs = list(vals.values())
    if not errs:
        raise SystemExit("no value in the window: unexpected for a window around a known formula")
    res = certify.finish(n, o, cname, E, len(cands), len(vals), errs, [None] * unres, t0, 0.0)
    res["workers"] = workers; res["fullmax"] = fullmax
    print(json.dumps(res, indent=1), flush=True)
    with open(outf, "a") as f:
        f.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
