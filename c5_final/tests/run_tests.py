"""Run the tests of motor X, its query rules and the independent evaluator lvl.

usage:  python run_tests.py [quick|n9|all] [--only name,name] [--timeout seconds]

Environment (see REPRODUCIR.md):
  PANDI_SRC      the certifier's src folder (motor B, motor A sources and the resolution files res_<n>_<o>.txt)
  PANDI_RESULTS  the certifier's results folder (certified.jsonl, c11.jsonl)
  QX_NCACHE      (n9 tests) a folder for the persistent NE-table cache of n = 9 (several GB; built on first use)
The code under test is ../codigo/c5 (motor X, qx, ivl, lvl, hitcheck) and ../codigo/cert (motor B helpers); xlib.so
must have been compiled in ../codigo/c5 (REPRODUCIR.md).

A test passes when it exits with code 0 and its summary line reports no problem (PROBLEMS 0, bad 0, missed 0,
UNCOVERED 0, differences 0, ...). Each test's full output goes to logs/<name>.log; the summary to resultados_tests.txt.
"""
import os, re, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
C5 = os.path.join(ROOT, "codigo", "c5")
CERT = os.path.join(ROOT, "codigo", "cert")

# (name, arguments, tier, what it checks: motorx_changes.md item or appendix section)
TESTS = [
    # --- the independent evaluator lvl (Appendix C.9; items 38-50 and patches 4, 5, 7, 7c)
    ("t_lvl", ["2000", "1"], "quick", "lvl against mpmath on random expressions"),
    ("t_lvlanc", ["2000", "1"], "quick", "lvl first-order anchors at exact powers"),
    ("t_lvlpw", ["3000", "1"], "quick", "lvl exact powers through prime-exponent vectors (L1)"),
    ("t_lvlpf", ["600", "1"], "quick", "lvl irrational exact offsets"),
    ("t_lvltw", ["600"], "quick", "lvl canonical towers (L1)"),
    ("t_lvltwc", ["800", "13"], "quick", "lvl powers of powers (L1)"),
    ("t_lvllin", ["300", "17"], "quick", "lvl exact linear forms that cancel (L2)"),
    ("t_lvloff", ["200"], "quick", "lvl X/(X +- c) anchored at 1 (patch 5)"),
    ("t_lvlanq", ["60"], "quick", "lvl powers of huge rationals (patch 5)"),
    ("t_lvlbin", [], "quick", "lvl differences with the true ln inside, exact zeros (patch 5)"),
    ("t_lvlbound", ["200", "5"], "quick", "lvl patches 7/7c: bounds, sign of forms, anchored differences, towers"),
    # --- query rules of motor X (Appendix C.4-C.6, C.9 items 1-37)
    ("t_fwd2", [], "quick", "forward enclosures (Lemma X2)"),
    ("t_pre2", [], "quick", "preimage windows (Lemma X3)"),
    ("t_level2", [], "quick", "second-level constraints"),
    ("t_c2prop", [], "quick", "second-level constraints passed down the recursion"),
    ("t_ancpow", [], "quick", "anchored bases with exact exponents"),
    ("t_ancrat", [], "quick", "rational anchors, exact 1 as 1/1"),
    ("t_offset", [], "quick", "exact offsets of exact powers"),
    ("t_onediff", [], "quick", "item 48: |1 - l| = u, second level of each branch"),
    ("t_unique", [], "quick", "one table entry per distinct statement"),
    ("t_index", [], "quick", "tiered and bounded indices of the tables"),
    ("t_xlibbits", ["@SRC/res_5_0.txt"], "quick", "xlib.so: bit identity of the vector exp/ln bounds (compare hashes)"),
    ("t_hc_smoke", [], "quick", "hitcheck on a small cell (values 828 inside-or-undecided 828 fail 0)"),
    ("t_c11multi", [], "quick", "both windows (C7 and C11) checked in one pass"),
    ("t_q", ["5", "0", "N"], "quick", "recursive queries against the materialized tables (n = 5)"),
    ("t_cover", ["5", "0"], "quick", "coverage of an approximate exhaustive enumeration (n = 5)"),
    ("t_pieces", ["0"], "quick", "a split in pieces gives the same hits as the whole split"),
    ("t_near1", ["5", "12345", "40"], "quick", "exact members of V(S) near 1"),
    ("t_exactv", ["7", "0", "40"], "quick", "exact zero and exact members in recursive V queries (n = 7)"),
    # --- tests that need the tables of n = 9 (res_9_0.txt and the NE-table cache; hours on first use)
    ("t_coinc", [], "n9", "C.6: exact coincidences written through an out-of-range node are found by the query"),
    ("t_exact1", ["3000"], "n9", "item 41: a base exactly 1 against a window with a known sign"),
    ("t_hcleaf", ["0"], "n9", "rule H: exact forms of the members of the tables"),
    ("t_release", ["4"], "n9", "emptying the caches between pieces changes no result"),
    ("t_members", ["10", "10", "50"], "n9", "exact members of every table class (invariant I)"),
    ("t_narrow", [], "n9", "exact narrowing of bases near 1"),
    ("t_vvmemo", ["9", "23456789", "57"], "n9", "repeated (V, V) resolutions answered once"),
    ("t_refine", ["9", "123456789", "1234", "4"], "n9", "refinement of bases near 1"),
    ("t_batchrec", ["9", "234689", "N", "20"], "n9", "item 42: recursive sub-queries in bounded batches"),
    ("t_recbatch", ["8", "234678", "N", "40"], "n9", "recursive sub-queries receive all windows at once"),
    ("t_nearx", ["9", "20", "200000"], "n9", "item 49: exact members of l near each member of s"),
]

# counts of failures in the summary lines; "exact differences N" (t_offset) is a count of cases, not of failures
BAD = re.compile(r"(PROBLEMS|problems|\bbad|BAD|MISSED|missed|UNCOVERED|uncovered|(?<!exact )differences|differ:|differs:|\bfail|FAIL)\s*:?\s*([1-9]\d*)")


def main():
    tier = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "quick"
    only = None
    timeout = 6 * 3600
    for i, a in enumerate(sys.argv):
        if a == "--only":
            only = set(sys.argv[i + 1].split(","))
        if a == "--timeout":
            timeout = int(sys.argv[i + 1])
    src = os.environ.get("PANDI_SRC")
    if not src or not os.path.isdir(src):
        sys.exit("set PANDI_SRC to the certifier's src folder (REPRODUCIR.md)")
    if not os.path.exists(os.path.join(C5, "xlib.dll" if os.name == "nt" else "xlib.so")):
        sys.exit("compile xlib.so in codigo/c5 first (REPRODUCIR.md)")
    env = dict(os.environ)
    env.update({"PYTHONPATH": os.pathsep.join([C5, CERT]), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                "PYTHONDONTWRITEBYTECODE": "1", "PYTHONIOENCODING": "utf-8"})
    env.setdefault("RESDIR", src)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    rows, nfail = [], 0
    for name, args, t, what in TESTS:
        if only is not None and name not in only:
            continue
        if only is None and tier != "all" and t != tier:
            continue
        args = [a.replace("@SRC", src) for a in args]
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, "-B", os.path.join(HERE, name + ".py"), *args], cwd=HERE, env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                               errors="replace", timeout=timeout)
            out, rc = r.stdout, r.returncode
        except subprocess.TimeoutExpired as ex:
            out, rc = (ex.stdout or "") if isinstance(ex.stdout, str) else "", "timeout"
        dt = time.time() - t0
        with open(os.path.join(HERE, "logs", name + ".log"), "w", encoding="utf-8") as f:
            f.write(out)
        lines = [l for l in out.strip().splitlines() if l.strip()]
        last = lines[-1] if lines else ""
        bad = [m.group(0) for m in BAD.finditer("\n".join(lines[-5:]))]
        ok = rc == 0 and not bad
        if name == "t_cover":
            # approximate coverage: measure.py's float values not inside an enclosure.  The first version of motor X
            # leaves the same 237 (n = 5), all with |ln u| >= 699 or measure.py's +-1e18 sentinel (float ghosts of
            # rounding residues, Appendix C.7); any other count is reported
            miss = sum(int(m) for m in re.findall(r"MISSING (\d+)", out))
            ok = rc == 0 and miss == 237
            last = f"MISSING {miss} in total (first version of motor X: 237)" + ("" if ok else " REVISAR")
        nfail += not ok
        row = f"{'OK   ' if ok else 'FALLO'} {name:12s} {dt:7.0f} s  rc {rc}  | {last[:220]}"
        print(row, flush=True)
        rows.append(row)
    with open(os.path.join(HERE, "resultados_tests.txt"), "a", encoding="utf-8") as f:
        f.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} tier {tier}{' only ' + ','.join(sorted(only)) if only else ''}\n")
        f.write("\n".join(rows) + "\n")
    print(f"{len(rows) - nfail}/{len(rows)} OK")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
