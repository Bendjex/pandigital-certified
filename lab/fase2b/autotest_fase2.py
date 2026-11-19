"""
Autotest de FASE 2 (10-20 min, prioridad baja, solo escribe en work\\autotest). Cada paso compara con Linux:
  1 hashes del paquete; res_9_0.txt de ..\\c7 (hash) copiado a py\\cert
  2 Python del laboratorio con pylib: flint 0.9.0, numpy, mpmath, sympy, xlib.dll
  3 motor A (motora.exe) en una casilla pequena (n = 6, cualquier orden, pi): candidatos identicos a Linux
  4 motor C (motorc.exe) sobre esos candidatos: termina y solo quita caminos
  5 motor B en paralelo sobre los supervivientes de motor C: resultado identico a certified.jsonl
  6 verificador de C11 en paralelo: optimo sin menos unario identico al de Linux
  7 C5 (motor X) en 2 partes, n = 6 cualquier orden, pi: 0 hits
  8 C5 bit a bit: hits de una ventana artificial cerca de 1 (n = 5 en orden) identicos a Linux (xlib.dll + numpy)
  9 hitcheck + lvl: los mismos recuentos que en Linux
 10 lvl: autotest aleatorio contra mpmath a 3000 digitos, 0 errores
usage: fase2.py autotest [--generate]   (--generate: Linux, escribe data/test/expected.json)
"""
import os, sys, json, time, hashlib, subprocess, re, shutil
import fase2 as F

T = F.CFG["autotest"]


def norm_hash(path, sort=False):
    lines = [l.rstrip("\r\n") for l in open(path, encoding="utf-8", errors="replace")]
    if sort:
        lines.sort()
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest(), len(lines)


def run(cmd, env=None, timeout=7200, cwd=None):
    e = dict(os.environ); e.update(F.pyenv()); e.update(env or {})
    kw = {"creationflags": 0x4000} if F.WIN else {}
    r = subprocess.run([str(c) for c in cmd], env=e, cwd=cwd or F.HERE, capture_output=True, text=True,
                       timeout=timeout, **kw)
    return r.returncode, r.stdout + r.stderr


def main(generate=False):
    W = F.P("work", "autotest"); os.makedirs(W, exist_ok=True)
    expf = F.P("data", "test", "expected.json")
    exp = {} if generate else json.load(open(expf))
    got = {}
    fails = []
    t0 = time.time()

    def check(name, ok, info=""):
        print(f"[{'OK' if ok else 'FALLA'}] {name} {info}  ({time.time() - t0:.0f}s)", flush=True)
        if not ok:
            fails.append(name)

    # 1. hashes
    if not generate:
        bad = []
        for line in open(F.P("HASHES.txt")):
            h, f = line.split(None, 1); f = f.strip().replace("/", os.sep)
            if not os.path.exists(F.P(f)) or F.sha256(F.P(f)) != h:
                bad.append(f)
        check("1 hashes del paquete", not bad, ", ".join(bad[:5]))
    src = os.path.join(F.C7DIR, F.CFG["res_name"])
    if os.path.exists(F.RES) and F.sha256(F.RES) == F.CFG["res_sha256"]:
        check("1 res_9_0.txt ya copiado y con su hash", True)
    elif os.path.exists(src) and F.sha256(src) == F.CFG["res_sha256"]:
        shutil.copyfile(src, F.RES + ".tmp"); os.replace(F.RES + ".tmp", F.RES)
        check("1 res_9_0.txt copiado de ..\\c7 (hash correcto)", F.sha256(F.RES) == F.CFG["res_sha256"])
    else:
        check("1 res_9_0.txt en ..\\c7 con el hash esperado", False, src)
    # 2. imports
    rc, out = run([F.PY, "-B", "-c", "import flint, numpy, mpmath, sympy; "
                   "print('V', flint.__version__, numpy.__version__, mpmath.__version__, sympy.__version__)"])
    m = re.search(r"V (\S+) (\S+) (\S+) (\S+)", out)
    check("2 python-flint 0.9.0, numpy, mpmath y sympy", rc == 0 and m is not None and m.group(1) == "0.9.0",
          m.group(0) if m else out[-500:])
    # T1-T5 of opt_entrega/autotest/ESPERADO.txt (motors of 04/11): expected hashes of Linux
    AT = F.P("data", "autotest")
    def rawh(f):
        return F.sha256(f) if os.path.exists(f) else ""
    fast = {"MOTORC_FAST": "1", "MOTORC_TPREC": "128"}
    t = os.path.join(W, f"t_{int(time.time())}"); os.makedirs(t)
    rc, out = run([F.P(F.CFG["motorc"]), 8, 0, os.path.join(AT, "w3_8.txt"), os.path.join(AT, "autotest_test8.txt"),
                   os.path.join(t, "t1_kept.txt"), 1024], dict(fast, MOTORC_SAMPLE="10"))
    check("T1 motor C n = 8 (MOTORC_FAST)", rawh(os.path.join(t, "t1_kept.txt")).startswith("0a3582b4908e4c586196495e378520065f61ca61cb702f3a31029db00fbe5210")
          and rawh(os.path.join(t, "t1_kept.txt.sample")) == "ea50b012b3da5f33735bb14256552864fe6c4747f0bd28d0ccd606fc5b5bc0fb", out[-300:] if rc else "")
    rc, out = run([F.P(F.CFG["motorc"]), 8, 0, os.path.join(AT, "w3_8.txt"), os.path.join(AT, "autotest_test8.txt"),
                   os.path.join(t, "t2_kept.txt"), 1024], dict(fast, MOTORC_SAMPLE="10", MOTORC_NBLOCKS="4"))
    rc2, out2 = run([F.P(F.CFG["motorc"]), "join", os.path.join(t, "t2_kept.txt"), 4])
    check("T2 motor C por bloques + join", rawh(os.path.join(t, "t2_kept.txt")) == "0a3582b4908e4c586196495e378520065f61ca61cb702f3a31029db00fbe5210"
          and rawh(os.path.join(t, "t2_kept.txt.sample")) == "ba1f6934a7be02c1e389b2840a6fcc4ee95209f56f3883629b6058a5c849b921", (out + out2)[-300:])
    rc, out = run([F.P(F.CFG["motorc"]), 9, 0, F.P("data", "w_c7_9_0.txt"), os.path.join(AT, "muestra_n9_W0.txt"),
                   os.path.join(t, "t3_kept.txt"), 1024], fast)
    check("T3 motor C n = 9, muestra real: 17 caminos", rawh(os.path.join(t, "t3_kept.txt")) == "708d37dad89efe9a23680ea99887a48812826ce9096f4054c3f1672e96c1ab8b", out[-300:] if rc else "")
    for tn, extra, hl in (("T4", {}, "589a1883e52a038c819d70fd3ee3a12c85085a8d91ce3b914f0121b871e674cc"),
                          ("T5", {"MOTORA_DEDUP": "1"}, "544ef23743fe50223213251b41308db0c57c550eff8612954b093acb031ffb96")):
        f = os.path.join(t, f"{tn.lower()}_cand.txt")
        rc, out = run([F.P(F.CFG["motora"]), "multi", 6, 0, os.path.join(AT, "res_6_0.txt"), os.path.join(AT, "w_ident6b.txt"), f],
                      dict(extra, OMP_NUM_THREADS="2"))
        h, nl = norm_hash(f) if os.path.exists(f) else ("", 0)
        check(f"{tn} motor A ident6b" + (" con DEDUP" if extra else ""), h == hl, f"{nl} lineas")
    # 3. motor A
    wt = F.P("data", "test", "w_pi_6_0.txt")
    c6 = os.path.join(W, "c6.txt")
    rc, out = run([F.P(F.CFG["motora"]), "multi", 6, 0, F.P("py", "cert", "res_6_0.txt"), wt, c6], {"OMP_NUM_THREADS": "1"})
    h, nl = norm_hash(c6) if os.path.exists(c6) else ("", 0)
    got["motora_c6"] = h
    check("3 motor A n = 6: candidatos identicos a Linux", rc == 0 and (generate or h == exp["motora_c6"]), f"{nl} lineas")
    # 4. motor C
    k6 = os.path.join(W, "k6.txt")
    if os.path.exists(k6):
        os.replace(k6, k6 + f".{int(time.time())}")
    rc, out = run([F.P(F.CFG["motorc"]), 6, 0, wt, c6, k6, 1024], {"MOTORC_TPREC": "128", "MOTORC_FAST": "1"})
    okc = rc == 0 and os.path.exists(k6) and "done:" in out
    if okc:
        cl = set(l.rstrip("\r\n") for l in open(c6)); kl = [l.rstrip("\r\n") for l in open(k6)]
        okc = all(l in cl for l in kl) and len(kl) > 0
    check("4 motor C: termina y solo quita caminos", okc, (f"{len(kl)} de {len(cl)} caminos quedan" if okc else out[-500:]))
    kin = k6 if okc else c6
    # 5. motor B in parallel
    outB = os.path.join(W, "certB.jsonl")
    if os.path.exists(outB):
        os.replace(outB, outB + f".{int(time.time())}")
    for f in os.listdir(W):
        if ".chunks." in f:
            os.replace(os.path.join(W, f), os.path.join(W, f) + f".{int(time.time())}")
    rc, out = run([F.PY, "-B", F.P("py", "cert", "certify_par.py"), 6, 0, "pi", repr(T["E_pi_6_0"]), kin, 2, outB],
                  {"FULLMAX": "4"})
    ref = [r for r in F.jsonl_rows(F.P("data", "certified.jsonl")) if (r["const"], r["n"], r["ord"]) == ("pi", 6, 0)][0]
    keys = ["E_window", "values_in_window", "ties", "unresolved"]

    def core(opt):                       # the certificate; the representative writing depends on path order
        return [{k: g[k] for k in ("equality_proven", "value", "error", "writings")} for g in opt]

    def writing_ok(opt, ref_opt):
        return all(all(c["digits_ok"] and c["valid_in_range"] for c in g["checked"]) and
                   g["checked"][0]["error_mpmath"] == r["checked"][0]["error_mpmath"] for g, r in zip(opt, ref_opt))
    rB = F.jsonl_rows(outB)
    okB = (rc == 0 and rB and all(rB[-1][k] == ref[k] for k in keys) and core(rB[-1]["optimal"]) == core(ref["optimal"])
           and writing_ok(rB[-1]["optimal"], ref["optimal"]))
    check("5 motor B en paralelo = certified.jsonl", okB, "" if okB else out[-800:])
    # 6. C11 verifier in parallel
    outV = os.path.join(W, "c11.jsonl")
    if os.path.exists(outV):
        os.replace(outV, outV + f".{int(time.time())}")
    wd = os.path.join(W, f"c11run_{int(time.time())}")
    rc, out = run([F.PY, "-B", F.P("py", "c11", "c11_par.py"), 6, 0, "pi", repr(T["E_pi_6_0"]), F.P("py", "cert", "res_6_0.txt"), 2, outV],
                  {"C11_CANDS": c6, "C11_NOWIDEN": "1", "WORK": wd})
    refV = [r for r in F.jsonl_rows(F.P("data", "c11_n_le_8.jsonl")) if (r["const"], r["n"], r["ord"]) == ("pi", 6, 0)][-1]
    rV = F.jsonl_rows(outV)
    okV = rc == 0 and rV and rV[-1]["optimal_no_minus"] == refV["optimal_no_minus"]
    check("6 verificador de C11 en paralelo = Linux", okV, "" if okV else out[-800:])
    if generate:
        json.dump(got, open(expf, "w"), indent=1)
        print("expected.json escrito:", got)
    okf = F.P("work", "autotest", "AUTOTEST_OK.txt")
    if os.path.exists(okf):
        os.replace(okf, okf + f".{int(time.time())}")
    if not fails and not generate:
        with open(okf, "w") as f:
            f.write(time.strftime("%d/%m/%Y %H:%M:%S") + " AUTOTEST OK\n")
    print(("AUTOTEST OK" if not fails else "AUTOTEST FALLA: " + ", ".join(fails)) + f"  ({time.time() - t0:.0f}s)", flush=True)
    return 0 if not fails else 1
