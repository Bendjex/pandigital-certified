"""
Autotest de FASE 2c (unos 10-20 min, prioridad baja, solo escribe en work\\autotest). Cada paso se compara con Linux:
  1  hashes del paquete; res_9_0.txt de ..\\c7 (hash) copiado a py\\cert
  2  Python del laboratorio con pylib: python-flint 0.9.0, mpmath, sympy
  T1-T5  motores de 04/11 (motor C y motor A) con las salidas de Linux
  3  identificacion de procesos (PID + hora de creacion), parar por PID y memoria libre de Windows
  4  la cadena REAL del orquestador en pequeno (n = 6, pi): motor A -> motor C 1024 por bloques -> join -> motor C 2048
     -> join -> verificador, con el vigilante de memoria forzado una vez (para y reanuda con un proceso menos);
     optimo sin menos unario identico al de Linux
  5  verificador de C11 con n = 9 sobre caminos reales (FULLMAX y reciclado de procesos de la configuracion):
     mismos recuentos que en Linux, y la velocidad en este ordenador
usage: fase2.py autotest      (Linux: python3 autotest_fase2.py --generate escribe data/test/expected.json)
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
    ts = time.strftime("%Y%m%d_%H%M%S")
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
    rc, out = run([F.PY, "-B", "-c", "import flint, mpmath, sympy; "
                   "print('V', flint.__version__, mpmath.__version__, sympy.__version__)"])
    m = re.search(r"V (\S+) (\S+) (\S+)", out)
    check("2 python-flint 0.9.0, mpmath y sympy", rc == 0 and m is not None and m.group(1) == "0.9.0",
          m.group(0) if m else out[-500:])
    # T1-T5 of opt_entrega/autotest/ESPERADO.txt (motors of 04/11): expected hashes of Linux
    AT = F.P("data", "autotest")

    def rawh(f):
        return F.sha256(f) if os.path.exists(f) else ""
    fast = {"MOTORC_FAST": "1", "MOTORC_TPREC": "128"}
    t = os.path.join(W, f"t_{ts}"); os.makedirs(t)
    rc, out = run([F.P(F.CFG["motorc"]), 8, 0, os.path.join(AT, "w3_8.txt"), os.path.join(AT, "autotest_test8.txt"),
                   os.path.join(t, "t1_kept.txt"), 1024], dict(fast, MOTORC_SAMPLE="10"))
    check("T1 motor C n = 8 (MOTORC_FAST)", rawh(os.path.join(t, "t1_kept.txt")) == "0a3582b4908e4c586196495e378520065f61ca61cb702f3a31029db00fbe5210"
          and rawh(os.path.join(t, "t1_kept.txt.sample")) == "ea50b012b3da5f33735bb14256552864fe6c4747f0bd28d0ccd606fc5b5bc0fb", out[-300:] if rc else "")
    rc, out = run([F.P(F.CFG["motorc"]), 8, 0, os.path.join(AT, "w3_8.txt"), os.path.join(AT, "autotest_test8.txt"),
                   os.path.join(t, "t2_kept.txt"), 1024], dict(fast, MOTORC_SAMPLE="10", MOTORC_NBLOCKS="4"))
    rc2, out2 = run([F.P(F.CFG["motorc"]), "join", os.path.join(t, "t2_kept.txt"), 4])
    check("T2 motor C por bloques + join", rawh(os.path.join(t, "t2_kept.txt")) == "0a3582b4908e4c586196495e378520065f61ca61cb702f3a31029db00fbe5210"
          and rawh(os.path.join(t, "t2_kept.txt.sample")) == "ba1f6934a7be02c1e389b2840a6fcc4ee95209f56f3883629b6058a5c849b921", (out + out2)[-300:])
    rc, out = run([F.P(F.CFG["motorc"]), 9, 0, F.P(F.CFG["w_c7"]), os.path.join(AT, "muestra_n9_W0.txt"),
                   os.path.join(t, "t3_kept.txt"), 1024], fast)
    check("T3 motor C n = 9, muestra real: 17 caminos", rawh(os.path.join(t, "t3_kept.txt")) == "708d37dad89efe9a23680ea99887a48812826ce9096f4054c3f1672e96c1ab8b", out[-300:] if rc else "")
    for tn, extra, hl in (("T4", {}, "589a1883e52a038c819d70fd3ee3a12c85085a8d91ce3b914f0121b871e674cc"),
                          ("T5", {"MOTORA_DEDUP": "1"}, "544ef23743fe50223213251b41308db0c57c550eff8612954b093acb031ffb96")):
        f = os.path.join(t, f"{tn.lower()}_cand.txt")
        rc, out = run([F.P(F.CFG["motora"]), "multi", 6, 0, os.path.join(AT, "res_6_0.txt"), os.path.join(AT, "w_ident6b.txt"), f],
                      dict(extra, OMP_NUM_THREADS="2"))
        h, nl = norm_hash(f) if os.path.exists(f) else ("", 0)
        check(f"{tn} motor A ident6b" + (" con DEDUP" if extra else ""), h == hl, f"{nl} lineas")
    # 3. process identity, stop by PID, free memory
    me = F.ptoken(os.getpid())
    p = F.launch([F.PY, "-B", "-c", "import time; time.sleep(600)"], os.path.join(t, "dormir.log"))
    tok = F.ptoken(p.pid)
    alive1 = F.ours(p.pid, tok)
    other = F.ours(p.pid, str(int(tok) + 1) if tok else "1")      # same pid, another creation time: not ours
    F.kill_tree(p.pid)
    try:
        p.wait(timeout=30)
    except Exception:
        pass
    alive2 = F.ours(p.pid, tok)
    ph, pc = F.mem_gb()
    check("3 procesos por PID + hora de creacion, parar por PID, memoria libre",
          me is not None and tok is not None and alive1 and not other and not alive2 and ph > 0 and pc > 0,
          f"(RAM libre {ph:.1f} GB, virtual libre {pc:.1f} GB)")
    # 4. the real chain of the orchestrator, n = 6, pi, any order
    base = os.path.join(W, f"cadena_{ts}")
    X = {"n": 6, "o": 0, "consts": ["pi"], "E": {"pi": T["E_pi_6_0"]}, "res": F.P("py", "cert", "res_6_0.txt"),
         "wlines": [l.strip() for l in open(F.P("data", "test", "w_pi_6_0.txt")) if l.strip()],
         "base": base, "logs": os.path.join(base, "logs"), "out": os.path.join(base, "c11_6_0.jsonl"),
         "nb1": 4, "nb2": 2, "parts": 2, "cores": 2, "workers": 2, "fullmax": F.CFG["fullmax"], "maxtasks": 2,
         "need": {k: 0.2 for k in F.CFG["need_gb"]}, "hours": F.CFG["hours"]}
    old_log = F.LOG
    F.LOG = os.path.join(base, "cadena.log")
    os.makedirs(X["logs"], exist_ok=True)
    steps = F.build_steps(X)
    S, procs = F.new_state(), {}
    guard = None
    try:
        tend = time.time() + 3600
        while time.time() < tend:
            fin = F.run_once(S, X, steps, procs)
            F.save_state(S, os.path.join(base, "estado.json"))
            st = S["steps"].get(steps[-1].name, {})
            if guard is None and st.get("status") == "en curso":
                os.environ["FASE2_TEST_MEM"] = "0.5,0.5"          # memory "runs out" while the verifier runs
                F.run_once(S, X, steps, procs)
                del os.environ["FASE2_TEST_MEM"]
                guard = (S["steps"][steps[-1].name]["status"], S["workers"].get("pi"))
            if fin:
                break
            time.sleep(0.5)
    finally:
        os.environ.pop("FASE2_TEST_MEM", None)
        F.LOG = old_log
    sts = {s.name: S["steps"].get(s.name, {}).get("status") for s in steps}
    cand = os.path.join(base, "pi", "cand.txt")
    k2 = os.path.join(base, "pi", "kept2.txt")
    h6, n6 = norm_hash(cand, sort=True) if os.path.exists(cand) else ("", 0)
    got["c6_dedup_sorted"] = h6
    subset = os.path.exists(k2) and set(l.rstrip("\r\n") for l in open(k2)) <= set(l.rstrip("\r\n") for l in open(cand))
    nk2 = sum(1 for _ in open(k2)) if os.path.exists(k2) else -1
    refV = [r for r in F.jsonl_rows(F.P("data", "c11_n_le_8.jsonl")) if (r["const"], r["n"], r["ord"]) == ("pi", 6, 0)][-1]
    rV = F.jsonl_rows(X["out"])
    okV = bool(rV) and rV[-1]["optimal_no_minus"] == refV["optimal_no_minus"]
    check("4a cadena n = 6: los 6 pasos terminados", all(v == "hecho" for v in sts.values()), str(sts))
    check("4b motor A con DEDUP: candidatos identicos a Linux", generate or h6 == exp["c6_dedup_sorted"], f"{n6} lineas")
    check("4c motor C: solo quita caminos", subset and nk2 > 0, f"{nk2} de {n6} quedan")
    check("4d vigilante de memoria: para el verificador y lo relanza con 1 proceso", guard == ("pendiente", 1), str(guard))
    check("4e verificador: optimo sin menos unario = Linux", okV,
          (rV[-1]["optimal_no_minus"][0]["formula_no_minus"] if okV else str(rV[-1:])[:500]))
    # 5. the verifier with n = 9 on real paths (C11 window of pi)
    s9 = F.P("data", "test", "muestra_c11_n9.txt")      # 80 real paths of C7 (pi), two of them with values
    out9 = os.path.join(W, f"n9_{ts}.jsonl")
    t9 = time.time()
    rc, out = run([F.PY, "-B", F.P("py", "c11", "c11_par.py"), 9, 0, "pi", repr(F.CFG["E_c11"]["pi"]), F.RES, 2, out9],
                  {"C11_CANDS": s9, "C11_NOWIDEN": "1", "WORK": os.path.join(W, f"n9run_{ts}"),
                   "FULLMAX": str(F.CFG["fullmax"]), "C11_MAXTASKS": "1"})
    t9 = time.time() - t9
    r9 = F.jsonl_rows(out9)
    key9 = None
    if r9:
        rd = r9[-1]["rounds"][-1]
        key9 = {"candidates": rd["candidates"], "values": rd["values"], "pos_values": rd["pos_values"],
                "optimal": [(o["error"], o["formula_no_minus"]) for o in r9[-1].get("optimal_no_minus", [])],
                "unresolved": r9[-1]["unresolved"]}
    got["n9"] = key9
    check("5 verificador n = 9 (caminos reales) = Linux", rc == 0 and key9 is not None and (generate or key9 == exp["n9"]),
          f"80 caminos en {t9:.0f} s con 2 procesos ({2 * t9 / 80:.2f} s por camino y proceso)"
          + ("" if rc == 0 and key9 is not None and (generate or key9 == exp["n9"]) else " " + str(key9) + out[-500:]))
    if generate:
        json.dump(got, open(expf, "w"), indent=1)
        print("expected.json escrito:", got)
    okf = F.P("work", "autotest", "AUTOTEST_OK.txt")
    if os.path.exists(okf):
        os.replace(okf, okf + f".{ts}")
    if not fails and not generate:
        with open(okf, "w") as f:
            f.write(time.strftime("%d/%m/%Y %H:%M:%S") + " AUTOTEST OK\n")
    print(("AUTOTEST OK" if not fails else "AUTOTEST FALLA: " + ", ".join(fails)) + f"  ({time.time() - t0:.0f}s)", flush=True)
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main("--generate" in sys.argv))
