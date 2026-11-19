"""
FASE 2 (laboratorio): orquestador de lo que queda de n = 9 en cualquier orden.
  C7: espera a motor A (..\\c7) -> motor C (1024 y 2048 bits) -> motor B en paralelo (certify_par)
      -> verificador de C11 sobre los mismos supervivientes (atajo: si el optimo se escribe sin menos unario)
  C5: motor X en partes independientes (cells_qpar), reanudable por separaciones; hitcheck si hay hits
  C11: motor A con la ventana de C11 (una constante por run) -> motor C -> verificador de C11 (c11_par)
Reglas: solo escribe dentro de esta carpeta; nunca borra (renombra); prioridad baja; ante algo inesperado
para la tarea afectada y lo escribe en fase2.log; el resto sigue. Todo es reanudable: relanzar continua.
usage: fase2.py run | estado | autotest [--generate]
"""
import sys, os, json, time, subprocess, shutil, hashlib, re, zipfile, glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WIN = os.name == "nt"
CFG = json.load(open(os.path.join(HERE, "fase2_config.json")))
P = lambda *a: os.path.normpath(os.path.join(HERE, *a))
WORK = P("work")
LOG = P("fase2.log")
STATE = P("fase2_state.json")
PY = sys.executable
IDLE, BELOW = (0x40, 0x4000)


# ---------------------------------------------------------------- system helpers (Windows and Linux)
RUNNING_RAM = [0.0]                      # estimated RAM of the jobs in progress (for the test hook)


def avail_gb():
    if os.environ.get("FASE2_TEST_AVAIL_GB"):    # tests only: simulated machine, free = total - estimates
        return float(os.environ["FASE2_TEST_AVAIL_GB"]) - RUNNING_RAM[0]
    if WIN:
        import ctypes
        class M(ctypes.Structure):
            _fields_ = [("l", ctypes.c_ulong), ("load", ctypes.c_ulong), ("tot", ctypes.c_ulonglong),
                        ("av", ctypes.c_ulonglong), ("tp", ctypes.c_ulonglong), ("ap", ctypes.c_ulonglong),
                        ("tv", ctypes.c_ulonglong), ("avv", ctypes.c_ulonglong), ("ave", ctypes.c_ulonglong)]
        m = M(); m.l = ctypes.sizeof(M)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.av / 2 ** 30
    for l in open("/proc/meminfo"):
        if l.startswith("MemAvailable"):
            return int(l.split()[1]) / 2 ** 20
    return 0.0


def alive(pid):
    if not pid:
        return False
    if WIN:
        import ctypes
        k = ctypes.windll.kernel32
        h = k.OpenProcess(0x1000, False, int(pid))
        if not h:
            return False
        code = ctypes.c_ulong()
        k.GetExitCodeProcess(h, ctypes.byref(code)); k.CloseHandle(h)
        return code.value == 259
    try:
        os.kill(int(pid), 0)
    except OSError:
        return False
    try:                                            # zombie children count as finished
        st = open(f"/proc/{int(pid)}/stat").read().split(")")[-1].split()[0]
        return st != "Z"
    except OSError:
        return False


def set_prio(pid, cls):
    try:
        if WIN:
            import ctypes
            k = ctypes.windll.kernel32
            h = k.OpenProcess(0x0200, False, int(pid))
            if h:
                k.SetPriorityClass(h, cls); k.CloseHandle(h)
        else:
            os.setpriority(os.PRIO_PROCESS, int(pid), 19 if cls == IDLE else 10)
    except Exception as e:
        log(f"aviso: no se pudo cambiar la prioridad de {pid}: {e}")


def kill_tree(pid):
    if WIN:
        subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"], capture_output=True)
    else:
        try:
            os.kill(int(pid), 15)
        except OSError:
            pass


def launch(cmd, logf, env=None, cwd=None, prio=BELOW):
    os.makedirs(os.path.dirname(logf), exist_ok=True)
    lf = open(logf, "a")
    lf.write(f"\n##### {time.strftime('%d/%m/%Y %H:%M:%S')} {' '.join(map(str, cmd))}\n"); lf.flush()
    e = dict(os.environ); e.update(env or {})
    kw = {}
    if WIN:
        kw["creationflags"] = 0x08000000 | prio
    else:
        kw["preexec_fn"] = (lambda: os.nice(19 if prio == IDLE else 10))
    p = subprocess.Popen([str(c) for c in cmd], stdout=lf, stderr=subprocess.STDOUT, env=e, cwd=cwd or HERE, **kw)
    return p


def log(msg):
    line = f"{time.strftime('%d/%m/%Y %H:%M:%S')}  {msg}"
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")
    print(line, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tail(path, nbytes=200000):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2); sz = f.tell(); f.seek(max(0, sz - nbytes))
            return f.read().decode("utf-8", "replace")
    except OSError:
        return ""


def pyenv():
    sep = os.pathsep
    return {"PYTHONPATH": sep.join([P("pylib"), P("py", "c5"), P("py", "cert"), P("py", "c11")]),
            "PYTHONHASHSEED": "0", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONIOENCODING": "utf-8",
            "PANDI_SRC": P("py", "cert"), "PANDI_RESULTS": P("data"), "OMP_NUM_THREADS": "1"}


# ---------------------------------------------------------------- the jobs
N, O = CFG["n"], CFG["ord"]
CONSTS = CFG["consts"]
C7DIR = P(CFG["c7_dir"])
RES = P("py", "cert", f"res_{N}_{O}.txt")
RAM = CFG["ram"]
LABEL = {"pi": "π", "gamma": "γ", "delta": "δ"}


def c7_progress():
    """state of the external C7 run: (% , done?)"""
    t = tail(os.path.join(C7DIR, "c7.log"), 400000)
    per = []
    for k in range(len(CONSTS)):
        if re.search(rf"^W{k} candidates:", t, re.M):
            per.append(100.0)
        else:
            m = re.findall(rf"^W{k} progress \d+/\d+ \(([\d.]+)%\)", t, re.M)
            per.append(float(m[-1]) if m else 0.0)
    done = os.path.exists(os.path.join(C7DIR, "c7_fin.txt")) and all(p == 100.0 for p in per)
    return sum(per) / len(per), done, per


class Job:
    resumable = False
    weight = 1.0
    ram = 1.0

    def __init__(self, name, deps=()):
        self.name, self.deps = name, list(deps)

    def logf(self):
        return P("work", "logs", self.name + ".log")

    def ready(self, S):
        return True

    def build(self, S):                  # -> (cmd, env, prio) or None for an internal job
        raise NotImplementedError

    def internal(self, S):
        raise NotImplementedError

    def ok(self, S):
        raise NotImplementedError

    def progress(self, S):
        return 0.0

    def outputs_partial(self):
        return []


def last_pct(text, pat):
    m = re.findall(pat, text)
    return float(m[-1]) if m else 0.0


class C5Part(Job):
    resumable = True

    def __init__(self, k):
        super().__init__(f"c5_p{k}")
        self.k, self.J = k, CFG["c5_parts"]
        self.out = P("work", "c5", f"qpar_{k}.jsonl")
        self.weight = CFG["w"]["c5"] / self.J
        self.ram = RAM["c5"]

    def build(self, S):
        env = pyenv(); env.update({f"E_{c}": repr(CFG["E_c7"][c]) for c in CONSTS})
        prio = IDLE                      # filler: only the CPU the rest leaves free; pre-empted for RAM
        return ([PY, "-B", P("py", "c5", "cells_qpar.py"), N, O, CFG["c5_tv"], CFG["c5_tn"], ",".join(CONSTS),
                 self.k, self.J, self.out], env, prio)

    def summaries(self):
        out = {}
        if os.path.exists(self.out):
            for l in open(self.out):
                try:
                    r = json.loads(l)
                except ValueError:
                    continue
                out[r["const"]] = r
        return out

    def ok(self, S):
        s = self.summaries()
        return all(c in s and s[c]["splits_done"] == s[c]["splits"] for c in CONSTS)

    def progress(self, S):
        if self.ok(S):
            return 100.0
        s = self.summaries()
        t = tail(self.logf(), 20000)
        cur = [c for c in CONSTS if not (c in s and s[c]["splits_done"] == s[c]["splits"])]
        done = len(CONSTS) - len(cur)
        p = 0.0
        if cur:
            m = re.findall(rf"\[progress\] {cur[0]} part {self.k}/{self.J}: split (\d+)/(\d+)", t)
            if m:
                a, b = map(int, m[-1]); p = a / max(1, b)
        return 100.0 * (done + p) / len(CONSTS)

    def detail(self, S):
        s = self.summaries()
        hits = sum(r["hits"] for r in s.values())
        return f"parte {self.k}/{self.J}, {len([c for c in CONSTS if c in s and s[c]['splits_done'] == s[c]['splits']])}/{len(CONSTS)} constantes, hits {hits}"


class Ext(Job):
    """external process job with a fixed command"""
    def __init__(self, name, deps, mk, okf, progf, weight, ram, resumable=False, ready=None, detail=None, partial=()):
        super().__init__(name, deps)
        self.mk, self.okf, self.progf, self.weight, self.ram, self.resumable = mk, okf, progf, weight, ram, resumable
        self._ready, self._detail, self.partial = ready, detail, list(partial)

    def outputs_partial(self):
        return self.partial

    def ready(self, S):
        return self._ready(S) if self._ready else True

    def build(self, S):
        return self.mk(S)

    def ok(self, S):
        return self.okf(S)

    def progress(self, S):
        return 100.0 if self.ok(S) else self.progf(S)

    def detail(self, S):
        return self._detail(S) if self._detail else ""


class Internal(Job):
    def __init__(self, name, deps, fn, okf, weight=0.0):
        super().__init__(name, deps)
        self.fn, self.okf, self.weight, self.ram = fn, okf, weight, 0.0

    def internal(self, S):
        self.fn(S)

    def ok(self, S):
        return self.okf(S)

    def progress(self, S):
        return 100.0 if self.ok(S) else 0.0


def jsonl_rows(path):
    out = []
    if os.path.exists(path):
        for l in open(path, encoding="utf-8"):
            try:
                out.append(json.loads(l))
            except ValueError:
                pass
    return out


def motorc_jobs(tag, deps, wfile, cand_src_fn, outdir, consts_fn):
    """motor C of 04/11 (opt_entrega, HANDOFF_OPT 5-6): block mode, resumable. Pass 1 at 1024 bits (J parts over
    NB1 blocks) -> join -> pass 2 at 2048 bits (J parts over NB2 blocks) -> join -> split by window"""
    J = CFG["motorc_parts"]
    jobs = []

    def stage(bits, nb, src_fn, kept, depn):
        parts = []
        for k in range(J):
            mine = [b for b in range(nb) if b % J == k]

            def mk(S, k=k, nb=nb, src_fn=src_fn, kept=kept, bits=bits):
                env = {"MOTORC_FAST": "1", "MOTORC_TPREC": "128", "MOTORC_SAMPLE": "1000", "MOTORC_PROGRESS": "600",
                       "MOTORC_NBLOCKS": str(nb), "MOTORC_PART": f"{k}/{J}"}
                return ([P(CFG["motorc"]), N, O, wfile, src_fn(S), kept, bits], env, BELOW)

            def okp(S, mine=mine, kept=kept):
                return all(os.path.exists(f"{kept}.b{b:04d}.done") for b in mine)

            name = f"{tag}_{bits}_{k}"

            def progp(S, mine=mine, kept=kept, name=name):
                d = sum(os.path.exists(f"{kept}.b{b:04d}.done") for b in mine)
                t = tail(P("work", "logs", name + ".log"), 20000)
                m = re.findall(r"\[\d+\] block \d+/\d+: ([\d.]+)% read", t)
                frac = float(m[-1]) / 100 if (m and d < len(mine)) else 0.0
                return 100.0 * (d + frac) / max(1, len(mine))
            parts.append(Ext(name, depn, mk, okp, progp, CFG["w"]["motorc"] * (1.0 if bits == 1024 else 0.1) / J,
                             RAM["motorc"], resumable=True))
        jn = f"{tag}_{bits}_join"

        def mkj(S, kept=kept, nb=nb):
            return ([P(CFG["motorc"]), "join", kept, nb], {}, BELOW)

        def okj(S, kept=kept, jn=jn):
            return os.path.exists(kept) and re.search(r"^joined \d+ blocks", tail(P("work", "logs", jn + ".log"), 5000), re.M) is not None
        jj = Ext(jn, [p.name for p in parts], mkj, okj, lambda S: 0.0, 0.0, 0.5, resumable=True)
        return parts + [jj]

    k1 = os.path.join(outdir, "kept1.txt")
    k2 = os.path.join(outdir, "kept2.txt")
    s1 = stage(1024, CFG["motorc_nblocks"], cand_src_fn, k1, deps)
    s2 = stage(2048, CFG["motorc_nblocks2"], lambda S: k1, k2, [s1[-1].name])
    marker = os.path.join(outdir, "split_done.txt")

    def split(S):
        cs = consts_fn(S)
        fs = {k: open(os.path.join(outdir, f"cand_{c}.txt.tmp"), "w") for k, c in enumerate(cs)}
        n = {c: 0 for c in cs}
        for l in open(k2):
            m = re.match(r"W(\d+) ", l)
            if m and int(m.group(1)) in fs:
                fs[int(m.group(1))].write(l.rstrip("\r\n") + "\n"); n[cs[int(m.group(1))]] += 1
        for k, c in enumerate(cs):
            fs[k].close(); os.replace(os.path.join(outdir, f"cand_{c}.txt.tmp"), os.path.join(outdir, f"cand_{c}.txt"))
        with open(marker, "w") as f:
            f.write(json.dumps(n) + "\n")
        log(f"{tag}: supervivientes de motor C (2048 bits) por constante: " + ", ".join(f"{LABEL.get(c, c)} {n[c]}" for c in cs))
    sp = Internal(f"{tag}_split", [s2[-1].name], split, lambda S: os.path.exists(marker))
    return s1 + s2 + [sp]


def build_jobs():
    J5 = CFG["c5_parts"]
    jobs = []
    # ---- C5
    c5 = [C5Part(k) for k in range(J5)]
    jobs += c5
    qdir = P("work", "c5")

    def hc_splits(S):
        sp = {c: set() for c in CONSTS}
        for j in c5:
            for c, r in j.summaries().items():
                sp[c] |= set(r.get("splits_with_hits", []))
        return sp

    def hc_mk(S):
        sp = hc_splits(S)
        env = pyenv(); env.update({f"E_{c}": repr(CFG["E_c7"][c]) for c in CONSTS})
        env["HC_SPLITS"] = ";".join(f"{c}:{','.join(map(str, sorted(v)))}" for c, v in sp.items() if v)
        cs = [c for c, v in sp.items() if v]
        return ([PY, "-B", P("py", "c5", "hitcheck.py"), N, O, ",".join(cs), P("work", "c5", "hitcheck.json")], env, BELOW)

    def hc_needed(S):
        return any(hc_splits(S).values())
    hc = Ext("c5_hitcheck", [j.name for j in c5], hc_mk,
             lambda S: (not hc_needed(S)) or os.path.exists(P("work", "c5", "hitcheck.json")),
             lambda S: 0.0, CFG["w"]["hitcheck"], RAM["hitcheck"], resumable=False)
    jobs.append(hc)
    # ---- C7: motor C, motor B, C11 shortcut
    w7 = P(CFG["w_c7"])
    mc7 = motorc_jobs("c7_motorc", ["c7"], w7, lambda S: os.path.join(C7DIR, "c7_cand.txt"), P("work", "c7"),
                      lambda S: CONSTS)
    jobs += mc7
    mb7_out = P("work", "c7", "cert_9_0.jsonl")
    prev = mc7[-1].name
    for i, c in enumerate(CONSTS):
        def mk(S, c=c, i=i):
            W = max(1, min(CFG["verif_workers"], int((avail_gb() - 2) / RAM["verif_worker"])))
            env = pyenv(); env["FULLMAX"] = str(CFG["fullmax"])
            return ([PY, "-B", P("py", "cert", "certify_par.py"), N, O, c, repr(CFG["E_c7"][c]),
                     P("work", "c7", f"cand_{c}.txt"), W, mb7_out, w7, i], env, BELOW)
        j = Ext(f"c7_motorB_{c}", [prev], mk,
                lambda S, c=c: any(r["const"] == c for r in jsonl_rows(mb7_out)),
                None, CFG["w"]["motorB"], RAM["verif_worker"] * CFG["verif_workers"], resumable=True)
        j.progf = (lambda S, j=j, c=c: last_pct(tail(j.logf(), 20000), rf"\[progress\] {c}: \d+/\d+ paths \(([\d.]+)%\)"))
        jobs.append(j); prev = j.name
    nv7_out = P("work", "c7", "c11_shortcut_9_0.jsonl")
    for c in CONSTS:
        def mk(S, c=c):
            W = max(1, min(CFG["verif_workers"], int((avail_gb() - 2) / RAM["verif_worker"])))
            env = pyenv(); env.update({"C11_CANDS": P("work", "c7", f"cand_{c}.txt"), "C11_NOWIDEN": "1",
                                       "WORK": P("work", "c7", "c11run"), "FULLMAX": str(CFG["fullmax"])})
            return ([PY, "-B", P("py", "c11", "c11_par.py"), N, O, c, repr(CFG["E_c7"][c]), RES, W, nv7_out], env, BELOW)
        j = Ext(f"c7_c11atajo_{c}", [prev], mk, lambda S, c=c: any(r["const"] == c for r in jsonl_rows(nv7_out)),
                None, CFG["w"]["nv"], RAM["verif_worker"] * CFG["verif_workers"], resumable=True)
        j.progf = (lambda S, j=j, c=c: last_pct(tail(j.logf(), 20000), rf"\[progress\] {c}: \d+/\d+ paths \(([\d.]+)%\)"))
        jobs.append(j); prev = j.name
    shortcut_done = prev

    def solved_by_shortcut(c):
        return any(r["const"] == c and r.get("optimal_no_minus") for r in jsonl_rows(nv7_out))
    # ---- C11: motor A per constant (window of C11), motor C, verifier
    wl = open(P(CFG["w_c11"])).read().split("\n")
    nv11_out = P("work", "c11", "c11_9_0.jsonl")
    prevA = mc7[-1].name                 # HANDOFF_OPT 6: C11 motor A right after C7 motor C (motor B shares the cores)
    for i, c in enumerate(CONSTS):
        d = P("work", "c11", c)
        wfile = os.path.join(d, "w.txt")

        def mkA(S, c=c, i=i, d=d, wfile=wfile):
            os.makedirs(d, exist_ok=True)
            with open(wfile, "w") as f:
                f.write(wl[i] + "\n")
            nc5 = sum(1 for k in range(J5) if S["jobs"].get(f"c5_p{k}", {}).get("status") == "running")
            env = {"OMP_NUM_THREADS": str(max(1, CFG["cores"] - nc5)), "MOTORA_DEDUP": "1"}
            return ([P(CFG["motora"]), "multi", N, O, RES, wfile, os.path.join(d, "cand.txt")], env, BELOW)

        def okA(S, c=c, d=d):
            return re.search(r"^W0 candidates:", tail(P("work", "logs", f"c11_motorA_{c}.log"), 50000), re.M) is not None

        def progA(S, c=c):
            return last_pct(tail(P("work", "logs", f"c11_motorA_{c}.log"), 50000), r"W0 progress \d+/\d+ \(([\d.]+)%\)")
        ja = Ext(f"c11_motorA_{c}", [prevA], mkA, okA, progA, CFG["w"]["motorA"], RAM["motora"],
                 ready=lambda S, c=c: not solved_by_shortcut(c), partial=[os.path.join(d, "cand.txt")], resumable=False)
        ja.skip_if = (lambda S, c=c: solved_by_shortcut(c))
        jobs.append(ja); prevA = ja.name
        mc = motorc_jobs(f"c11_motorC_{c}", [ja.name], wfile, lambda S, d=d: os.path.join(d, "cand.txt"), d, lambda S, c=c: [c])
        for j in mc:
            j.skip_if = (lambda S, c=c: solved_by_shortcut(c))
        jobs += mc

        def mkV(S, c=c, d=d):
            W = max(1, min(CFG["verif_workers"], int((avail_gb() - 2) / RAM["verif_worker"])))
            env = pyenv(); env.update({"C11_CANDS": os.path.join(d, f"cand_{c}.txt"), "C11_NOWIDEN": "1",
                                       "WORK": os.path.join(d, "run"), "FULLMAX": str(CFG["fullmax"])})
            return ([PY, "-B", P("py", "c11", "c11_par.py"), N, O, c, repr(CFG["E_c11"][c]), RES, W, nv11_out], env, BELOW)
        jv = Ext(f"c11_verif_{c}", [mc[-1].name], mkV, lambda S, c=c: any(r["const"] == c for r in jsonl_rows(nv11_out)),
                 None, CFG["w"]["nv"], RAM["verif_worker"] * CFG["verif_workers"], resumable=True)
        jv.progf = (lambda S, j=jv, c=c: last_pct(tail(j.logf(), 20000), rf"\[progress\] {c}: \d+/\d+ paths \(([\d.]+)%\)"))
        jv.skip_if = (lambda S, c=c: solved_by_shortcut(c))
        jobs.append(jv)
    # ---- final package
    def pack(S):
        z = P("resultados_fase2.zip")
        with zipfile.ZipFile(z + ".tmp", "w", zipfile.ZIP_DEFLATED) as Z:
            for pat in ["fase2.log", "fase2_state.json", "fase2_config.json", "work/logs/*.log", "work/**/*.jsonl",
                        "work/**/*.json", "work/**/split_done.txt", "work/**/kept2_*.txt.sample", "work/c7/cand_*.txt",
                        "work/c11/*/cand_*.txt"]:
                for f in glob.glob(P(pat), recursive=True):
                    if os.path.getsize(f) < CFG["pack_max_mb"] * 2 ** 20:
                        Z.write(f, os.path.relpath(f, HERE))
        os.replace(z + ".tmp", z)
        log(f"resultados empaquetados en {z}")
    main = [j.name for j in jobs if not j.name.startswith("c5_")]
    c5n = [j.name for j in jobs if j.name.startswith("c5_")]
    TERM = ("hecho", "omitido", "fallo", "bloqueado")

    def pack_main(S):
        pack(S)
        with open(P("fase2_principal_fin.txt"), "w") as f:
            f.write(time.strftime("%d/%m/%Y %H:%M:%S") + " C7 y C11 terminados: devolver resultados_fase2.zip\n")
        log("C7 y C11 TERMINADOS: resultados_fase2.zip listo (C5 sigue de relleno si no ha acabado)")
    jp = Internal("empaquetar", [], pack_main, lambda S: os.path.exists(P("resultados_fase2.zip")))
    jp.ready = (lambda S: all(S["jobs"].get(n, {}).get("status") in TERM for n in main))
    jobs.append(jp)

    def pack_c5(S):
        z = P("resultados_c5.zip")
        with zipfile.ZipFile(z + ".tmp", "w", zipfile.ZIP_DEFLATED) as Z:
            for f in glob.glob(P("work", "c5", "*")) + glob.glob(P("work", "logs", "c5_*.log")):
                if os.path.isfile(f) and os.path.getsize(f) < CFG["pack_max_mb"] * 2 ** 20:
                    Z.write(f, os.path.relpath(f, HERE))
        os.replace(z + ".tmp", z)
        log(f"C5 terminado: resultados en {z}")
    jc = Internal("empaquetar_c5", [], pack_c5, lambda S: os.path.exists(P("resultados_c5.zip")))
    jc.ready = (lambda S: all(S["jobs"].get(n, {}).get("status") in TERM for n in c5n))
    jobs.append(jc)
    return jobs


# ---------------------------------------------------------------- state, scheduler, status
def load_state():
    if os.path.exists(STATE):
        return json.load(open(STATE))
    return {"jobs": {}, "started": time.time(), "c7_done": False}


def save_state(S):
    tmp = STATE + ".tmp"
    json.dump(S, open(tmp, "w"), indent=1)
    os.replace(tmp, STATE)


def status_block(S, jobs, write=True):
    c7p, c7d, per = c7_progress()
    W = CFG["w"]
    tot = W["c7"] + sum(j.weight for j in jobs)
    acc = W["c7"] * c7p / 100
    lines = [f"C7  motor A (en c7\\)       {c7p:5.1f} %  " + ", ".join(f"W{k} {p:.0f} %" for k, p in enumerate(per)) +
             ("  hecho" if c7d else "")]
    groups = {}
    for j in jobs:
        st = S["jobs"].get(j.name, {}).get("status", "pendiente")
        p = 100.0 if st in ("hecho", "omitido") else j.progress(S)
        acc += j.weight * p / 100
        g = re.sub(r"_(\d+|1024_\d+|2048_\d+|cat1|split)$", "", j.name)
        groups.setdefault(g, []).append((j, st, p))
    for g, lst in groups.items():
        wsum = sum(j.weight for j, _, _ in lst) or 1
        p = sum(j.weight * p for j, _, p in lst) / wsum if any(j.weight for j, _, _ in lst) else (100.0 if all(s == "hecho" for _, s, _ in lst) else 0.0)
        sts = sorted(set(s for _, s, _ in lst))
        det = "; ".join(d for d in (getattr(j, "detail", lambda S: "")(S) for j, s, _ in lst if s == "en curso") if d)
        lines.append(f"{g:27s} {p:5.1f} %  {'/'.join(sts)}" + (f"  ({det})" if det else ""))
    c5w = sum(j.weight for j in jobs if j.name.startswith("c5_"))
    c5a = sum(j.weight * (100.0 if S["jobs"].get(j.name, {}).get("status") in ("hecho", "omitido") else j.progress(S)) / 100
              for j in jobs if j.name.startswith("c5_"))
    gp = 100 * (acc - c5a) / max(1e-9, tot - c5w)
    left_h = (tot - c5w - (acc - c5a)) / CFG["cores"]
    warn = S.get("warnings", [])[-5:]
    head = (f"=== {time.strftime('%d/%m/%Y %H:%M')}  FASE 2 (C7 + C11): {gp:.1f} %, ~{left_h / 24:.1f} dias restantes (estimacion)"
            f"  |  C5 n=9 (relleno): {100 * c5a / max(1e-9, c5w):.1f} % ===")
    txt = "\n".join([head] + lines + [f"RAM libre {avail_gb():.1f} GB", "Avisos: " + ("; ".join(warn) if warn else "ninguno"), ""])
    if write:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(txt + "\n")
    return txt


def warn(S, msg):
    S.setdefault("warnings", []).append(f"{time.strftime('%d/%m %H:%M')} {msg}")
    log("AVISO: " + msg)


PIDF = P("fase2_pid.txt")


def orchestrator_alive():
    try:
        return alive(int(open(PIDF).read().split()[0]))
    except (OSError, ValueError, IndexError):
        return False


def parar():
    """stop the orchestrator and every process it launched, by PID only (never by window title)"""
    S = load_state()
    pids = []
    try:
        pids.append(int(open(PIDF).read().split()[0]))
    except (OSError, ValueError, IndexError):
        pass
    pids += [int(st["pid"]) for st in S["jobs"].values() if st.get("status") == "en curso" and st.get("pid")]
    for pid in pids:
        if alive(pid):
            kill_tree(pid); log(f"parar: proceso {pid} parado")
    log("FASE 2 parada a peticion (lo hecho se conserva; lanzar_fase2.bat reanuda)")


def run():
    if orchestrator_alive():
        print("el orquestador ya esta en marcha"); return
    with open(PIDF, "w") as f:
        f.write(f"{os.getpid()}\n")
    S = load_state()
    jobs = build_jobs()
    for d in [("work", "logs"), ("work", "c5"), ("work", "c7"), ("work", "c11")] + [("work", "c11", c) for c in CONSTS]:
        os.makedirs(P(*d), exist_ok=True)
    byname = {j.name: j for j in jobs}
    log("orquestador arrancado" + (" (reanudando)" if S["jobs"] else ""))
    procs = {}
    last_status = 0
    while True:
        c7p, c7d, _ = c7_progress()
        if c7d and not S["c7_done"]:
            S["c7_done"] = True
            log("C7 (motor A) terminado: empieza motor C de C7 (C5 sigue de relleno, con prioridad inactiva)")
        done_names = {"c7"} if S["c7_done"] else set()
        # refresh running jobs
        for j in jobs:
            st = S["jobs"].setdefault(j.name, {"status": "pendiente", "tries": 0})
            if st["status"] == "en curso":
                p = procs.get(j.name)
                rc = p.poll() if p else (None if alive(st.get("pid")) else "?")
                if rc is None:
                    continue
                st["end"] = time.time(); st["rc"] = rc
                if j.ok(S):
                    st["status"] = "hecho"; log(f"{j.name}: terminado")
                elif st["tries"] < CFG["retries"]:
                    st["status"] = "pendiente"
                    warn(S, f"{j.name} termino sin completar (rc {rc}); se " + ("reanudara" if j.resumable else "relanzara desde el principio"))
                else:
                    st["status"] = "fallo"; warn(S, f"{j.name} FALLO (rc {rc}); ver work/logs/{j.name}.log")
                procs.pop(j.name, None)
        for j in jobs:
            if S["jobs"][j.name]["status"] in ("hecho", "omitido"):
                done_names.add(j.name)
        RUNNING_RAM[0] = sum(j.ram for j in jobs if S["jobs"][j.name]["status"] == "en curso")
        # memory guard: never let Windows page; the newest C5 part yields (it resumes later)
        if avail_gb() < CFG["ram_min_free"]:
            run5 = [j for j in jobs if isinstance(j, C5Part) and S["jobs"][j.name]["status"] == "en curso"]
            if run5:
                j = max(run5, key=lambda j: S["jobs"][j.name].get("start", 0))
                kill_tree(S["jobs"][j.name]["pid"]); S["jobs"][j.name]["status"] = "pendiente"; procs.pop(j.name, None)
                S["jobs"][j.name]["tries"] -= 1
                warn(S, f"poca RAM libre ({avail_gb():.1f} GB): {j.name} parado, se reanudara")
        # launch
        starved = False
        for j in sorted(jobs, key=lambda j: isinstance(j, C5Part)):    # main jobs first, C5 (filler) last
            st = S["jobs"][j.name]
            if st["status"] != "pendiente":
                continue
            if any(S["jobs"].get(d, {}).get("status") == "fallo" for d in j.deps):
                st["status"] = "bloqueado"; warn(S, f"{j.name} bloqueado: fallo una tarea previa"); continue
            if not all(d in done_names for d in j.deps):
                continue
            if getattr(j, "skip_if", None) and j.skip_if(S):
                st["status"] = "omitido"; log(f"{j.name}: omitido (C11 resuelto por el atajo)"); done_names.add(j.name); continue
            if not j.ready(S):
                continue
            if not isinstance(j, Internal) and j.ok(S):
                st["status"] = "hecho"; done_names.add(j.name); log(f"{j.name}: ya estaba hecho"); continue
            if isinstance(j, Internal):
                try:
                    j.internal(S); st["status"] = "hecho"; done_names.add(j.name); log(f"{j.name}: hecho")
                except Exception as e:
                    st["status"] = "fallo"; warn(S, f"{j.name} FALLO: {e!r}")
                continue
            if avail_gb() < j.ram + CFG["ram_margin"]:
                if isinstance(j, C5Part):
                    continue
                run5 = [x for x in jobs if isinstance(x, C5Part) and S["jobs"][x.name]["status"] == "en curso"]
                if run5 and avail_gb() + RAM["c5"] >= j.ram + CFG["ram_margin"]:
                    x = max(run5, key=lambda x: S["jobs"][x.name].get("start", 0))
                    kill_tree(S["jobs"][x.name]["pid"]); S["jobs"][x.name]["status"] = "pendiente"
                    S["jobs"][x.name]["tries"] -= 1; procs.pop(x.name, None); RUNNING_RAM[0] -= x.ram
                    log(f"{x.name} pausado para dejar RAM a {j.name} (C5 se reanuda despues)")
                    time.sleep(10)
                starved = True
                continue
            if isinstance(j, C5Part) and starved:
                continue
            cmd, env, prio = j.build(S)
            for f in j.outputs_partial():
                if os.path.exists(f):
                    k = 1
                    while os.path.exists(f"{f}.interrumpido_{k}"):
                        k += 1
                    os.replace(f, f"{f}.interrumpido_{k}"); log(f"{j.name}: salida parcial renombrada a {os.path.basename(f)}.interrumpido_{k}")
            p = launch(cmd, j.logf(), env, prio=prio)
            procs[j.name] = p
            st.update({"status": "en curso", "pid": p.pid, "start": time.time(), "tries": st["tries"] + 1})
            RUNNING_RAM[0] += j.ram
            log(f"{j.name}: lanzado (pid {p.pid}, intento {st['tries']})")
            time.sleep(5)
        save_state(S)
        if time.time() - last_status > CFG["status_every_s"]:
            status_block(S, jobs); last_status = time.time()
        if all(S["jobs"][j.name]["status"] in ("hecho", "omitido", "fallo", "bloqueado") for j in jobs):
            status_block(S, jobs)
            bad = [j.name for j in jobs if S["jobs"][j.name]["status"] in ("fallo", "bloqueado")]
            log("FASE 2 TERMINADA" + (" con fallos: " + ", ".join(bad) if bad else " sin fallos"))
            with open(P("fase2_fin.txt"), "w") as f:
                f.write(time.strftime("%d/%m/%Y %H:%M:%S") + (" con fallos\n" if bad else " sin fallos\n"))
            save_state(S)
            return
        time.sleep(CFG["poll_s"])


def estado():
    S = load_state()
    jobs = build_jobs()
    print(status_block(S, jobs, write=False))
    run_ = [(n, s.get("pid")) for n, s in S["jobs"].items() if s.get("status") == "en curso"]
    print("procesos en curso: " + (", ".join(f"{n} (pid {p}, {'vivo' if alive(p) else 'NO VIVO'})" for n, p in run_) or "ninguno"))
    if os.path.exists(P("fase2_fin.txt")):
        print("TERMINADO: " + open(P("fase2_fin.txt")).read())


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "estado"
    if cmd == "run":
        run()
    elif cmd == "estado":
        estado()
    elif cmd == "vivo":
        sys.exit(0 if orchestrator_alive() else 1)
    elif cmd == "parar":
        parar()
    elif cmd == "autotest":
        import autotest_fase2
        sys.exit(autotest_fase2.main("--generate" in sys.argv))
