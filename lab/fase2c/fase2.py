"""
FASE 2c (laboratorio): C11 con n = 9 en cualquier orden, sin atajo. C7 ya esta certificado (fase2b).
Para cada constante, EN SERIE y un paso pesado cada vez (cada motor ya usa los 4 nucleos):
  1 motor A (ventana de C11, MOTORA_DEDUP=1, 4 hilos)
  2 motor C 1024 bits (4 procesos por bloques)   3 join
  4 motor C 2048 bits (4 procesos por bloques)   5 join
  6 verificador de C11 (c11_par.py) sobre los supervivientes
Reglas: solo escribe dentro de esta carpeta; nunca borra (renombra); prioridad baja; reanudable (relanzar continua):
motor C y el verificador siguen desde sus bloques/trozos hechos; un motor A a medias se repite (solo esa constante).
Memoria: un paso no se lanza sin la RAM que necesita (espera y lo dice en el log). El verificador se vigila: si la
RAM libre o la memoria virtual libre (commit) bajan del minimo, se para y se relanza con un proceso menos (sigue
desde su checkpoint). Los procesos se identifican por PID + hora de creacion: nunca se mata un programa ajeno.
usage: fase2.py run | estado | vivo | parar | autotest
"""
import sys, os, json, time, subprocess, hashlib, re, zipfile, glob, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WIN = os.name == "nt"
CFG = json.load(open(os.path.join(HERE, "fase2_config.json")))
P = lambda *a: os.path.normpath(os.path.join(HERE, *a))
LOG = P("fase2.log")
STATE = P("fase2_state.json")
PIDF = P("fase2_pid.txt")
FIN = P("fase2_fin.txt")
PY = sys.executable
BELOW = 0x4000
C7DIR = P(CFG["c7_dir"])
RES = P("py", "cert", CFG["res_name"])
LABEL = {"pi": "π", "gamma": "γ", "delta": "δ"}


# ---------------------------------------------------------------- system helpers (Windows and Linux)
def mem_gb():
    """(free physical RAM, free commit = what RAM + page file can still give) in GB"""
    if os.environ.get("FASE2_TEST_MEM"):                 # autotest only: "phys,commit"
        a, b = os.environ["FASE2_TEST_MEM"].split(",")
        return float(a), float(b)
    if WIN:
        import ctypes
        class M(ctypes.Structure):
            _fields_ = [("l", ctypes.c_ulong), ("load", ctypes.c_ulong), ("tot", ctypes.c_ulonglong),
                        ("av", ctypes.c_ulonglong), ("tp", ctypes.c_ulonglong), ("ap", ctypes.c_ulonglong),
                        ("tv", ctypes.c_ulonglong), ("avv", ctypes.c_ulonglong), ("ave", ctypes.c_ulonglong)]
        m = M(); m.l = ctypes.sizeof(M)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.av / 2 ** 30, m.ap / 2 ** 30
    info = {}
    for l in open("/proc/meminfo"):
        k, v = l.split(":")[0], l.split()[1]
        info[k] = int(v) / 2 ** 20
    return info.get("MemAvailable", 0.0), info.get("MemAvailable", 0.0) + info.get("SwapFree", 0.0)


_K32 = None


def _k32():
    global _K32
    if _K32 is None:
        import ctypes
        from ctypes import wintypes as W
        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.OpenProcess.argtypes = [W.DWORD, W.BOOL, W.DWORD]; k.OpenProcess.restype = W.HANDLE
        k.GetExitCodeProcess.argtypes = [W.HANDLE, ctypes.POINTER(W.DWORD)]; k.GetExitCodeProcess.restype = W.BOOL
        k.GetProcessTimes.argtypes = [W.HANDLE] + [ctypes.POINTER(W.FILETIME)] * 4; k.GetProcessTimes.restype = W.BOOL
        k.CloseHandle.argtypes = [W.HANDLE]; k.CloseHandle.restype = W.BOOL
        _K32 = (ctypes, W, k)
    return _K32


def ptoken(pid):
    """identity of a LIVE process: its creation time (None if it does not exist or has ended). With the pid it names
    one process only, so a pid that Windows reuses for another program is never taken for ours"""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return None
    if pid <= 0:
        return None
    if WIN:
        ctypes, W, k = _k32()
        h = k.OpenProcess(0x1000, False, pid)            # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return None
        try:
            code = W.DWORD()
            if not k.GetExitCodeProcess(h, ctypes.byref(code)) or code.value != 259:      # STILL_ACTIVE
                return None
            t = [W.FILETIME() for _ in range(4)]
            if not k.GetProcessTimes(h, *[ctypes.byref(x) for x in t]):
                return None
            return str((t[0].dwHighDateTime << 32) | t[0].dwLowDateTime)
        finally:
            k.CloseHandle(h)
    try:
        f = open(f"/proc/{pid}/stat").read().rsplit(")", 1)[1].split()
    except (OSError, IndexError):
        return None
    return None if f[0] in ("Z", "X") else f[19]         # field 22 of stat: start time


def ours(pid, tok):
    return tok is not None and ptoken(pid) == tok


def kill_tree(pid):
    if WIN:
        subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"], capture_output=True)
    else:
        try:
            os.killpg(os.getpgid(int(pid)), 9)
        except OSError:
            try:
                os.kill(int(pid), 9)
            except OSError:
                pass


def launch(cmd, logf, env=None):
    os.makedirs(os.path.dirname(logf), exist_ok=True)
    lf = open(logf, "a")
    lf.write(f"\n##### {time.strftime('%d/%m/%Y %H:%M:%S')} {' '.join(map(str, cmd))}\n"); lf.flush()
    e = dict(os.environ); e.update(env or {})
    e.pop("FASE2_TEST_MEM", None)
    kw = {"creationflags": 0x08000000 | BELOW} if WIN else {"preexec_fn": lambda: (os.setpgrp(), os.nice(10))}
    p = subprocess.Popen([str(c) for c in cmd], stdout=lf, stderr=subprocess.STDOUT, env=e, cwd=HERE, **kw)
    lf.close()
    return p


def _append(path, txt):
    for _ in range(20):
        try:
            with open(path, "a", encoding="utf-8") as f:
                f.write(txt)
            return
        except OSError:
            time.sleep(0.5)


def log(msg):
    line = f"{time.strftime('%d/%m/%Y %H:%M:%S')}  {msg}"
    _append(LOG, line + "\n")
    try:
        print(line, flush=True)
    except Exception:
        pass


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
    return {"PYTHONPATH": os.pathsep.join([P("pylib"), P("py", "cert"), P("py", "c11")]),
            "PYTHONHASHSEED": "0", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONIOENCODING": "utf-8",
            "PANDI_SRC": P("py", "cert"), "PANDI_RESULTS": P("data"), "OMP_NUM_THREADS": "1"}


def jsonl_rows(path):
    out = []
    if os.path.exists(path):
        for l in open(path, encoding="utf-8"):
            try:
                out.append(json.loads(l))
            except ValueError:
                pass
    return out


def last_pct(text, pat):
    m = re.findall(pat, text)
    return float(m[-1]) if m else 0.0


# ---------------------------------------------------------------- the chain of one constant
def prod_ctx():
    """the laboratory run (n = 9, any order); the autotest builds the same chain for n = 6"""
    return {"n": CFG["n"], "o": CFG["ord"], "consts": CFG["consts"], "E": CFG["E_c11"], "res": RES,
            "wlines": [l.strip() for l in open(P(CFG["w_c11"])) if l.strip()],
            "base": P("work", "c11"), "logs": P("work", "logs"), "out": P("work", "c11", f"c11_{CFG['n']}_{CFG['ord']}.jsonl"),
            "nb1": CFG["motorc_nblocks"], "nb2": CFG["motorc_nblocks2"], "parts": CFG["motorc_parts"],
            "cores": CFG["cores"], "workers": CFG["verif_workers"], "fullmax": CFG["fullmax"],
            "maxtasks": CFG["verif_maxtasks"], "need": CFG["need_gb"], "hours": CFG["hours"]}


class Step:
    resumable = True
    kind = ""

    def __init__(self, X, c, num, label):
        self.X, self.c, self.num = X, c, num
        self.name = f"c11_{c}_{num}_{label}"
        self.d = os.path.join(X["base"], c)
        self.hours = X["hours"].get(self.kind, 0.05)

    def logf(self, k=0):
        return os.path.join(self.X["logs"], f"{self.name}{'' if k == 0 else '_' + str(k)}.log")

    def need(self, S):                              # GB of free RAM needed to start
        return self.X["need"].get(self.kind, 1.0)

    def cmds(self, S):                              # -> list of (cmd, env, logfile)
        raise NotImplementedError

    def ok(self, S):
        raise NotImplementedError

    def progress(self, S):
        return 100.0 if self.ok(S) else 0.0

    def partial(self):
        return []


class MotorA(Step):
    resumable = False
    kind = "motorA"

    def __init__(self, X, c, i):
        super().__init__(X, c, 1, "motorA")
        self.i = i

    def cmds(self, S):
        os.makedirs(self.d, exist_ok=True)
        with open(os.path.join(self.d, "w.txt"), "w") as f:
            f.write(self.X["wlines"][self.i] + "\n")
        env = {"OMP_NUM_THREADS": str(self.X["cores"]), "MOTORA_DEDUP": "1"}
        return [([P(CFG["motora"]), "multi", self.X["n"], self.X["o"], self.X["res"], os.path.join(self.d, "w.txt"),
                  os.path.join(self.d, "cand.txt")], env, self.logf())]

    def ok(self, S):
        return os.path.exists(os.path.join(self.d, "cand.txt")) and \
            re.search(r"^W0 candidates:", tail(self.logf(), 50000), re.M) is not None

    def progress(self, S):
        return 100.0 if self.ok(S) else last_pct(tail(self.logf(), 50000), r"W0 progress \d+/\d+ \(([\d.]+)%\)")

    def partial(self):
        return [os.path.join(self.d, "cand.txt")]


class MotorC(Step):
    def __init__(self, X, c, num, bits, nb, src, kept):
        self.kind = "motorC1024" if bits == 1024 else "motorC2048"
        super().__init__(X, c, num, f"motorC{bits}")
        self.bits, self.nb, self.src, self.kept = bits, nb, src, kept

    def need(self, S):
        return self.X["need"]["motorC_part"] * self.X["parts"]

    def cmds(self, S):
        J = self.X["parts"]
        out = []
        for k in range(J):
            env = {"MOTORC_FAST": "1", "MOTORC_TPREC": "128", "MOTORC_SAMPLE": "1000", "MOTORC_PROGRESS": "600",
                   "MOTORC_NBLOCKS": str(self.nb), "MOTORC_PART": f"{k}/{J}"}
            out.append(([P(CFG["motorc"]), self.X["n"], self.X["o"], os.path.join(self.d, "w.txt"), self.src,
                         self.kept, self.bits], env, self.logf(k + 1)))
        return out

    def ok(self, S):
        return all(os.path.exists(f"{self.kept}.b{b:04d}.done") for b in range(self.nb))

    def progress(self, S):
        J = self.X["parts"]
        tot = 0.0
        for k in range(J):
            mine = [b for b in range(self.nb) if b % J == k]
            d = sum(os.path.exists(f"{self.kept}.b{b:04d}.done") for b in mine)
            m = re.findall(r"\[\d+\] block \d+/\d+: ([\d.]+)% read", tail(self.logf(k + 1), 20000))
            tot += d + (float(m[-1]) / 100 if (m and d < len(mine)) else 0.0)
        return min(100.0, 100.0 * tot / self.nb)


class Join(Step):
    kind = "join"

    def __init__(self, X, c, num, kept, nb):
        super().__init__(X, c, num, "join")
        self.kept, self.nb = kept, nb

    def cmds(self, S):
        return [([P(CFG["motorc"]), "join", self.kept, self.nb], {}, self.logf())]

    def ok(self, S):
        return os.path.exists(self.kept) and \
            re.search(r"^joined \d+ blocks", tail(self.logf(), 5000), re.M) is not None


class Verif(Step):
    kind = "verif"

    def __init__(self, X, c, kept):
        super().__init__(X, c, 6, "verificador")
        self.kept = kept

    def workers(self, S):
        return int(S.setdefault("workers", {}).get(self.c, self.X["workers"]))

    def need(self, S):
        return self.X["need"]["verif_base"] + self.X["need"]["verif_worker"] * self.workers(S)

    def cmds(self, S):
        env = pyenv()
        env.update({"C11_CANDS": self.kept, "C11_NOWIDEN": "1", "WORK": os.path.join(self.d, "run"),
                    "FULLMAX": str(self.X["fullmax"]), "C11_MAXTASKS": str(self.X["maxtasks"])})
        return [([PY, "-B", P("py", "c11", "c11_par.py"), self.X["n"], self.X["o"], self.c, repr(self.X["E"][self.c]),
                  self.X["res"], self.workers(S), self.X["out"]], env, self.logf())]

    def ok(self, S):
        return any(r.get("const") == self.c and r.get("n") == self.X["n"] for r in jsonl_rows(self.X["out"]))

    def progress(self, S):
        return 100.0 if self.ok(S) else \
            last_pct(tail(self.logf(), 20000), rf"\[progress\] {self.c}: \d+/\d+ paths \(([\d.]+)%\)")


def build_steps(X):
    steps = []
    for i, c in enumerate(X["consts"]):
        d = os.path.join(X["base"], c)
        k1, k2 = os.path.join(d, "kept1.txt"), os.path.join(d, "kept2.txt")
        steps += [MotorA(X, c, i),
                  MotorC(X, c, 2, 1024, X["nb1"], os.path.join(d, "cand.txt"), k1),
                  Join(X, c, 3, k1, X["nb1"]),
                  MotorC(X, c, 4, 2048, X["nb2"], k1, k2),
                  Join(X, c, 5, k2, X["nb2"]),
                  Verif(X, c, k2)]
    return steps


# ---------------------------------------------------------------- state and the main loop
def new_state():
    return {"steps": {}, "workers": {}, "warnings": [], "started": time.time()}


def load_state():
    if os.path.exists(STATE):
        for _ in range(20):
            try:
                with open(STATE) as f:
                    return json.load(f)
            except (OSError, ValueError):
                time.sleep(0.5)
        raise RuntimeError("no se puede leer fase2_state.json")
    return new_state()


def save_state(S, path=None):
    """atomic; Windows can refuse the replace for a moment (antivirus, an editor reading the file): retry"""
    path = path or STATE
    tmp = path + ".tmp"
    for _ in range(60):
        try:
            with open(tmp, "w") as f:
                json.dump(S, f, indent=1)
            os.replace(tmp, path)
            return True
        except OSError:
            time.sleep(1)
    log("aviso: no se pudo guardar el estado en 60 s (se reintenta en la siguiente vuelta)")
    return False


def warn(S, msg):
    S["warnings"].append(f"{time.strftime('%d/%m %H:%M')} {msg}")
    log("AVISO: " + msg)


def rename_partial(s):
    for f in s.partial():
        if os.path.exists(f):
            k = 1
            while os.path.exists(f"{f}.interrumpido_{k}"):
                k += 1
            os.replace(f, f"{f}.interrumpido_{k}")
            log(f"{s.name}: salida parcial renombrada a {os.path.basename(f)}.interrumpido_{k}")


def live_pids(st):
    return [pid for pid, tok in zip(st.get("pids", []), st.get("tokens", [])) if ours(pid, tok)]


def stop_step(st, procs):
    for pid in live_pids(st):
        kill_tree(pid)
    for p in procs or []:
        try:
            p.wait(timeout=60)
        except Exception:
            pass


def after_step(S, X, s, st):
    if isinstance(s, Join) and s.num == 5:
        n = sum(1 for _ in open(s.kept, "rb")) if os.path.exists(s.kept) else 0
        st["kept"] = n
        log(f"{s.name}: {n} caminos sobreviven a motor C (2048 bits) para {LABEL.get(s.c, s.c)}")
    if isinstance(s, Verif):
        r = [r for r in jsonl_rows(X["out"]) if r.get("const") == s.c and r.get("n") == X["n"]][-1]
        o = r.get("optimal_no_minus") or []
        if o:
            log(f"C11 {LABEL.get(s.c, s.c)} n = {X['n']}: optimo sin menos unario, error {o[0]['error']}, "
                f"{o[0]['formula_no_minus']}" + (f", {len(o)} valores empatados" if len(o) > 1 else "") +
                f"; sin decidir {r.get('unresolved')}")
        else:
            warn(S, f"C11 {LABEL.get(s.c, s.c)}: ningun valor sin menos unario en la ventana (no deberia pasar): avisar a Claude")


def run_once(S, X, steps, procs):
    """one turn of the main loop; True when every step has ended (done, failed or blocked)"""
    cur = next((s for s in steps if S["steps"].get(s.name, {}).get("status") == "en curso"), None)
    if cur is not None:
        st = S["steps"][cur.name]
        mine = procs.get(cur.name)
        if mine is not None:
            running = any(p.poll() is None for p in mine)
            rcs = [p.poll() for p in mine]
        else:                                         # after a restart of the orchestrator
            running = bool(live_pids(st))
            rcs = ["?"]
        if running and isinstance(cur, Verif):        # memory guard of the verifier
            ph, pc = mem_gb()
            if ph < CFG["ram_min_free"] or pc < CFG["commit_min_free"]:
                stop_step(st, mine)
                w = max(1, cur.workers(S) - 1)
                S["workers"][cur.c] = w
                st["status"] = "pendiente"; st["tries"] -= 1; procs.pop(cur.name, None)
                warn(S, f"{cur.name}: memoria baja (RAM libre {ph:.1f} GB, virtual libre {pc:.1f} GB): parado; se "
                        f"relanza con {w} proceso(s) y sigue desde su checkpoint")
                return False
        if running:
            return False
        st["end"] = time.time(); st["rc"] = rcs
        procs.pop(cur.name, None)
        if cur.ok(S):
            st["status"] = "hecho"
            log(f"{cur.name}: terminado ({(st['end'] - st['start']) / 3600:.2f} h)")
            after_step(S, X, cur, st)
        elif st["tries"] < CFG["retries"]:
            st["status"] = "pendiente"
            warn(S, f"{cur.name} termino sin completar (rc {rcs}); se " +
                 ("reanuda" if cur.resumable else "repite desde el principio"))
        else:
            st["status"] = "fallo"
            warn(S, f"{cur.name} FALLO (rc {rcs}); ver work/logs/{cur.name}*.log")
        return False
    for s in steps:                                   # the next step (a failure blocks only its own constant)
        st = S["steps"].setdefault(s.name, {"status": "pendiente", "tries": 0})
        if st["status"] in ("hecho", "fallo", "bloqueado"):
            continue
        prev = [x for x in steps[:steps.index(s)] if x.c == s.c]
        if any(S["steps"].get(x.name, {}).get("status") in ("fallo", "bloqueado") for x in prev):
            st["status"] = "bloqueado"; warn(S, f"{s.name} bloqueado: fallo un paso anterior de {s.c}")
            continue
        if s.resumable and s.ok(S):
            st["status"] = "hecho"; log(f"{s.name}: ya estaba hecho"); return False
        ph, pc = mem_gb()
        need = s.need(S)
        if ph < need or pc < need + CFG["commit_min_free"]:
            if time.time() - st.get("waitlog", 0) > CFG["status_every_s"]:
                st["waitlog"] = time.time()
                warn(S, f"{s.name} espera memoria: necesita {need:.1f} GB libres (RAM libre {ph:.1f} GB, virtual "
                        f"libre {pc:.1f} GB). Si dura, cerrar otros programas pesados")
            st["esperando"] = True
            return False
        st["esperando"] = False
        if not s.resumable:
            rename_partial(s)
        ps = []
        for cmd, env, lf in s.cmds(S):
            ps.append(launch(cmd, lf, env))
        procs[s.name] = ps
        st.update({"status": "en curso", "pids": [p.pid for p in ps], "tokens": [ptoken(p.pid) for p in ps],
                   "start": time.time(), "tries": st["tries"] + 1})
        log(f"{s.name}: lanzado (pids {[p.pid for p in ps]}, intento {st['tries']}" +
            (f", {s.workers(S)} procesos" if isinstance(s, Verif) else "") + ")")
        return False
    return True


def status_block(S, steps, write=True):
    lines, left = [], 0.0
    for s in steps:
        st = S["steps"].get(s.name, {})
        sts = st.get("status", "pendiente")
        p = 100.0 if sts == "hecho" else s.progress(S)
        if sts == "en curso":
            el = (time.time() - st.get("start", time.time())) / 3600
            left += (el * (100 - p) / p) if p > 3 else max(0.25 * s.hours, s.hours - el)
        elif sts == "pendiente":
            left += s.hours
        extra = ""
        if isinstance(s, Verif) and sts == "en curso":
            extra = f"  ({s.workers(S)} procesos)"
        if sts == "pendiente" and st.get("esperando"):
            extra = "  (esperando memoria libre)"
        if isinstance(s, Join) and "kept" in st:
            extra = f"  ({st['kept']} caminos)"
        lines.append(f"  {s.name:28s} {p:5.1f} %  {sts}{extra}")
    done = sum(1 for s in steps if S["steps"].get(s.name, {}).get("status") == "hecho")
    ph, pc = mem_gb()
    head = (f"=== {time.strftime('%d/%m/%Y %H:%M')}  FASE 2c (C11, n = {CFG['n']}): {done}/{len(steps)} pasos hechos; "
            f"quedan unas {left:.0f} h (estimacion) ===")
    txt = "\n".join([head] + lines + [f"  RAM libre {ph:.1f} GB, memoria virtual libre {pc:.1f} GB",
                                      "  Avisos: " + (" | ".join(S["warnings"][-5:]) if S["warnings"] else "ninguno"), ""])
    if write:
        _append(LOG, txt + "\n")
    return txt


def pack():
    z = P("resultados_fase2c.zip")
    pats = ["fase2.log", "fase2_state.json", "fase2_config.json", "fase2_fin.txt", "fase2_consola.log",
            "work/logs/*.log", "work/c11/*.jsonl", "work/c11/*/kept2.txt", "work/c11/*/kept*.txt.sample",
            "work/c11/*/run/*.jsonl"]
    with zipfile.ZipFile(z + ".tmp", "w", zipfile.ZIP_DEFLATED) as Z:
        for pat in pats:
            for f in sorted(glob.glob(P(pat))):
                if os.path.getsize(f) < CFG["pack_max_mb"] * 2 ** 20:
                    Z.write(f, os.path.relpath(f, HERE))
    os.replace(z + ".tmp", z)
    log(f"resultados empaquetados en {z}")


def orchestrator():
    """(pid, token) of the running orchestrator, or None"""
    try:
        pid, tok = open(PIDF).read().split()[:2]
    except (OSError, ValueError):
        return None
    return (int(pid), tok) if ours(pid, tok) else None


def parar():
    o = orchestrator()
    if o:
        kill_tree(o[0]); log(f"parar: orquestador {o[0]} parado")
    S = load_state()
    for name, st in S["steps"].items():
        if st.get("status") == "en curso":
            for pid in live_pids(st):
                kill_tree(pid); log(f"parar: {name}, proceso {pid} parado")
            st["status"] = "pendiente"; st["tries"] = max(0, st.get("tries", 1) - 1)   # a requested stop is no failure
    save_state(S)
    log("FASE 2c parada a peticion (lo hecho se conserva; lanzar_fase2.bat reanuda)")


def run():
    if orchestrator():
        print("el orquestador ya esta en marcha"); return
    with open(PIDF, "w") as f:
        f.write(f"{os.getpid()} {ptoken(os.getpid())}\n")
    X = prod_ctx()
    S = load_state()
    for k, v in new_state().items():
        S.setdefault(k, v)
    steps = build_steps(X)
    for c in X["consts"]:
        os.makedirs(os.path.join(X["base"], c), exist_ok=True)
    os.makedirs(X["logs"], exist_ok=True)
    log("orquestador arrancado" + (" (reanudando)" if S["steps"] else "") + f", pid {os.getpid()}")
    procs, last_status = {}, 0.0
    while True:
        try:
            finished = run_once(S, X, steps, procs)
            save_state(S)
            if time.time() - last_status > CFG["status_every_s"]:
                status_block(S, steps); last_status = time.time()
            if finished:
                bad = [s.name for s in steps if S["steps"][s.name]["status"] != "hecho"]
                msg = "con fallos: " + ", ".join(bad) if bad else "sin fallos"
                status_block(S, steps)
                with open(FIN, "w") as f:
                    f.write(time.strftime("%d/%m/%Y %H:%M:%S") + " " + msg + "\n")
                pack()
                log(f"FASE 2c TERMINADA {msg}: devolver resultados_fase2c.zip")
                save_state(S)
                return
        except Exception:
            log("ERROR en el orquestador (sigue): " + traceback.format_exc().replace("\n", " | ")[-1500:])
        time.sleep(CFG["poll_s"])


def estado():
    S = load_state()
    print(status_block(S, build_steps(prod_ctx()), write=False))
    for name, st in S["steps"].items():
        if st.get("status") == "en curso":
            print(f"  en curso: {name}, procesos " +
                  ", ".join(f"{p} ({'vivo' if ours(p, t) else 'NO VIVO'})" for p, t in zip(st.get("pids", []), st.get("tokens", []))))
    print("  orquestador: " + ("vivo" if orchestrator() else "NO esta en marcha"))
    if os.path.exists(FIN):
        print("  TERMINADO: " + open(FIN).read().strip())


if __name__ == "__main__":
    for _s in (sys.stdout, sys.stderr):        # a pipe on Windows is cp1252: never die on a "π" in a message
        try:
            _s.reconfigure(errors="replace")
        except Exception:
            pass
    cmd = sys.argv[1] if len(sys.argv) > 1 else "estado"
    if cmd == "run":
        run()
    elif cmd == "estado":
        estado()
    elif cmd == "vivo":
        sys.exit(0 if orchestrator() else 1)
    elif cmd == "parar":
        parar()
    elif cmd == "autotest":
        import autotest_fase2
        sys.exit(autotest_fase2.main())
    else:
        print(__doc__)
