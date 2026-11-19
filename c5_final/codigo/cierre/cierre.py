"""
C5 (teorema 2 del motor X) para n = 8 y n = 9 en cualquier orden, constantes pi, gamma y delta: la ejecucion que
cierra las casillas, para el optimo con menos unario (C7 / certified.jsonl) y, con E11 en la configuracion, tambien
para el de la variante sin menos unario (C11): una sola consulta con la ventana mas ancha, cada valor comprobado contra
las dos.  Para cada casilla:
  1. el archivo de resolucion res_<n>_0.txt (copia con el SHA-256 comprobado)
  2. las tablas N del modo con procedencia (varios procesos, cache en work\\ncache<n>)
  3. c5run.py por partes: cada proceso lleva las separaciones i con i mod J = k, las tres constantes a la vez y cada
     separacion trozo a trozo (consulta de motor X con procedencia + comprobacion independiente lvl de cada hit)
  4. el cierre: la casilla queda CERRADA para una constante si estan todas sus separaciones, con 0 valores sin decidir
     y 0 operaciones del motor B sin decidir
Antes de la primera casilla, una prueba rapida (xlib.dll bit a bit contra Linux + un hitcheck pequeno).
Reglas: solo escribe dentro de esta carpeta; nunca borra; prioridad baja; reanudable (relanzar continua por el trozo
pendiente); si falta memoria en el equipo para el proceso que mas ocupa, lo anota y lo repite solo al final.
usage: cierre.py run | lanzar | estado | revisar | parar | enmarcha | prueba_windows | prueba_lanzador | comprobar
(13/11, parche 4) 4 partes; guardas de memoria por niveles sobre el estado real de Windows (memoria comprometible
libre, RAM disponible, crecimiento del archivo de paginacion, escrituras a disco, espacio en C:); parada limpia entre
trozos; una parte que falla se relanza; lanzar crea una tarea programada (independiente de VS Code); estado da un
veredicto TODO BIEN / REVISAR y revisar prepara el zip para el director.
"""
import sys, os, json, time, subprocess, hashlib, shutil, zipfile, glob, re
HERE = os.path.dirname(os.path.abspath(__file__))
WIN = os.name == "nt"
CFG = json.load(open(os.environ.get("CIERRE_CONFIG", os.path.join(HERE, "cierre_config.json"))))
# (13/11, parche 4) cierre_config_p4.json: values that replace those of cierre_config.json (4 parts, the memory guards);
# cierre_config.json itself is not changed
_P4 = os.environ.get("CIERRE_CONFIG_P4", os.path.join(HERE, "cierre_config_p4.json"))
if os.path.exists(_P4):
    CFG.update(json.load(open(_P4, encoding="utf-8")))
P = lambda *a: os.path.normpath(os.path.join(HERE, *a))
WORK = P(CFG.get("work", "work"))
LOG = P("cierre.log")
PY = sys.executable
CONSTS = CFG["consts"]


def log(msg):
    line = f"{time.strftime('%d/%m/%Y %H:%M:%S')}  {msg}"
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception:
        pass                                             # (13/11) a log line never stops the run
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


# ---------------------------------------------------------------- memory, CPU time (Windows API through ctypes; Linux /proc)
def avail_gb():
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


def commit_avail_gb():
    """memory the system can still commit (Windows: RAM + page file still available; Linux: MemAvailable + SwapFree):
    the limit that made the earlier measurement fail, not only the free RAM"""
    if WIN:
        import ctypes
        class M(ctypes.Structure):
            _fields_ = [("l", ctypes.c_ulong), ("load", ctypes.c_ulong), ("tot", ctypes.c_ulonglong),
                        ("av", ctypes.c_ulonglong), ("tp", ctypes.c_ulonglong), ("ap", ctypes.c_ulonglong),
                        ("tv", ctypes.c_ulonglong), ("avv", ctypes.c_ulonglong), ("ave", ctypes.c_ulonglong)]
        m = M(); m.l = ctypes.sizeof(M)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ap / 2 ** 30
    d = {}
    for l in open("/proc/meminfo"):
        k, v = l.split(":")[0], l.split()[1]
        d[k] = int(v)
    return (d.get("MemAvailable", 0) + d.get("SwapFree", 0)) / 2 ** 20


def proc_stats(pid):
    """see _proc_stats; None (no data) if the system call fails: the job itself is not affected"""
    try:
        return _proc_stats(pid)
    except Exception:
        return None


def _proc_stats(pid):
    """(memory now GB, peak memory GB, CPU seconds) of a process: Windows private bytes / peak commit; Linux RSS / HWM"""
    if WIN:
        import ctypes
        from ctypes import wintypes
        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        k = ctypes.windll.kernel32
        k.OpenProcess.restype = wintypes.HANDLE
        k.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
        k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.c_void_p] * 4
        k.CloseHandle.argtypes = [wintypes.HANDLE]
        h = k.OpenProcess(0x1000 | 0x0010, False, int(pid))
        if not h:
            return None
        try:
            pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
            if not k.K32GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb):
                return None
            ft = [wintypes.FILETIME() for _ in range(4)]
            k.GetProcessTimes(h, *[ctypes.byref(x) for x in ft])
            cpu = sum((x.dwHighDateTime << 32 | x.dwLowDateTime) for x in ft[2:]) / 1e7
            return pmc.PagefileUsage / 2 ** 30, pmc.PeakPagefileUsage / 2 ** 30, cpu
        finally:
            k.CloseHandle(h)
    try:
        st = open(f"/proc/{pid}/status").read()
        rss = int(re.search(r"VmRSS:\s+(\d+)", st).group(1)) / 2 ** 20
        hwm = int(re.search(r"VmHWM:\s+(\d+)", st).group(1)) / 2 ** 20
        f = open(f"/proc/{pid}/stat").read().split(")")[-1].split()
        cpu = (int(f[11]) + int(f[12])) / os.sysconf("SC_CLK_TCK")
        return rss, hwm, cpu
    except (OSError, AttributeError):
        return None


# ---------------------------------------------------------------- (13/11, parche 4) the state of the machine
def _memstatus():
    """Windows MEMORYSTATUSEX: (RAM available, commit available, commit limit) in GB.  The commit limit is RAM + the
    page file: it grows when Windows enlarges the page file (what filled C: with the C11 shortcuts)"""
    import ctypes
    class M(ctypes.Structure):
        _fields_ = [("l", ctypes.c_ulong), ("load", ctypes.c_ulong), ("tot", ctypes.c_ulonglong),
                    ("av", ctypes.c_ulonglong), ("tp", ctypes.c_ulonglong), ("ap", ctypes.c_ulonglong),
                    ("tv", ctypes.c_ulonglong), ("avv", ctypes.c_ulonglong), ("ave", ctypes.c_ulonglong)]
    m = M(); m.l = ctypes.sizeof(M)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
        return None, None, None
    return m.av / 2 ** 30, m.ap / 2 ** 30, m.tp / 2 ** 30


def sysdrive_free_gb():
    """free space of the system drive (C:), GB; None if unknown"""
    try:
        root = (os.environ.get("SystemDrive", "C:") + "\\") if WIN else "/"
        return shutil.disk_usage(root).free / 2 ** 30
    except Exception:
        return None


def sysmem():
    """{ram, commit_free, commit_limit, cfree}: RAM available, memory that can still be committed, the commit limit
    and the free space of C:, GB (None when a reading fails: the guards then use the others).  CIERRE_FAKE_MEM (tests
    only): a JSON file whose values replace the readings"""
    m = {"ram": None, "commit_free": None, "commit_limit": None, "cfree": None}
    try:
        if WIN:
            m["ram"], m["commit_free"], m["commit_limit"] = _memstatus()
        else:
            d = {}
            for l in open("/proc/meminfo"):
                d[l.split(":")[0]] = int(l.split()[1])
            m["ram"] = d.get("MemAvailable", 0) / 2 ** 20
            m["commit_free"] = (d.get("MemAvailable", 0) + d.get("SwapFree", 0)) / 2 ** 20
            m["commit_limit"] = (d.get("MemTotal", 0) + d.get("SwapTotal", 0)) / 2 ** 20
    except Exception:
        pass
    m["cfree"] = sysdrive_free_gb()
    fake = os.environ.get("CIERRE_FAKE_MEM")
    if fake:
        try:
            m.update({k: v for k, v in json.load(open(fake)).items() if k in m})
        except Exception:
            pass
    return m


class PageOut:
    """(13/11) pages written to disk by Windows per second (performance counter \\Memory\\Pages Output/sec, by its
    English name, so it works in a Spanish Windows): sustained writes mean that Windows is paging out for lack of
    memory.  Reading the D: tables does not count (those are reads).  If the counter cannot be opened, rate() is
    None and the other guards are used."""
    def __init__(self):
        self.err = None
        self.q = self.c = None
        self.hist = []
        if not WIN:
            self.err = "no es Windows"
            return
        try:
            import ctypes
            from ctypes import wintypes
            pdh = ctypes.windll.pdh
            pdh.PdhOpenQueryW.argtypes = [wintypes.LPCWSTR, ctypes.c_size_t, ctypes.POINTER(ctypes.c_void_p)]
            pdh.PdhAddEnglishCounterW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR, ctypes.c_size_t, ctypes.POINTER(ctypes.c_void_p)]
            pdh.PdhCollectQueryData.argtypes = [ctypes.c_void_p]
            q, c = ctypes.c_void_p(), ctypes.c_void_p()
            r = pdh.PdhOpenQueryW(None, 0, ctypes.byref(q))
            if r != 0:
                raise OSError(f"PdhOpenQueryW {r:#x}")
            r = pdh.PdhAddEnglishCounterW(q, "\\Memory\\Pages Output/sec", 0, ctypes.byref(c))
            if r != 0:
                raise OSError(f"PdhAddEnglishCounterW {r:#x}")
            pdh.PdhCollectQueryData(q)
            self.pdh, self.q, self.c = pdh, q, c
        except Exception as ex:
            self.err = repr(ex)

    def rate(self):
        """pages written per second since the previous call, or None"""
        if self.q is None:
            return None
        try:
            import ctypes
            from ctypes import wintypes
            class V(ctypes.Structure):
                _fields_ = [("CStatus", wintypes.DWORD), ("doubleValue", ctypes.c_double)]
            self.pdh.PdhGetFormattedCounterValue.argtypes = [ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p, ctypes.POINTER(V)]
            if self.pdh.PdhCollectQueryData(self.q) != 0:
                return None
            v = V()
            if self.pdh.PdhGetFormattedCounterValue(self.c, 0x00000200, None, ctypes.byref(v)) != 0 or v.CStatus not in (0, 1):
                return None
            return float(v.doubleValue)
        except Exception:
            return None

    def sustained(self, now, window):
        """the mean rate over the last `window` seconds (None until the window is covered)"""
        r = self.rate()
        fake = os.environ.get("CIERRE_FAKE_MEM")
        if fake:
            try:
                f = json.load(open(fake))
                if "pageout" in f:
                    r = float(f["pageout"])
            except Exception:
                pass
        if r is not None:
            self.hist.append((now, r))
        self.hist = [(t, x) for (t, x) in self.hist if t >= now - window]
        if not self.hist or now - self.hist[0][0] < 0.9 * window:
            return None
        return sum(x for _, x in self.hist) / len(self.hist)


def proc_table():
    """{pid: (parent pid, exe name)} of all processes (Windows: Toolhelp snapshot; Linux: /proc); {} if unknown"""
    out = {}
    try:
        if WIN:
            import ctypes
            from ctypes import wintypes
            class PE(ctypes.Structure):
                _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ProcessID", wintypes.DWORD),
                            ("th32DefaultHeapID", ctypes.c_void_p), ("th32ModuleID", wintypes.DWORD),
                            ("cntThreads", wintypes.DWORD), ("th32ParentProcessID", wintypes.DWORD),
                            ("pcPriClassBase", ctypes.c_long), ("dwFlags", wintypes.DWORD), ("szExeFile", ctypes.c_wchar * 260)]
            k = ctypes.windll.kernel32
            k.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
            k.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
            k.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PE)]
            k.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PE)]
            k.CloseHandle.argtypes = [wintypes.HANDLE]
            h = k.CreateToolhelp32Snapshot(2, 0)
            try:
                pe = PE(); pe.dwSize = ctypes.sizeof(PE)
                ok = k.Process32FirstW(h, ctypes.byref(pe))
                while ok:
                    out[int(pe.th32ProcessID)] = (int(pe.th32ParentProcessID), pe.szExeFile)
                    ok = k.Process32NextW(h, ctypes.byref(pe))
            finally:
                k.CloseHandle(h)
        else:
            for st in glob.glob("/proc/[0-9]*/stat"):
                try:
                    t = open(st).read()
                    out[int(st.split("/")[2])] = (int(t.split(")")[-1].split()[1]), t[t.index("(") + 1:t.rindex(")")])
                except (OSError, ValueError, IndexError):
                    pass
    except Exception:
        pass
    return out


def ancestors(pid):
    """the names of the ancestors of a process, nearest first"""
    t = proc_table(); out = []; seen = set()
    p = t.get(int(pid), (None, None))[0]
    while p and p in t and p not in seen and len(out) < 30:
        seen.add(p)
        out.append(t[p][1]); p = t[p][0]
    return out


def code_fingerprint():
    """the fingerprint of the checker (py/c5/lvl.py + hitcheck.py), the same as c5run.CODE"""
    h = hashlib.sha256()
    for nm in ("lvl.py", "hitcheck.py"):
        h.update(open(P("py", "c5", nm), "rb").read().replace(b"\r\n", b"\n"))
    return h.hexdigest()[:16]


def pyenv(cell):
    sep = os.pathsep
    pylib = [P(x) for x in CFG["pylib"]]
    n = cell["n"]
    env = {"PYTHONPATH": sep.join(pylib + [P("py", "c5"), P("py", "cert")]),
           "PYTHONHASHSEED": "0", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONIOENCODING": "utf-8",
           "PANDI_SRC": P("py", "cert"), "PANDI_RESULTS": P("data"), "RESDIR": os.path.join(WORK, "res"),
           "QX_NCACHE": os.path.join(WORK, f"ncache{n}"), "OMP_NUM_THREADS": "1", "QX_CHUNK": str(CFG.get("chunk", 1000000))}
    env.update({f"E_{c}": repr(float(e)) for c, e in cell.get("E", {}).items()})
    # (11/11) the optimum of the variant without unary minus (C11) of each constant: c5run checks its window too
    env.update({f"E11_{c}": repr(float(e)) for c, e in cell.get("E11", {}).items()})
    return env


def job_limit(pid, gb):
    """(14/11, parche 5) Windows: puts the process in a job object of its own with a limit on its committed memory
    (gb): an allocation above the limit fails inside that process (MemoryError, the part ends with an error and is
    relaunched or left as FALLO) instead of pushing Windows to enlarge the page file.  Returns the handle of the job
    (kept open while the process lives), or None if it could not be set (the other guards still apply)."""
    if not WIN or not gb:
        return None
    try:
        import ctypes
        from ctypes import wintypes
        k = ctypes.windll.kernel32

        class IO(ctypes.Structure):
            _fields_ = [(n, ctypes.c_ulonglong) for n in ("ro", "wo", "oo", "rt", "wt", "ot")]

        class BASIC(ctypes.Structure):
            _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                        ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                        ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                        ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD),
                        ("SchedulingClass", wintypes.DWORD)]

        class EXT(ctypes.Structure):
            _fields_ = [("BasicLimitInformation", BASIC), ("IoInfo", IO), ("ProcessMemoryLimit", ctypes.c_size_t),
                        ("JobMemoryLimit", ctypes.c_size_t), ("PeakProcessMemoryUsed", ctypes.c_size_t),
                        ("PeakJobMemoryUsed", ctypes.c_size_t)]

        k.CreateJobObjectW.restype = wintypes.HANDLE
        k.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        k.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        k.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        k.OpenProcess.restype = wintypes.HANDLE
        k.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k.CloseHandle.argtypes = [wintypes.HANDLE]
        h = k.CreateJobObjectW(None, None)
        if not h:
            return None
        info = EXT()
        info.BasicLimitInformation.LimitFlags = 0x00000100          # JOB_OBJECT_LIMIT_PROCESS_MEMORY
        info.ProcessMemoryLimit = int(gb * 2 ** 30)
        if not k.SetInformationJobObject(h, 9, ctypes.byref(info), ctypes.sizeof(info)):   # extended limits
            k.CloseHandle(h)
            return None
        ph = k.OpenProcess(0x0001 | 0x0100, False, int(pid))        # PROCESS_TERMINATE | PROCESS_SET_QUOTA
        ok = bool(ph) and bool(k.AssignProcessToJobObject(h, ph))
        if ph:
            k.CloseHandle(ph)
        if not ok:
            k.CloseHandle(h)
            return None
        return h
    except Exception:
        return None


def popen(cmd, logf, env):
    os.makedirs(os.path.dirname(logf), exist_ok=True)
    e = dict(os.environ); e.update(env)
    lf = open(logf, "a", encoding="utf-8")
    lf.write(f"\n##### {time.strftime('%d/%m/%Y %H:%M:%S')} {' '.join(map(str, cmd))}\n"); lf.flush()
    kw = {"creationflags": 0x08000000 | 0x4000} if WIN else {"preexec_fn": lambda: os.nice(10)}
    return subprocess.Popen([str(c) for c in cmd], stdout=lf, stderr=subprocess.STDOUT, env=e, cwd=HERE, **kw)


def kill(p):
    if WIN:
        subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"], capture_output=True,
                       creationflags=0x08000000)
    else:
        p.kill()
    p.wait()


def mark(name):
    return os.path.join(WORK, "hecho", name + ".txt")


def done(name):
    return os.path.exists(mark(name))


def set_done(name, text):
    os.makedirs(os.path.dirname(mark(name)), exist_ok=True)
    with open(mark(name), "w", encoding="utf-8") as f:
        f.write(text + "\n")


def n_splits(n):
    from math import comb
    return sum(comb(n, r) for r in range(1, n // 2 + 1)) - (comb(n, n // 2) // 2 if n % 2 == 0 else 0)


def guard_all(live, tag, t0):
    """the memory of all the live processes of this program {name: Popen}: appended to work\\memoria\\total.tsv; returns
    a reason to stop EVERYTHING (never to retry by itself) when their sum passes max_total_gb, None otherwise"""
    tot = 0.0; n = 0
    for p in live.values():
        s = proc_stats(p.pid)
        if s:
            tot += s[0]; n += 1
    os.makedirs(os.path.join(WORK, "memoria"), exist_ok=True)
    with open(os.path.join(WORK, "memoria", "total.tsv"), "a", encoding="utf-8") as f:
        f.write(f"{time.strftime('%d/%m %H:%M:%S')}\t{tag}\t{n}\t{tot:.2f}\t{avail_gb():.1f}\t{commit_avail_gb():.1f}\n")
    cap = CFG.get("max_total_gb", 20)
    if tot > cap:
        return f"la memoria de todos los procesos de este programa paso de {cap} GB ({tot:.1f} GB)"
    return None


def stop_all(what, reason):
    """everything has been stopped for memory: say so and leave (what is done is kept; nothing is relaunched)"""
    log(f"PARADO TODO (nivel 3) en {what}: {reason}. No se relanza solo. Lo hecho se conserva; ejecutar revisar.bat "
        f"y pasar el zip al director antes de volver a lanzarlo")
    try:
        pack()
    except Exception:
        pass
    raise SystemExit(2)


# ---------------------------------------------------------------- steps
def step_res(cell):
    n = cell["n"]
    dst = os.path.join(WORK, "res", f"res_{n}_0.txt")
    if os.path.exists(dst) and sha256(dst) == cell["res_sha256"]:
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    for src in cell["res_from"]:
        src = P(src)
        if os.path.exists(src) and sha256(src) == cell["res_sha256"]:
            shutil.copyfile(src, dst + ".tmp"); os.replace(dst + ".tmp", dst)
            log(f"res_{n}_0.txt copiado de {src} (SHA-256 correcto)")
            return
    raise SystemExit(f"FALLO: no se encuentra res_{n}_0.txt con el SHA-256 esperado en: " + ", ".join(cell["res_from"]))


def step_check():
    if done("prueba"):
        return
    log("prueba rapida: xlib.dll bit a bit y un hitcheck pequeno (unos minutos)")
    lf = os.path.join(WORK, "logs", "prueba.log")
    out1 = os.path.join(WORK, "logs", "prueba_xlibbits.txt")
    os.makedirs(os.path.dirname(out1), exist_ok=True)
    env = dict(os.environ); env.update(pyenv({"n": 5}))
    with open(out1, "w", encoding="utf-8") as f:
        rc1 = subprocess.run([PY, "-B", P("py", "c5", "t_xlibbits.py"), P("py", "cert", "res_5_0.txt")], stdout=f,
                             stderr=subprocess.STDOUT, env=env, cwd=HERE).returncode
    keep = lambda t: [l.split() for l in t.splitlines() if re.match(r"^(exp_lo|exp_hi|ln_lo|ln_hi|tables_n5) ", l)]
    got = keep(open(out1, encoding="utf-8").read())
    want = keep(open(P("data", "xlibbits_linux.txt"), encoding="utf-8").read())
    if rc1 != 0 or got != want:
        raise SystemExit(f"FALLO: xlib.dll no coincide bit a bit con Linux (ver {out1}); no se lanza nada mas")
    p = popen([PY, "-B", P("py", "c5", "t_hc_smoke.py")], lf, pyenv({"n": 5}))
    rc2 = p.wait()
    m = re.findall(r"values (\d+) inside-or-undecided (\d+) fail (\d+)", tail(lf, 20000))
    want2 = re.findall(r"values (\d+) inside-or-undecided (\d+) fail (\d+)", open(P("data", "hc_smoke_linux.txt")).read())
    if rc2 != 0 or not m or m[-1] != want2[-1]:
        raise SystemExit(f"FALLO: el hitcheck pequeno no da lo mismo que en Linux (ver {lf}); no se lanza nada mas")
    set_done("prueba", "ok")
    log("prueba rapida OK (xlib.dll y el hitcheck pequeno dan lo mismo que en Linux)")


def step_tables(cell):
    n = cell["n"]
    if done(f"tablas{n}"):
        return
    J = CFG.get("tables_workers", 3)
    log(f"tablas N de n = {n}: {J} procesos en paralelo (cache en work\\ncache{n}; si ya estaban, es rapido)")
    env = pyenv(cell); env["QX_ALLPROV"] = "1"      # the tables of the provenance mode (c5run / hitcheck)
    procs = []
    for k in range(J):
        lf = os.path.join(WORK, "logs", f"tablas{n}_{k}.log")
        procs.append((popen([PY, "-B", P("py", "c5", "ntables.py"), n, 0, 4, 5, k, J], lf, env), lf))
    t0 = time.time()
    while any(p.poll() is None for p, _ in procs):
        time.sleep(CFG.get("poll_s", 15))
        reason = guard_all({k: p for k, (p, _) in enumerate(procs) if p.poll() is None}, f"tablas{n}", t0)
        if reason:
            for p, _ in procs:
                if p.poll() is None:
                    kill(p)
            stop_all(f"tablas de n = {n}", reason)
    bad = [lf for p, lf in procs if p.returncode != 0 or not re.search(r"^done \d+ tables", tail(lf, 20000), re.M)]
    if bad:
        raise SystemExit("FALLO en las tablas N (ver " + ", ".join(bad) + "); relanzar continua")
    for p, lf in procs:
        # (12/11) the count is read from the line "done N tables ...; motor B undecided operations: K" itself, wherever
        # it is in the log: in Windows the process can write more lines after it (the check of 11/11 wanted it last)
        m = re.findall(r"^done \d+ tables.*motor B undecided operations: (\d+)", tail(lf, 20000), re.M)
        if not m or m[-1] != "0":
            raise SystemExit(f"FALLO: el motor B dejo operaciones sin decidir en las tablas (ver {lf})")
    set_done(f"tablas{n}", f"{time.time() - t0:.0f} s")
    log(f"tablas N de n = {n} terminadas ({(time.time() - t0) / 3600:.2f} h)")


def part_out(n, k):
    return os.path.join(WORK, "res_c5", f"c5_{n}_parte{k}.jsonl")


def part_files(n):
    return sorted(glob.glob(os.path.join(WORK, "res_c5", f"c5_{n}_parte*.jsonl")))


def all_records(n):
    """(13/11) every record of the cell, from the files of all its parts (whatever the number of parts that wrote them),
    oldest first: the records without a time (before parche 4) first, then by time; for the same piece or split the
    newest record is the last one (the same order as c5run.load_records)"""
    recs = []
    for fi, f in enumerate(part_files(n)):
        try:
            with open(f, encoding="utf-8") as fh:
                for li, l in enumerate(fh):
                    try:
                        r = json.loads(l)
                    except ValueError:
                        continue
                    if isinstance(r, dict) and "const" in r and "split" in r:
                        recs.append((float(r.get("t", 0)), fi, li, r))
        except OSError:
            continue
    recs.sort(key=lambda x: x[:3])
    return [r for (_, _, _, r) in recs]


def records(n, J=None):
    """(split records, piece records) of the cell: the newest record of each split and each piece"""
    sp, ch = {}, {}
    for r in all_records(n):
        if "chunk" in r:
            ch[(r["const"], r["split"], r["chunk"])] = r
        elif "S1" in r:
            sp[(r["const"], r["split"])] = r
    return sp, ch


def stale_undecided(r, code):
    """a record with undecided values made by another checker: c5run makes it again (C5RUN_REDO_UNDECIDED)"""
    return bool(r.get("undecided") or r.get("undecided_c11") or r.get("inside_c11")) and r.get("code") != code


def nube_reservados(n):
    """(14/11) the splits of n made by the cloud sessions (nube_config.json), minus those given back to the lab
    (nube_liberados.json, liberar_nube): the parts of the lab leave them out (C5RUN_SKIP)"""
    try:
        nc = json.load(open(P("nube_config.json"), encoding="utf-8"))
    except (OSError, ValueError):
        return set()
    if int(nc.get("n", 9)) != n:
        return set()
    res = {int(i) for v in nc.get("sesiones", {}).values() for i in v}
    try:
        res -= {int(i) for i in json.load(open(os.path.join(WORK, "nube_liberados.json"), encoding="utf-8"))}
    except (OSError, ValueError):
        pass
    return res


def part_complete(n, k, J, sp, code=None):
    skip = nube_reservados(n)
    return all((c, i) in sp and not (code and stale_undecided(sp[(c, i)], code))
               for c in CONSTS for i in range(n_splits(n)) if i % J == k and i not in skip)


def _write_json(path, obj):
    try:
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=1)
        os.replace(tmp, path)
    except Exception:
        pass


def _append(path, line):
    """an auxiliary log line (memory): a failure never stops the run"""
    try:
        with open(path, "a", encoding="utf-8") as f:
            f.write(line)
    except Exception:
        pass


def _count_lines(path):
    try:
        with open(path, "rb") as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


def _last_run_text(lf):
    """the text of the log of a part after its last launch header"""
    t = tail(lf, 400000)
    return t.rsplit("\n##### ", 1)[-1]


STATE_F = lambda: os.path.join(WORK, "estado_partes.json")
STOP_CODE = 3


def guards():
    g = {"workers": 4, "poll_s": 2, "mem_log_s": 10, "launch_gap_s": 30,
         "l1_commit_free_gb": 8, "l1_ram_avail_gb": 4,
         "l2_commit_free_gb": 5, "l2_ram_avail_gb": 2.5, "l2_job_gb": 14,
         "l3_commit_free_gb": 3, "l3_ram_avail_gb": 1.5, "l3_pagefile_growth_gb": 0.5, "l3_cdrive_free_gb": 1.0,
         "l3_pageout_per_s": 500, "l3_pageout_window_s": 120, "l3_pageout_ram_gb": 4, "max_total_gb": 30,
         "resume_commit_free_gb": 12, "resume_ram_avail_gb": 8, "resume_after_s": 300, "retry_after_s": 30, "l2_gap_s": 15,
         "l3_pagefile_sustain_s": 60, "job_mem_gb": 12,
         "max_launches": 25}
    g.update({k: CFG[k] for k in g if k in CFG})
    return g


_PF_SINCE = [None]


def level3(G, m, ours, base_limit, po, now=None):
    """a reason to stop EVERYTHING (Windows is paging or about to page, or C: is almost full), or None.  (14/11) The
    growth of the page file counts when it lasts l3_pagefile_sustain_s seconds: on 14/11 at 01:13 Windows enlarged it
    by 0.9 GB for a moment (and released it) with 10 GB of commit still free, and everything was stopped"""
    if m["commit_free"] is not None and m["commit_free"] < G["l3_commit_free_gb"]:
        return f"memoria comprometible libre {m['commit_free']:.1f} GB (< {G['l3_commit_free_gb']} GB)"
    if m["ram"] is not None and m["ram"] < G["l3_ram_avail_gb"]:
        return f"RAM disponible {m['ram']:.1f} GB (< {G['l3_ram_avail_gb']} GB): Windows empezaria a paginar a disco"
    grown = bool(base_limit and m["commit_limit"] is not None and m["commit_limit"] - base_limit > G["l3_pagefile_growth_gb"])
    now = time.time() if now is None else now
    if not grown:
        _PF_SINCE[0] = None
    else:
        _PF_SINCE[0] = _PF_SINCE[0] or now
        if now - _PF_SINCE[0] >= G["l3_pagefile_sustain_s"]:
            return (f"el archivo de paginacion lleva {now - _PF_SINCE[0]:.0f} s ampliado ({m['commit_limit'] - base_limit:.1f} GB "
                    f"sobre el arranque): Windows esta usando disco como memoria")
    if m["cfree"] is not None and m["cfree"] < G["l3_cdrive_free_gb"]:
        return f"C: tiene {m['cfree']:.2f} GB libres (< {G['l3_cdrive_free_gb']} GB)"
    if po is not None and po > G["l3_pageout_per_s"] and (m["ram"] is None or m["ram"] < G["l3_pageout_ram_gb"]):
        return (f"Windows escribe {po:.0f} paginas/s a disco (media de {G['l3_pageout_window_s']} s) con poca RAM "
                f"disponible: esta paginando")
    if ours > G["max_total_gb"]:
        return f"la memoria de los procesos de este programa paso de {G['max_total_gb']} GB ({ours:.1f} GB)"
    return None


def level2(G, m, st):
    """(part, reason, runaway) to stop AT ONCE, or None"""
    for k, s in st.items():
        if s[0] > G["l2_job_gb"]:
            return k, f"el proceso paso de {G['l2_job_gb']} GB ({s[0]:.1f} GB): desbocado", True
    low = []
    if m["commit_free"] is not None and m["commit_free"] < G["l2_commit_free_gb"]:
        low.append(f"memoria comprometible libre {m['commit_free']:.1f} GB")
    if m["ram"] is not None and m["ram"] < G["l2_ram_avail_gb"]:
        low.append(f"RAM disponible {m['ram']:.1f} GB")
    if low and st:
        k = max(st, key=lambda k_: st[k_][0])
        return k, "; ".join(low) + f"; se para ya la parte que mas ocupa ({st[k][0]:.1f} GB)", False
    return None


def level1(G, m):
    """a reason to stop one part CLEANLY (after its current piece), or None"""
    low = []
    if m["commit_free"] is not None and m["commit_free"] < G["l1_commit_free_gb"]:
        low.append(f"memoria comprometible libre {m['commit_free']:.1f} GB (< {G['l1_commit_free_gb']} GB)")
    if m["ram"] is not None and m["ram"] < G["l1_ram_avail_gb"]:
        low.append(f"RAM disponible {m['ram']:.1f} GB (< {G['l1_ram_avail_gb']} GB)")
    return "; ".join(low) or None


def healthy(G, m):
    return ((m["commit_free"] is None or m["commit_free"] >= G["resume_commit_free_gb"]) and
            (m["ram"] is None or m["ram"] >= G["resume_ram_avail_gb"]))


def step_parts(cell):
    """(13/11, parche 4) the J parts of the cell in parallel, watched every poll_s seconds:
      level 3 (Windows paging or about to page, C: almost full): everything is stopped and nothing is relaunched;
      level 2 (memory low, or a process above l2_job_gb): the largest part is stopped at once and resumed later;
      level 1 (memory getting low): the largest part stops cleanly after its current piece (no work lost) and is
        resumed when the machine has had room again for resume_after_s;
      a part that ends with an error is relaunched; if it fails again without having recorded anything new, it stays
        stopped (FALLO in cierre.log) and the others go on."""
    n = cell["n"]; G = guards(); J = int(cell.get("workers", G["workers"]))
    os.makedirs(os.path.join(WORK, "res_c5"), exist_ok=True); os.makedirs(os.path.join(WORK, "memoria"), exist_ok=True)
    env = pyenv(cell); env["C5RUN_CHUNKS"] = "1"
    env["C5RUN_REDO_UNDECIDED"] = "1"            # pieces left undecided by another checker are made again
    env["C5RUN_RECORDS_GLOB"] = os.path.join(WORK, "res_c5", f"c5_{n}_parte*.jsonl")   # the records of every part
    skip = nube_reservados(n)
    if skip:
        env["C5RUN_SKIP"] = ",".join(map(str, sorted(skip)))   # (14/11) made by the cloud sessions
        log(f"n = {n}: {len(skip)} separaciones reservadas a las sesiones de la nube (no las hace el laboratorio)")
    code = code_fingerprint()
    sp, _ = records(n)
    todo = [k for k in range(J) if not part_complete(n, k, J, sp, code)]
    if not todo:
        return
    log(f"n = {n}: {len(todo)} partes de {J} ({', '.join(CONSTS)} a la vez); cada trozo terminado queda guardado; "
        f"guardas: nivel 1 (parada limpia) con menos de {G['l1_commit_free_gb']} GB comprometibles o "
        f"{G['l1_ram_avail_gb']} GB de RAM; nivel 2 (parada inmediata) con menos de {G['l2_commit_free_gb']} / "
        f"{G['l2_ram_avail_gb']} GB o un proceso de mas de {G['l2_job_gb']} GB; nivel 3 (se para todo) con menos de "
        f"{G['l3_commit_free_gb']} / {G['l3_ram_avail_gb']} GB, si crece el archivo de paginacion, si Windows pagina "
        f"o si C: baja de {G['l3_cdrive_free_gb']} GB")
    t0 = time.time()
    m0 = sysmem(); base_limit = m0["commit_limit"]
    po = PageOut()
    if po.err:
        log(f"(contador de escrituras a disco no disponible: {po.err}; se usan las demas guardas)")
    parts = {k: {"state": "pending", "since": t0, "strikes": 0, "launches": 0, "recs_at_fail": None, "reason": ""} for k in todo}
    live = {}                                            # k -> (Popen, log file, stop token)
    stopfile = lambda k: os.path.join(WORK, f"parar_parte{k}.txt")
    last = {"launch": 0.0}
    jobs = {}                                            # k -> job handle (Windows), kept while the part lives

    def launch(k):
        token = f"{os.getpid()}-{k}-{time.time():.3f}"
        e = dict(env); e["C5RUN_STOP_FILE"] = stopfile(k); e["C5RUN_STOP_TOKEN"] = token
        lf = os.path.join(WORK, "logs", f"c5_{n}_parte{k}.log")
        p = popen([PY, "-B", P("py", "c5", "c5run.py"), n, 0, ",".join(CONSTS), k, J, part_out(n, k)], lf, e)
        jh = job_limit(p.pid, G["job_mem_gb"])           # (14/11) a part can never commit more than job_mem_gb
        if WIN and jh is None and not launch.warned:
            launch.warned = True
            log(f"(no se pudo poner el tope de {G['job_mem_gb']} GB por proceso; siguen las demas guardas)")
        jobs[k] = jh
        live[k] = (p, lf, token)
        parts[k].update(state="live", since=time.time(), launches=parts[k]["launches"] + 1, reason="")
        last["launch"] = time.time()
    launch.last = last
    launch.warned = False

    def request_stop(k, text):
        try:
            with open(stopfile(k), "w", encoding="utf-8") as f:
                f.write(text)
        except OSError:
            pass

    def save_state(extra=None):
        d = {"updated": time.strftime("%d/%m/%Y %H:%M:%S"), "n": n, "workers": J, "parts": parts}
        if extra:
            d.update(extra)
        _write_json(STATE_F(), d)

    try:
        _watch_loop(n, G, J, t0, parts, live, launch, request_stop, save_state, stopfile, base_limit, po)
    except BaseException:
        # (13/11) whatever stops the watching (an error, a level-3 stop, Ctrl+C), no part is left running unwatched:
        # a relaunch would otherwise start a second copy next to it
        for p, _, _ in list(live.values()):
            try:
                kill(p)
            except Exception:
                pass
        raise
    save_state()
    failed = sorted(k for k, d in parts.items() if d["state"] == "failed")
    if failed:
        raise SystemExit(f"n = {n}: partes paradas por un fallo repetido: {failed}. Pasar revisar_*.zip (revisar.bat) al "
                         f"director; relanzar lanzar_cierre.bat continua por los trozos pendientes")


def _watch_loop(n, G, J, t0, parts, live, launch, request_stop, save_state, stopfile, base_limit, po):
    """the watching loop of step_parts (see there)"""
    last_l1 = last_pause = last_l2 = 0.0
    healthy_since = None
    last_memlog = last_state = 0.0
    while True:
        now = time.time()
        # launch the pending parts one at a time (their tables load at the same moment otherwise)
        pend = [k for k, d in parts.items() if d["state"] == "pending" or (d["state"] == "retry" and now - d["since"] >= G["retry_after_s"])]
        m = sysmem()
        if pend and now - launch.last["launch"] >= G["launch_gap_s"] and level1(G, m) is None:
            launch(pend[0])
        if not live and not any(d["state"] in ("pending", "retry", "paused", "stopping") for d in parts.values()):
            break
        time.sleep(G["poll_s"])
        now = time.time()
        m = sysmem()
        st = {}
        for k, (p, lf, tok) in live.items():
            s = proc_stats(p.pid)
            if s:
                st[k] = s
        fake = os.environ.get("CIERRE_FAKE_MEM")
        if fake:                                         # tests only: the memory of some parts replaced
            try:
                for k_, gb in (json.load(open(fake)).get("jobs") or {}).items():
                    if int(k_) in st:
                        st[int(k_)] = (float(gb),) + tuple(st[int(k_)][1:])
            except Exception:
                pass
        ours = sum(s[0] for s in st.values())
        pavg = po.sustained(now, G["l3_pageout_window_s"])
        # ---- level 3: stop everything
        r3 = level3(G, m, ours, base_limit, pavg, now)
        if r3:
            for p, _, _ in live.values():
                try:
                    kill(p)
                except Exception:
                    pass
            for k in live:
                parts[k]["state"] = "stopped (nivel 3)"
            save_state({"nivel3": r3})
            stop_all(f"n = {n}", r3)
        # ---- the parts that ended
        for k in [k for k, (p, _, _) in live.items() if p.poll() is not None]:
            p, lf, tok = live.pop(k)
            txt = _last_run_text(lf)
            if p.returncode == 0 and re.search(r"^DONE ", txt, re.M):
                parts[k].update(state="done", since=now)
                log(f"n = {n}, parte {k}: terminada ({(now - t0) / 3600:.2f} h desde el arranque de esta ejecucion)")
            elif p.returncode == STOP_CODE or re.search(r"^STOPPED CLEANLY", txt, re.M):
                parts[k].update(state="paused", since=now); last_pause = now; healthy_since = None
                log(f"n = {n}, parte {k}: parada limpia por memoria al acabar su trozo ({parts[k]['reason']}); se reanudara "
                    f"cuando haya memoria")
            else:
                recs = _count_lines(part_out(n, k))
                if parts[k]["recs_at_fail"] is not None and recs > parts[k]["recs_at_fail"]:
                    parts[k]["strikes"] = 0              # it recorded new pieces since the last failure
                parts[k]["strikes"] += 1; parts[k]["recs_at_fail"] = recs
                err = re.findall(r"^\w*(?:Error|Exception)[^\n]*", txt, re.M)
                if parts[k]["strikes"] >= 2 or parts[k]["launches"] >= G["max_launches"]:
                    parts[k].update(state="failed", since=now, reason=(err[-1] if err else f"codigo {p.returncode}"))
                    log(f"n = {n}, parte {k}: FALLO: termino con error dos veces sin avanzar (codigo {p.returncode}"
                        f"{'; ' + err[-1][:200] if err else ''}); queda parada, las demas siguen; ver {lf}")
                else:
                    parts[k].update(state="retry", since=now)
                    log(f"n = {n}, parte {k}: termino con error (codigo {p.returncode}{'; ' + err[-1][:200] if err else ''}); "
                        f"se relanza en {G['retry_after_s']} s desde sus trozos guardados")
        # ---- level 2: stop one part at once
        r2 = level2(G, m, st)
        # (13/11) after a stop for low memory, l2_gap_s seconds before the next one: the memory of a stopped process takes
        # a moment to come back, and without the gap one low reading would stop several parts (a runaway process: at once)
        if r2 and (r2[2] or now - last_l2 >= G["l2_gap_s"]):
            k, reason, runaway = r2
            last_l2 = now
            p, lf, tok = live.pop(k)
            try:
                kill(p)
            except Exception:
                pass
            last_pause = now; healthy_since = None
            if runaway:
                recs = _count_lines(part_out(n, k))
                if parts[k]["recs_at_fail"] is not None and recs > parts[k]["recs_at_fail"]:
                    parts[k]["strikes"] = 0
                parts[k]["strikes"] += 1; parts[k]["recs_at_fail"] = recs
                if parts[k]["strikes"] >= 2:
                    parts[k].update(state="failed", since=now, reason=reason)
                    log(f"n = {n}, parte {k}: FALLO: {reason}, dos veces sin avanzar; queda parada, las demas siguen")
                else:
                    parts[k].update(state="paused", since=now, reason=reason)
                    log(f"n = {n}, parte {k}: PARADA inmediata ({reason}); se reanudara cuando haya memoria")
            else:
                parts[k].update(state="paused", since=now, reason=reason)
                log(f"n = {n}, parte {k}: PARADA inmediata por memoria ({reason}); se reanudara cuando haya memoria")
        # ---- level 1: one part stops cleanly after its piece
        r1 = level1(G, m)
        stopping = [k for k, d in parts.items() if d["state"] == "stopping" and k in live]
        if r1 and not stopping and len(live) > 1 and now - last_l1 >= 60 and now - last_pause >= 60:
            k = max(live, key=lambda k_: st.get(k_, (0,))[0])
            request_stop(k, live[k][2])
            parts[k].update(state="stopping", reason=r1); last_l1 = now
            log(f"n = {n}, parte {k}: memoria justa ({r1}); se parara limpiamente al acabar su trozo")
        elif stopping and r1 is None and healthy(G, m) and now - last_l1 >= G["resume_after_s"]:
            for k in stopping:                          # the memory came back before the piece ended: it goes on
                request_stop(k, "-")
                parts[k].update(state="live", reason="")
                log(f"n = {n}, parte {k}: vuelve a haber memoria; sigue sin pararse")
        # ---- resume a paused part when the machine has had room for resume_after_s
        if healthy(G, m) and r1 is None:
            healthy_since = healthy_since or now
        else:
            healthy_since = None
        paused = sorted((d["since"], k) for k, d in parts.items() if d["state"] == "paused")
        if paused and healthy_since and now - healthy_since >= G["resume_after_s"] and now - last_pause >= G["resume_after_s"]:
            k = paused[0][1]
            log(f"n = {n}, parte {k}: se reanuda (memoria comprometible libre {m['commit_free'] or 0:.1f} GB, "
                f"RAM disponible {m['ram'] or 0:.1f} GB)")
            launch(k); healthy_since = now
        # ---- memory logs and state file
        if now - last_memlog >= G["mem_log_s"]:
            last_memlog = now
            for k, s in st.items():
                last = re.findall(r"\[chunks\][^\n]*", tail(live[k][1], 4000)) if k in live else []
                _append(os.path.join(WORK, "memoria", f"c5_{n}_parte{k}.tsv"),
                        f"{now - t0:.0f}\t{s[0]:.2f}\t{s[1]:.2f}\t{s[2]:.0f}\t{(last[-1][:150] if last else '')}\n")
            _append(os.path.join(WORK, "memoria", "total.tsv"),
                    f"{time.strftime('%d/%m %H:%M:%S')}\tc5_{n}\t{len(st)}\t{ours:.2f}\t{m['ram'] or 0:.1f}\t"
                    f"{m['commit_free'] or 0:.1f}\t{m['commit_limit'] or 0:.1f}\t{m['cfree'] or 0:.2f}\t"
                    f"{'' if pavg is None else round(pavg)}\n")
        if now - last_state >= 10:
            last_state = now
            save_state({"memoria": m, "nuestra_gb": round(ours, 2)})


def closure(cell):
    """the state of the cell from its records: per constant, CERRADA iff every split is recorded with 0 undecided values
    and motor B left 0 operations undecided"""
    n = cell["n"]; J = int(cell.get("workers", CFG.get("workers", 3)))
    sp, ch = records(n, J)
    ns = n_splits(n)
    # CPU of the pieces: a piece made for several constants together is counted once (under its first constant)
    cpu = sum(r["cpu"] for k, r in ch.items() if r.get("together", [k[0]])[0] == k[0])
    out = {"n": n, "splits": ns, "cpu_h": round(cpu / 3600, 2), "consts": {}}
    for c in CONSTS:
        rs = [sp[(c, i)] for i in range(ns) if (c, i) in sp]
        d = {"splits_done": len(rs), "hits": sum(r["hits"] for r in rs), "undecided": sum(r["undecided"] for r in rs),
             "motorB_unresolved": max([r.get("motorB_unresolved", 0) for r in rs] or [0]),
             "pieces_done": sum(1 for k in ch if k[0] == c)}
        d["closed"] = d["splits_done"] == ns and d["undecided"] == 0 and d["motorB_unresolved"] == 0
        if c in cell.get("E11", {}):
            # C11: every split recorded with its window checked (E11 in the record), 0 values undecided for it
            d["E11"] = float(cell["E11"][c])
            d["splits_c11"] = sum(1 for r in rs if r.get("E11") is not None)
            d["undecided_c11"] = sum(r.get("undecided_c11") or 0 for r in rs)
            d["inside_c11"] = sum(r.get("inside_c11") or 0 for r in rs)
            d["closed_c11"] = d["splits_c11"] == ns and d["undecided_c11"] == 0 and d["motorB_unresolved"] == 0
        out["consts"][c] = d
    return out


def closure_text(cl):
    lines = [f"  n = {cl['n']}: CPU de los trozos hechos {cl['cpu_h']} h"]
    for c, d in cl["consts"].items():
        if d["closed"]:
            st = "CERRADA"
        elif d["splits_done"] < cl["splits"]:
            st = "en curso / incompleta"
        else:
            st = "ABIERTA"
        lines.append(f"  n = {cl['n']}, {c}, C7 (con menos unario): {st}: {d['splits_done']}/{cl['splits']} separaciones, "
                     f"{d['hits']} hits comprobados, {d['undecided']} valores sin decidir, motor B sin decidir {d['motorB_unresolved']}")
        if "closed_c11" in d:
            st11 = ("CERRADA" if d["closed_c11"] else
                    ("en curso / incompleta" if d["splits_c11"] < cl["splits"] else "ABIERTA"))
            lines.append(f"  n = {cl['n']}, {c}, C11 (sin menos unario, E11 = {d['E11']:.6g}): {st11}: {d['splits_c11']}/{cl['splits']} "
                         f"separaciones, {d['undecided_c11']} valores sin decidir (de ellos {d['inside_c11']} con certeza dentro de "
                         f"la ventana: hay que ver si necesitan el menos unario)")
    return "\n".join(lines)


def compare_expected(cell):
    """the hits of every (constant, split) against those of the same run in Linux, when the package brings them (n = 8):
    information only (the closing does not depend on it)"""
    f = cell.get("expected")
    if not f or not os.path.exists(P(f)):
        return
    exp = json.load(open(P(f)))
    sp, _ = records(cell["n"], int(cell.get("workers", CFG.get("workers", 3))))
    got = {f"{c} {i}": r["hits"] for (c, i), r in sp.items()}
    diff = [k for k in exp if got.get(k) != exp[k]]
    if not diff:
        log(f"n = {cell['n']}: los hits de las {len(exp)} parejas (constante, separacion) coinciden con los de Linux")
    else:
        log(f"n = {cell['n']}: AVISO: {len(diff)} parejas (constante, separacion) con distinto numero de hits que en Linux "
            f"(no afecta al cierre; anotar): " + ", ".join(f"{k}: {got.get(k)} (Linux {exp[k]})" for k in diff[:12]))


def summary():
    lines = [f"RESUMEN C5 {time.strftime('%d/%m/%Y %H:%M:%S')}"]
    for cell in CFG["cells"]:
        lines.append(closure_text(closure(cell)))
    return "\n".join(lines)


def pack():
    zp = P("resultados_cierre.zip")
    with zipfile.ZipFile(zp + ".tmp", "w", zipfile.ZIP_DEFLATED) as z:
        for f in [LOG] + glob.glob(os.path.join(WORK, "logs", "*.log")) + glob.glob(os.path.join(WORK, "logs", "*.txt")) + \
                glob.glob(os.path.join(WORK, "res_c5", "*.json*")) + glob.glob(os.path.join(WORK, "memoria", "*.tsv")) + \
                glob.glob(os.path.join(WORK, "hecho", "*.txt")) + glob.glob(os.path.join(WORK, "cierre_*.json")):
            if os.path.getsize(f) < CFG.get("pack_max_mb", 200) * 2 ** 20:
                z.write(f, os.path.relpath(f, HERE))
        z.writestr("RESUMEN.txt", summary())
    try:
        os.replace(zp + ".tmp", zp)
    except OSError:                                      # (14/11) resultados_cierre.zip open elsewhere: another name
        os.replace(zp + ".tmp", P(f"resultados_cierre_{time.strftime('%Y%m%d_%H%M%S')}.zip"))


def _created(pid):
    """creation time of a process (Windows FILETIME ticks; Linux start time in clock ticks), None if it does not exist"""
    try:
        if WIN:
            import ctypes
            from ctypes import wintypes
            k = ctypes.windll.kernel32
            k.OpenProcess.restype = wintypes.HANDLE
            k.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
            k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.c_void_p] * 4
            k.CloseHandle.argtypes = [wintypes.HANDLE]
            h = k.OpenProcess(0x1000, False, int(pid))
            if not h:
                return None
            try:
                ft = [wintypes.FILETIME() for _ in range(4)]
                if not k.GetProcessTimes(h, *[ctypes.byref(x) for x in ft]):
                    return None
                return ft[0].dwHighDateTime << 32 | ft[0].dwLowDateTime
            finally:
                k.CloseHandle(h)
        return int(open(f"/proc/{int(pid)}/stat").read().split(")")[-1].split()[19])
    except (OSError, ValueError, IndexError):
        return None


PIDF = os.path.join(WORK, "cierre_pid.json")


def stop():
    """stops ONLY this orchestrator and its children: the PID it recorded at start, and only if the process with that
    PID has the same creation time (a PID reused by another program after a restart is never touched)"""
    if not os.path.exists(PIDF):
        print("No hay registro de un cierre en marcha (work\\cierre_pid.json): no se para nada.")
        return
    d = json.load(open(PIDF))
    c = _created(d["pid"])
    if c is None or c != d["created"]:
        print(f"El proceso {d['pid']} ya no es el del cierre (terminado o PID reutilizado): no se para nada.")
        return
    if WIN:
        r = subprocess.run(["taskkill", "/PID", str(d["pid"]), "/T", "/F"], capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr.strip())
    else:
        root = int(d["pid"]); kids = {root}; grew = True
        while grew:                                      # the descendants of the orchestrator (/proc ppid)
            grew = False
            for st in glob.glob("/proc/[0-9]*/stat"):
                try:
                    f = open(st).read().split(")")[-1].split()
                    pid_, ppid_ = int(st.split("/")[2]), int(f[1])
                except (OSError, ValueError, IndexError):
                    continue
                if ppid_ in kids and pid_ not in kids:
                    kids.add(pid_); grew = True
        for pid_ in sorted(kids, reverse=True):
            try:
                os.kill(pid_, 9)
            except OSError:
                pass
    log(f"cierre: parado a peticion (PID {d['pid']}); lo hecho se conserva, lanzar_cierre.bat continua")


def running():
    if not os.path.exists(PIDF):
        return False
    d = json.load(open(PIDF))
    c = _created(d["pid"])
    return c is not None and c == d["created"] and int(d["pid"]) != os.getpid()


def split_sizes(n):
    """|S1| of each top-level split in the order of QX.table_side: sizes 1 .. n // 2, comb(n, r) splits each (half of
    them for 2 r = n)"""
    from math import comb
    out = []
    for r in range(1, n // 2 + 1):
        out += [r] * (comb(n, r) // (2 if 2 * r == n else 1))
    return out


def _state():
    try:
        return json.load(open(STATE_F(), encoding="utf-8"))
    except Exception:
        return {}


def verdict():
    """(TODO BIEN | REVISAR | TERMINADO, reasons to look at, notes)"""
    reasons, notes = [], []
    try:
        code = code_fingerprint()
    except Exception as ex:
        code = None; reasons.append(f"no se puede leer el verificador (py/c5/lvl.py, hitcheck.py): {ex!r}")
    fin = os.path.exists(P("cierre_fin.txt")); run = running()
    for cell in CFG["cells"]:
        n = cell["n"]
        sp, ch = records(n)
        bad, stale, inside = {}, {}, {}
        with_pieces = {(c, i) for (c, i, _) in ch}
        items = [((c, i, pc), r) for (c, i, pc), r in ch.items()] + \
                [((c, i, "-"), r) for (c, i), r in sp.items() if (c, i) not in with_pieces]
        for (c, i, pc), r in items:
            u = max(r.get("undecided") or 0, r.get("undecided_c11") or 0)
            ins = r.get("inside_c11") or 0
            if u:
                (stale if code and stale_undecided(r, code) else bad).setdefault((i, pc), {})[c] = u
            if ins:
                inside.setdefault((i, pc), {})[c] = ins
        fmt = lambda d: ", ".join(f"split {i} trozo {pc} (" + ", ".join(f"{c} {u}" for c, u in sorted(v.items())) + ")"
                                  for (i, pc), v in sorted(d.items())[:12]) + (f" y {len(d) - 12} trozos mas" if len(d) > 12 else "")
        if bad:
            reasons.append(f"n = {n}: valores sin decidir con el verificador actual: {fmt(bad)}")
        if inside:
            reasons.append(f"n = {n}: valores CON CERTEZA DENTRO de la ventana de C11: {fmt(inside)}")
        if stale:
            notes.append(f"n = {n}: {len(stale)} trozos con sin decidir de un verificador anterior; se rehacen solos con el "
                         f"actual (no hay que hacer nada): {fmt(stale)}")
    st = _state()
    now = time.time()
    for k, d in (st.get("parts") or {}).items():
        if d.get("state") == "failed":
            reasons.append(f"parte {k} parada por un fallo repetido: {d.get('reason', '')[:200]}")
        elif d.get("state") == "paused" and now - d.get("since", now) > 6 * 3600 and run:
            reasons.append(f"parte {k} pausada por memoria desde hace {(now - d['since']) / 3600:.1f} h ({d.get('reason', '')[:150]})")
        elif d.get("state") == "paused":
            notes.append(f"parte {k} pausada por memoria; se reanudara sola ({d.get('reason', '')[:150]})")
    if st.get("nivel3"):
        reasons.append(f"PARADO TODO (nivel 3): {st['nivel3']}")
    if not run and not fin:
        last = [l for l in tail(LOG, 60000).splitlines() if "PARADO" in l or "FALLO" in l]
        reasons.append("el cierre NO esta en marcha y no ha terminado" + (f" (ultimo aviso: {last[-1][:300]})" if last else ""))
    if run:
        for k, d in (st.get("parts") or {}).items():
            if d.get("state") in ("live", "stopping"):
                for cell in CFG["cells"]:
                    lf = os.path.join(WORK, "logs", f"c5_{cell['n']}_parte{k}.log")
                    if cell["n"] == st.get("n") and os.path.exists(lf):
                        age = now - os.path.getmtime(lf)
                        if age > 4 * 3600:
                            reasons.append(f"parte {k}: su log no cambia desde hace {age / 3600:.1f} h")
    cf = sysdrive_free_gb()
    if cf is not None and cf < 2:
        reasons.append(f"C: tiene {cf:.2f} GB libres: cerrar y abrir VS Code (el calculo escribe en D:, pero con C: lleno "
                       f"Windows no puede agrandar el archivo de paginacion; con menos de 1 GB se para todo)")
    head = "TERMINADO" if fin and not reasons else ("REVISAR" if reasons else "TODO BIEN")
    return head, reasons, notes


def progress_lines(cell):
    """splits done by size of the table side, in the last 24 h, and a rough estimate of the time left"""
    n = cell["n"]; G = guards(); J = int(cell.get("workers", G["workers"]))
    sp, ch = records(n)
    sizes = split_sizes(n); ns = len(sizes)
    done = [i for i in range(ns) if all((c, i) in sp for c in CONSTS)]
    now = time.time()
    last24 = sum(1 for i in done if max(float(sp[(c, i)].get("t", 0)) for c in CONSTS) > now - 86400)
    by = {}
    for r_ in sorted(set(sizes)):
        idx = [i for i in range(ns) if sizes[i] == r_]
        by[r_] = (sum(1 for i in idx if i in done), len(idx))
    cpu = {}
    for (c, i, pc), r in ch.items():
        if r.get("together", [c])[0] == c:
            cpu[i] = cpu.get(i, 0.0) + float(r.get("cpu", 0))
    est, unknown = 0.0, []
    skip = nube_reservados(n)
    for r_, (d, tot) in by.items():
        xs = [cpu[i] for i in done if sizes[i] == r_ and i in cpu]
        left_lab = sum(1 for i in range(ns) if sizes[i] == r_ and i not in done and i not in skip)
        if not left_lab:
            continue
        if len(xs) >= 3:
            est += left_lab * sum(xs) / len(xs)
        else:
            unknown.append(r_)
    nube_line = []
    if skip:
        nube_line = [f"  n = {n}: reservadas a la nube {len(skip)} separaciones; de ellas ya importadas y comprobadas "
                     f"{sum(1 for i in skip if i in done)} (importar_nube.bat)"]
    lines = [f"  n = {n}: {len(done)}/{ns} separaciones terminadas (" +
             ", ".join(f"|S1| = {r_}: {d}/{tot}" for r_, (d, tot) in by.items()) + f"); en las ultimas 24 h: {last24}"] + nube_line
    if len(done) < ns:
        if est > 0:
            lines.append(f"  n = {n}: tiempo restante orientativo: {est / 3600 / J:.0f} h de reloj con {J} partes"
                         + (f", sin contar las de |S1| = {', '.join(map(str, unknown))} (aun sin medidas)" if unknown else "")
                         + "; puede cambiar bastante")
        else:
            lines.append(f"  n = {n}: tiempo restante: aun sin medidas suficientes para estimarlo")
    return lines


def estado():
    try:
        head, reasons, notes = verdict()
    except Exception as ex:
        head, reasons, notes = "REVISAR", [f"no se pudo calcular el veredicto: {ex!r}"], []
    print(f"VEREDICTO: {head}")
    for r in reasons:
        print(f"  REVISAR: {r}")
    for r in notes:
        print(f"  nota: {r}")
    print(summary())
    st = _state()
    for cell in CFG["cells"]:
        try:
            for l in progress_lines(cell):
                print(l)
        except Exception as ex:
            print(f"  (progreso no disponible: {ex!r})")
        n = cell["n"]
        if st.get("n") != n:
            continue
        for k, d in sorted((st.get("parts") or {}).items()):
            mem = os.path.join(WORK, "memoria", f"c5_{n}_parte{k}.tsv")
            last = tail(mem, 3000).strip().splitlines() if os.path.exists(mem) else []
            t, now_, pk, cpu, prog = ((last[-1].split("\t") if last else []) + [""] * 5)[:5]
            print(f"  n = {n}, parte {k}: {d.get('state')}"
                  + (f"; memoria {now_} GB (pico {pk} GB), CPU {float(cpu or 0) / 3600:.2f} h" if d.get("state") in ("live", "stopping") else "")
                  + (f"; {d.get('reason')[:120]}" if d.get("reason") else "")
                  + (f"\n      {prog}" if prog and d.get("state") in ("live", "stopping") else ""))
    m = sysmem()
    print(f"  equipo: RAM disponible {m['ram'] or 0:.1f} GB; memoria comprometible libre {m['commit_free'] or 0:.1f} GB de "
          f"{m['commit_limit'] or 0:.1f} GB; C: libre {m['cfree'] or 0:.2f} GB; en marcha: {'SI' if running() else 'NO'}")
    if head == "REVISAR":
        print("  -> ejecutar revisar.bat y pasar al director esta salida y el zip revisar_<fecha>.zip que crea")


def revisar():
    """revisar_<fecha>.zip: what the director needs when estado says REVISAR (read only for the calculation)"""
    import io, contextlib
    zp = P(f"revisar_{time.strftime('%Y%m%d_%H%M')}.zip")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            estado()
        except Exception as ex:
            print(f"(estado fallo: {ex!r})")
    with zipfile.ZipFile(zp + ".tmp", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("estado.txt", buf.getvalue())
        for f, nb in [(LOG, 400000), (P("cierre_consola.log"), 200000), (P("cierre_inicio.txt"), 20000), (STATE_F(), 200000)]:
            if os.path.exists(f):
                z.writestr(os.path.basename(f), tail(f, nb))
        for cell in CFG["cells"]:
            lines = []
            for r in all_records(cell["n"]):
                if r.get("undecided") or r.get("undecided_c11") or r.get("inside_c11"):
                    lines.append(json.dumps(r))
            z.writestr(f"registros_sin_decidir_{cell['n']}.jsonl", "\n".join(lines) + "\n")
        for f in glob.glob(os.path.join(WORK, "logs", "c5_*_parte*.log")):
            z.writestr("logs/" + os.path.basename(f), tail(f, 300000))
        for f in glob.glob(os.path.join(WORK, "memoria", "*.tsv")):
            z.writestr("memoria/" + os.path.basename(f), tail(f, 200000))
        hs = []
        for f in sorted(glob.glob(P("py", "c5", "*.py"))) + [P("cierre.py")]:
            hs.append(f"{sha256(f)}  {os.path.relpath(f, HERE)}")
        z.writestr("SHA256_codigo.txt", "\n".join(hs) + "\n")
        for f in (P("cierre_config.json"), _P4):
            if os.path.exists(f):
                z.write(f, os.path.basename(f))
    os.replace(zp + ".tmp", zp)
    print(f"Creado {zp} ({os.path.getsize(zp)} bytes). Pasarlo al director junto con la salida de estado_cierre.bat.")


# ---------------------------------------------------------------- (13/11) launching as a scheduled task (outside VS Code)
TASK = "C5cierre"


def _psq(s):
    return "'" + str(s).replace("'", "''") + "'"


def _pythonw():
    w = os.path.join(os.path.dirname(sys.executable), "pythonw.exe")
    return w if os.path.exists(w) else None


NOWIN = {"creationflags": 0x08000000} if WIN else {}       # (13/11) no console window for helper commands


def schedule_and_start(name, exe, args):
    """registers (or replaces) a scheduled task of the current user that runs `exe args` (no time limit), and starts it.
    The process is created by the Task Scheduler, not by VS Code: closing VS Code does not stop it.  (ok, message)"""
    script = ("$ErrorActionPreference = 'Stop'; "
              f"$a = New-ScheduledTaskAction -Execute {_psq(exe)} -Argument {_psq(args)} -WorkingDirectory {_psq(HERE)}; "
              "$s = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -AllowStartIfOnBatteries "
              "-DontStopIfGoingOnBatteries -MultipleInstances IgnoreNew -Priority 7; "
              f"Register-ScheduledTask -TaskName {_psq(name)} -Action $a -Settings $s -Force | Out-Null; "
              f"Start-ScheduledTask -TaskName {_psq(name)}; 'TAREA OK'")
    try:
        r = subprocess.run(["powershell.exe", "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-Command", script],
                           capture_output=True, text=True, timeout=180, **NOWIN)
    except Exception as ex:
        return False, repr(ex)
    return (r.returncode == 0 and "TAREA OK" in r.stdout), (r.stdout + r.stderr).strip()[-1500:]


def unschedule(name):
    try:
        subprocess.run(["powershell.exe", "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-Command",
                        f"Unregister-ScheduledTask -TaskName {_psq(name)} -Confirm:$false"], capture_output=True, text=True,
                       timeout=120, **NOWIN)
    except Exception:
        pass


def lanzar():
    """(13/11) launches the closing in the background, as a scheduled task of Windows (independent of VS Code)"""
    if os.path.exists(P("cierre_fin.txt")):
        print("Ya terminado: ver cierre_fin.txt y resultados_cierre.zip"); return 1
    if running():
        print("Ya esta en marcha (estado_cierre.bat)."); return 1
    with open(P("cierre_inicio.txt"), "a", encoding="utf-8") as f:
        f.write(f"Lanzado {time.strftime('%d/%m/%Y %H:%M:%S')} (cierre.py lanzar)\n")
    if WIN:
        pw = _pythonw()
        if pw:      # no window: cierre.py run_bg writes its console output to cierre_consola.log
            ok, msg = schedule_and_start(TASK, pw, f'-B "{P("cierre.py")}" run_bg')
        else:       # a console window (minimised by nobody): do not close it
            ok, msg = schedule_and_start(TASK, os.environ.get("ComSpec", "cmd.exe"),
                                         f'/s /c ""{sys.executable}" -B "{P("cierre.py")}" run >> "{P("cierre_consola.log")}" 2>&1"')
        if not ok:
            print("FALLO: no se pudo crear o arrancar la tarea programada:\n" + msg)
            print("No se ha lanzado nada. (Alternativa manual: lanzar_cierre_directo.bat; entonces no cerrar VS Code.)")
            return 2
    else:
        with open(P("cierre_consola.log"), "a") as lf:
            subprocess.Popen([sys.executable, "-B", P("cierre.py"), "run"], stdout=lf, stderr=subprocess.STDOUT,
                             cwd=HERE, start_new_session=True, env=dict(os.environ))
    for _ in range(60):
        time.sleep(2)
        if running():
            pid = json.load(open(PIDF))["pid"]
            anc = ancestors(pid)
            print(f"EN MARCHA (PID {pid}). Lanzado por: {' <- '.join(anc[:4]) or '?'}"
                  + (" -- AVISO: cuelga de VS Code" if any(a.lower() in ("code.exe",) for a in anc) else
                     " (independiente de VS Code)"))
            print("Estado: estado_cierre.bat   Log: cierre.log")
            return 0
    print("FALLO: no consta en marcha tras 2 minutos; ver cierre_consola.log y cierre.log")
    return 3


def prueba_windows():
    """(13/11) dry run on the lab machine: reads every measure the guards use and prints it; launches nothing"""
    ok = True
    m = sysmem()
    print(f"RAM disponible: {m['ram']} GB; memoria comprometible libre: {m['commit_free']} GB; limite de memoria "
          f"comprometida (RAM + archivo de paginacion): {m['commit_limit']} GB; C: libre: {m['cfree']} GB")
    for k in ("ram", "commit_free", "commit_limit", "cfree"):
        if m[k] is None or not m[k] > 0:
            ok = False; print(f"FALLO: lectura {k} = {m[k]}")
    if m["commit_limit"] and m["ram"] and m["commit_limit"] < m["ram"]:
        ok = False; print("FALLO: el limite de memoria comprometida es menor que la RAM disponible")
    po = PageOut()
    if po.err:
        print(f"AVISO: contador de escrituras a disco no disponible ({po.err}); se usaran las demas guardas")
    else:
        time.sleep(1.5); r1 = po.rate(); time.sleep(1.5); r2 = po.rate()
        print(f"escrituras a disco (paginas/s): {r1}, {r2}")
        if r1 is None or r2 is None:
            print("AVISO: el contador de escrituras a disco no da valores; se usaran las demas guardas")
    s = proc_stats(os.getpid())
    print(f"memoria de este proceso (GB, pico, CPU s): {s}")
    if not s:
        ok = False; print("FALLO: no se lee la memoria de un proceso")
    anc = ancestors(os.getpid())
    print(f"antecesores de este proceso: {' <- '.join(anc[:6]) or '?'}")
    if not proc_table():
        ok = False; print("FALLO: no se lee la tabla de procesos")
    try:
        print(f"huella del verificador: {code_fingerprint()}")
    except Exception as ex:
        ok = False; print(f"FALLO: {ex!r}")
    g = guards()
    print(f"guardas: {json.dumps(g)}")
    try:
        p4 = json.load(open(_P4, encoding="utf-8"))
        miss = [k for k in g if k not in p4 or p4[k] != g[k]]
        if miss:
            ok = False; print(f"FALLO: guardas que no vienen de {os.path.basename(_P4)}: {miss}")
        if g["workers"] not in (3, 4) or g["max_total_gb"] != 30:
            ok = False; print("FALLO: se esperaban 3 o 4 partes y un tope de 30 GB")
    except Exception as ex:
        ok = False; print(f"FALLO: no se lee {_P4}: {ex!r}")
    if WIN:
        # (14/11) the memory limit per process: a child limited to 1 GB must fail to take 2 GB and manage 0.5 GB
        code = ("import sys, time\ntime.sleep(3)\ntry:\n    b = bytearray(int(sys.argv[1]) * 2 ** 20)\n"
                "    print('ALLOC OK')\nexcept MemoryError:\n    print('ALLOC MEMORYERROR')\n")
        res = {}
        for mb in (2048, 512):
            try:
                c = subprocess.Popen([sys.executable, "-c", code, str(mb)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, creationflags=0x08000000)
                jh = job_limit(c.pid, 1.0)                  # the child waits 3 s before allocating
                out, _ = c.communicate(timeout=120)
                res[mb] = (jh is not None, out.strip()[-60:])
            except Exception as ex:
                res[mb] = (False, repr(ex))
        print(f"tope de memoria por proceso: 2048 MB con tope 1 GB -> {res.get(2048)}; 512 MB -> {res.get(512)}")
        if not (res.get(2048, (0, ""))[0] and "MEMORYERROR" in res[2048][1] and "ALLOC OK" in res.get(512, (0, ""))[1]):
            ok = False; print("FALLO: el tope de memoria por proceso no funciona como se espera")
    print("PRUEBA WINDOWS: OK" if ok else "PRUEBA WINDOWS: FALLO")
    return 0 if ok else 1


def prueba_tarea(out):
    """run by the scheduled test task: writes its PID and ancestors, waits 20 s, writes the end"""
    _write_json(out, {"pid": os.getpid(), "ancestros": ancestors(os.getpid()), "inicio": time.time()})
    time.sleep(20)
    d = json.load(open(out)); d["fin"] = time.time(); _write_json(out, d)


def prueba_lanzador():
    """(13/11) tries the launcher with a harmless task (C5cierre_prueba) that only writes a file and waits 20 s: the task
    must start, must not hang from VS Code and must live its 20 s; the test task is then removed"""
    os.makedirs(WORK, exist_ok=True)
    out = os.path.join(WORK, "prueba_lanzador.json")
    _write_json(out, {"estado": "esperando"})
    name = TASK + "_prueba"
    if WIN:
        pw = _pythonw() or sys.executable
        ok, msg = schedule_and_start(name, pw, f'-B "{P("cierre.py")}" prueba_tarea "{out}"')
        if not ok:
            print("PRUEBA LANZADOR: FALLO al crear la tarea programada:\n" + msg); return 1
    else:
        subprocess.Popen([sys.executable, "-B", P("cierre.py"), "prueba_tarea", out], cwd=HERE, start_new_session=True)
    d = {}
    for _ in range(60):
        time.sleep(1)
        try:
            d = json.load(open(out))
        except Exception:
            d = {}
        if "pid" in d:
            break
    for _ in range(45):
        if "fin" in d:
            break
        time.sleep(1)
        try:
            d = json.load(open(out))
        except Exception:
            pass
    if WIN:
        unschedule(name)
    anc = d.get("ancestros", [])
    vs = [a for a in anc if a.lower() in ("code.exe",)]
    print(f"tarea de prueba: PID {d.get('pid')}; lanzada por: {' <- '.join(anc[:5]) or '?'}; vivio sus 20 s: {'fin' in d}")
    ok = "pid" in d and "fin" in d and not vs
    print("PRUEBA LANZADOR: OK" if ok else "PRUEBA LANZADOR: FALLO" + (" (cuelga de VS Code)" if vs else ""))
    return 0 if ok else 1


def comprobar():
    """(13/11) the tests of parche 4 with this machine's Python (a few minutes; the closing must be stopped): the new
    rules of lvl and hitcheck, emptying the caches, and c5run (normal run, clean stop, records of other parts,
    pieces made again only with another checker) on one split of n = 8"""
    os.makedirs(os.path.join(WORK, "prueba_p4"), exist_ok=True)
    lf = os.path.join(WORK, "logs", "comprobar_p4.log")
    os.makedirs(os.path.dirname(lf), exist_ok=True)
    n9 = [c for c in CFG["cells"] if c["n"] == 9]
    env9 = dict(os.environ); env9.update(pyenv(n9[0] if n9 else CFG["cells"][-1]))
    env5 = dict(os.environ); env5.update(pyenv({"n": 5}))
    results = []

    def run(name, cmd, env, timeout, check):
        t1 = time.time()
        with open(lf, "a", encoding="utf-8") as f:
            f.write(f"\n##### {time.strftime('%d/%m/%Y %H:%M:%S')} {name}: {' '.join(map(str, cmd))}\n"); f.flush()
            try:
                r = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env,
                                   cwd=HERE, timeout=timeout, text=True, encoding="utf-8", errors="replace")
                out, rc = r.stdout, r.returncode
            except subprocess.TimeoutExpired as ex:
                out, rc = (ex.stdout or "") if isinstance(ex.stdout, str) else "", "timeout"
            f.write(out)
        ok, why = check(rc, out)
        results.append((name, ok, why, time.time() - t1))
        print(f"  {'OK   ' if ok else 'FALLO'} {name} ({time.time() - t1:.0f} s){'' if ok else ': ' + why}", flush=True)
        return ok

    # (13/11) c5run on one cheap split of n = 8 (the last one, |S1| = 4: its tables are in work\ncache8): n = 5 in any
    # order is not usable (motor A leaves derivations of size 4 unresolved with res_5_0)
    n8 = [c for c in CFG["cells"] if c["n"] == 8]
    env8 = dict(os.environ); env8.update(pyenv(n8[0] if n8 else {"n": 8}))
    J8 = n_splits(8)
    c5 = lambda *a: [PY, "-B", P("py", "c5", "c5run.py"), 8, 0, "pi,gamma,delta", *a]
    ok_problems = lambda rc, out: ((rc == 0 and re.search(r"PROBLEMS 0\b", out, re.I) is not None), f"codigo {rc}; " + out[-300:])
    run("t_lvllin (sumas que se cancelan)", [PY, "-B", P("py", "c5", "t_lvllin.py"), 300, 17], env5, 1800, ok_problems)
    run("t_lvltwc (potencias de potencias)", [PY, "-B", P("py", "c5", "t_lvltwc.py"), 800, 13], env5, 1800, ok_problems)
    run("t_lvltw (formas canonicas)", [PY, "-B", P("py", "c5", "t_lvltw.py"), 600], env5, 1800, ok_problems)
    run("t_lvloff (X / (X +- c), parche 5)", [PY, "-B", P("py", "c5", "t_lvloff.py"), 200], env5, 1800, ok_problems)
    run("t_lvlanq (potencias de racionales enormes, parche 5)", [PY, "-B", P("py", "c5", "t_lvlanq.py"), 60], env5, 1800,
        ok_problems)
    want = re.findall(r"values (\d+) inside-or-undecided (\d+) fail (\d+)", open(P("data", "hc_smoke_linux.txt")).read())
    run("t_hc_smoke (hitcheck pequeno, igual que en Linux)", [PY, "-B", P("py", "c5", "t_hc_smoke.py")], env5, 1800,
        lambda rc, out: (rc == 0 and re.findall(r"values (\d+) inside-or-undecided (\d+) fail (\d+)", out)[-1:] == want[-1:],
                         f"codigo {rc}; " + out[-300:]))
    if n9:
        run("t_hcleaf (formas exactas de los miembros de tabla, n = 9)", [PY, "-B", P("py", "c5", "t_hcleaf.py"), 0], env9, 3600, ok_problems)
        run("t_release (vaciado de caches, n = 9)", [PY, "-B", P("py", "c5", "t_release.py"), 4], env9, 3600, ok_problems)
        run("t_exact1 (base exactamente 1 frente a ventana con signo, parche 5)", [PY, "-B", P("py", "c5", "t_exact1.py"), 3000],
            env9, 3600, ok_problems)
    d = os.path.join(WORK, "prueba_p4")
    tag = time.strftime("%Y%m%d_%H%M%S")
    ra, rb, rc_, rd = (os.path.join(d, f"{x}_{tag}.jsonl") for x in "abcd")
    e = dict(env8); e["C5RUN_CHUNKS"] = "1"
    def has_done(rc, out):
        return (rc == 0 and re.search(r"^DONE ", out, re.M) is not None), f"codigo {rc}; " + out[-300:]
    if run("c5run n = 8, ultima separacion (normal)", c5(J8 - 1, J8, ra), e, 3600, has_done):
        recs = [json.loads(l) for l in open(ra)]
        good = recs and all("t" in r and "code" in r for r in recs)
        results.append(("registros con hora y huella", bool(good), "faltan campos t/code", 0))
        print(f"  {'OK   ' if good else 'FALLO'} registros con hora y huella ({len(recs)} registros)")
        # a clean stop
        tok = "prueba-" + tag
        sf = os.path.join(d, f"parar_{tag}.txt")
        with open(sf, "w") as f:
            f.write(tok)
        eb = dict(e); eb["C5RUN_STOP_FILE"] = sf; eb["C5RUN_STOP_TOKEN"] = tok
        run("c5run parada limpia", c5(J8 - 1, J8, rb), eb, 1800,
            lambda rc, out: (rc == STOP_CODE and "STOPPED CLEANLY" in out, f"codigo {rc}; " + out[-300:]))
        # the records of another part: nothing is made again
        ec = dict(e); ec["C5RUN_RECORDS_GLOB"] = ra
        run("c5run con los registros de otra parte (no repite nada)", c5(J8 - 1, J8, rc_), ec, 1800,
            lambda rc, out: (rc == 0 and "DONE" in out and _count_lines(rc_) == 0, f"codigo {rc}; lineas nuevas {_count_lines(rc_)}"))
        # a piece left undecided by another checker is made again; by this checker, not
        first = next(r for r in recs if "chunk" in r)
        key = (first["const"], first["split"])
        recs2 = []
        for r0 in recs:
            r = dict(r0)
            if (r["const"], r["split"]) == key and ("chunk" not in r or r.get("chunk") == first["chunk"]):
                r["undecided"] = 3; r["code"] = "viejo"     # the piece and its split, as a record of the old checker
            recs2.append(r)
        with open(rd, "w") as f:
            f.write("".join(json.dumps(r) + "\n" for r in recs2))
        ed = dict(e); ed["C5RUN_REDO_UNDECIDED"] = "1"
        n0 = _count_lines(rd)
        run("c5run rehace un trozo de otro verificador", c5(J8 - 1, J8, rd), ed, 1800,
            lambda rc, out: (rc == 0 and _count_lines(rd) > n0, f"codigo {rc}; lineas nuevas {_count_lines(rd) - n0}"))
        n1 = _count_lines(rd)
        run("c5run no lo rehace otra vez", c5(J8 - 1, J8, rd), ed, 1800,
            lambda rc, out: (rc == 0 and _count_lines(rd) == n1, f"codigo {rc}; lineas nuevas {_count_lines(rd) - n1}"))
    # (14/11, parche 5) the query on the deepest split of n = 8 (|S1| = 1: recursion over 7 digits, chunks of item 42,
    # bases exactly 1 of item 41): the same hits as the reference run (data/esperado_8.json, also the lab's result)
    ex8 = P("data", "esperado_8.json")
    if os.path.exists(ex8):
        want8 = json.load(open(ex8))
        re_ = os.path.join(d, f"e_{tag}.jsonl")
        def same_as_ref(rc, out):
            if rc != 0 or not re.search(r"^DONE ", out, re.M):
                return False, f"codigo {rc}; " + out[-300:]
            got = {f"{r['const']} {r['split']}": r["hits"] for r in map(json.loads, open(re_)) if "S1" in r}
            bad_ = [k for k in got if want8.get(k) != got[k]]
            return (len(got) == 3 and not bad_), f"hits {got} frente a {[(k, want8.get(k)) for k in got]}"
        run("c5run n = 8, separacion 0 (|S1| = 1): mismos hits que la referencia", c5(0, J8, re_), e, 3600, same_as_ref)
    bad = [r for r in results if not r[1]]
    print(f"COMPROBACION PARCHE 5: {'OK' if not bad else 'FALLO'} ({len(results) - len(bad)}/{len(results)} pruebas bien; "
          f"detalle en {lf})")
    return 0 if not bad else 1


def importar_nube(args):
    """(14/11) the results of the cloud sessions (resultados_nube_XX.zip in the folder nube_resultados, or the paths
    given): checked by verificar_nube.py (code, radii, every split of each session with 0 undecided, the control split)
    and, only if all is right, copied to work/res_c5 as c5_9_parte_nubeXX.jsonl.  Nothing is copied otherwise."""
    zips = args or sorted(glob.glob(P("nube_resultados", "resultados_nube_*.zip")))
    if not zips:
        print("no hay resultados: dejar los resultados_nube_XX.zip en la carpeta nube_resultados"); return 1
    lf = os.path.join(WORK, "logs", "importar_nube.log")
    os.makedirs(os.path.dirname(lf), exist_ok=True)
    r = subprocess.run([PY, "-B", P("verificar_nube.py"), HERE] + zips + ["--copiar", os.path.join(WORK, "res_c5")],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", cwd=HERE)
    with open(lf, "a", encoding="utf-8") as f:
        f.write(f"\n##### {time.strftime('%d/%m/%Y %H:%M:%S')} {zips}\n{r.stdout}")
    print(r.stdout)
    log(f"importar_nube: {'OK' if r.returncode == 0 else 'FALLO'} ({len(zips)} archivos); detalle en {lf}")
    return r.returncode


def liberar_nube(args):
    """(14/11) the splits of a cloud session that will not finish (e.g. the session was lost), given back to the lab:
    liberar_nube XX (all the splits of session XX not yet imported) or liberar_nube XX i j ... (only those)"""
    nc = json.load(open(P("nube_config.json"), encoding="utf-8"))
    ses = int(args[0])
    want = [int(x) for x in args[1:]] or nc["sesiones"][str(ses)]
    sp, _ = records(int(nc.get("n", 9)))
    want = [i for i in want if not all((c, i) in sp for c in CONSTS)]
    fn = os.path.join(WORK, "nube_liberados.json")
    try:
        old = json.load(open(fn, encoding="utf-8"))
    except (OSError, ValueError):
        old = []
    new = sorted(set(old) | set(want))
    _write_json(fn, new)
    print(f"sesion {ses:02d}: {len(want)} separaciones devueltas al laboratorio: {want}. Las hara al relanzar el cierre "
          f"(lanzar_cierre.bat).")
    log(f"liberar_nube {ses:02d}: {want}")
    return 0


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "estado"
    if cmd == "parar":
        stop()
        return
    if cmd == "enmarcha":
        print("SI" if running() else "NO")
        return
    if cmd == "estado":
        estado()
        return
    if cmd == "revisar":
        revisar()
        return
    if cmd == "lanzar":
        raise SystemExit(lanzar())
    if cmd == "prueba_windows":
        raise SystemExit(prueba_windows())
    if cmd == "prueba_lanzador":
        raise SystemExit(prueba_lanzador())
    if cmd == "prueba_tarea":
        prueba_tarea(sys.argv[2])
        return
    if cmd == "comprobar":
        if running():
            raise SystemExit("el cierre esta en marcha: comprobar solo con el cierre parado")
        raise SystemExit(comprobar())
    if cmd == "importar_nube":
        raise SystemExit(importar_nube(sys.argv[2:]))
    if cmd == "liberar_nube":
        raise SystemExit(liberar_nube(sys.argv[2:]))
    if cmd == "run_bg":
        # (13/11) run from the scheduled task with pythonw (no window): the console output goes to cierre_consola.log
        f = open(P("cierre_consola.log"), "a", encoding="utf-8", buffering=1)
        sys.stdout = sys.stderr = f
        cmd = "run"
    if cmd != "run":
        raise SystemExit(__doc__)
    os.makedirs(WORK, exist_ok=True)
    if running():
        raise SystemExit("el cierre ya esta en marcha (estado_cierre.bat): no se lanza otro")
    with open(PIDF, "w") as f:
        json.dump({"pid": os.getpid(), "created": _created(os.getpid())}, f)
    _write_json(STATE_F(), {"updated": time.strftime("%d/%m/%Y %H:%M:%S"), "arranque": True})
    log("cierre: arranca" + (" (reanudando)" if os.path.isdir(os.path.join(WORK, "hecho")) else ""))
    try:
        for cell in CFG["cells"]:
            step_res(cell)
        step_check()
        for cell in CFG["cells"]:
            n = cell["n"]
            step_tables(cell)
            step_parts(cell)
            cl = closure(cell)
            with open(os.path.join(WORK, f"cierre_{n}.json"), "w", encoding="utf-8") as f:
                json.dump(cl, f, indent=1)
            log(f"n = {n}: todas las partes terminadas\n" + closure_text(cl))
            compare_expected(cell)
    except SystemExit as ex:
        # (12/11) a stop of the closing is also written in cierre.log (before, only in the console log)
        if ex.code not in (None, 0):
            log(f"PARADO: {ex.code}")
        raise
    except Exception as ex:
        import traceback
        log("PARADO por un error inesperado: " + "".join(traceback.format_exception_only(type(ex), ex)).strip())
        raise
    finally:
        try:
            pack()
        except Exception as ex:
            log(f"(no se pudo empaquetar los resultados: {ex!r})")
    with open(P("cierre_fin.txt"), "w", encoding="utf-8") as f:
        f.write(summary() + "\n")
    log("cierre: TERMINADO. Devolver resultados_cierre.zip\n" + summary())


if __name__ == "__main__":
    main()
