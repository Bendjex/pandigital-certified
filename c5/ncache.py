"""persistent cache of motor X tables, as one .npy file per field in a directory per table, read back with
mmap (read-only): every process reads the same pages of the operating system's file cache, so parallel C5
processes share the memory of the tables.  Only input/output: the content is exactly the computed arrays
(the cache name, made in motorx.ncache_path, holds a fingerprint of everything the table depends on).
A table is complete when its file DONE exists (written last, after an atomic rename of the directory)."""
import os, shutil, numpy as np

MMAP = os.environ.get("QX_MMAP", "1") == "1"


def load(d, fields, extra=()):
    """dict field -> array (memory-mapped, read-only) of a complete table directory, or None"""
    if not os.path.exists(os.path.join(d, "DONE")):
        return None
    out = {}
    for k in tuple(fields) + tuple(extra):
        f = os.path.join(d, k + ".npy")
        if not os.path.exists(f):
            if k in extra:
                continue
            return None
        out[k] = np.load(f, mmap_mode="r" if MMAP else None)
    return out


def save(d, arrays):
    """write the arrays (dict name -> array) as the table directory d, atomically; if another process has
    written the same table meanwhile, its copy is kept (both are the same deterministic computation)"""
    if os.path.exists(os.path.join(d, "DONE")):
        return
    tmp = f"{d}.tmp{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    for k, a in arrays.items():
        np.save(os.path.join(tmp, k + ".npy"), np.ascontiguousarray(a))
    with open(os.path.join(tmp, "DONE"), "w") as f:
        f.write("ok\n")
    try:
        os.replace(tmp, d)
    except OSError:                       # the directory exists already (another process): keep that one
        shutil.rmtree(tmp, ignore_errors=True)
