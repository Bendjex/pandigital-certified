"""ctypes wrapper of xlib.so: motor A tables as ln-enclosures, rigorous vector exp/ln bounds"""
import ctypes, numpy as np, os
_L = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "xlib.dll" if os.name == "nt" else "xlib.so"))
_dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_up = np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS")
_L.x_init.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
_L.x_count.argtypes = [ctypes.c_int]; _L.x_count.restype = ctypes.c_long
_L.x_dump.argtypes = [ctypes.c_int, _dp, _dp, _up]
_ip = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
_L.x_rats.argtypes = [ctypes.c_int, _ip, _ip]
def rats(mask):
    n = _L.x_count(mask)
    p = np.zeros(n, np.int64); q = np.zeros(n, np.int64)
    if n: _L.x_rats(mask, p, q)
    return p, q
for f in ("v_exp_lo", "v_exp_hi", "v_ln_lo", "v_ln_hi", "v_pred", "v_succ"):
    getattr(_L, f).argtypes = [_dp, _dp, ctypes.c_long]

def init(nd, ordered, res, upto):
    _L.x_init(nd, int(ordered), res.encode(), upto)

def table(mask):
    n = _L.x_count(mask)
    lo = np.empty(n); hi = np.empty(n); fl = np.empty(n, np.uint8)
    if n: _L.x_dump(mask, lo, hi, fl)
    return lo, hi, fl

def _v(name, x):
    x = np.ascontiguousarray(x, dtype=np.float64); y = np.empty_like(x)
    getattr(_L, name)(x, y, x.size)
    return y
exp_lo = lambda x: _v("v_exp_lo", x)
exp_hi = lambda x: _v("v_exp_hi", x)
ln_lo = lambda x: _v("v_ln_lo", x)
ln_hi = lambda x: _v("v_ln_hi", x)
pred = lambda x: _v("v_pred", x)
succ = lambda x: _v("v_succ", x)
