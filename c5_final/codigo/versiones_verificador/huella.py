"""Fingerprint of a checker version, as written in the field "code" of every record (c5run.py, _code()):
the first 16 hex digits of SHA-256( lvl.py + hitcheck.py ), the two files concatenated in that order, with CRLF
line ends turned into LF.   usage: python huella.py <folder with lvl.py and hitcheck.py> ..."""
import hashlib, os, sys
for d in sys.argv[1:]:
    h = hashlib.sha256()
    for nm in ("lvl.py", "hitcheck.py"):
        h.update(open(os.path.join(d, nm), "rb").read().replace(b"\r\n", b"\n"))
    print(h.hexdigest()[:16], d)
