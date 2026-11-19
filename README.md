# Certified pandigital approximations of π, γ, δ and e — code and data

This folder contains everything needed to recompute and check the certified optima of the
pandigital problem (digits 1..n, each used once; + − × / ^, parentheses, decimal point, concatenation,
unary minus; values and every subexpression 0 or in [1e-300, 1e300]). The proofs are in the appendices of the article (A: certifier; B: the constant e; C: motor X); Appendix A.10 maps each lemma to the function that implements it. Version 3.0.

Archived in Zenodo: https://doi.org/10.5281/zenodo.24356239 · Repository: https://github.com/Bendjex/pandigital-certified

## Contents

| Path | What it is |
|---|---|
| `src/motora_*.c`, `src/motora.c`, `src/build_*.sh` | Motor A (C): tables with rigorous double enclosures, window queries, candidate paths |
| `src/coremath/pow.c`, `src/roundeven_shim.c` | CORE-MATH correctly rounded `pow` (MIT licence); `roundeven` for MinGW |
| `src/motorb.py`, `src/exact.py`, `src/canon.py`, `src/exprparse.py` | Motor B (Python + python-flint/Arb): exact values, window queries, canonical forms |
| `src/resolve.py`, `src/build_res.sh` | Resolution of motor A's singular derivations by motor B (staged table build) |
| `src/verify.py`, `src/certify.py`, `src/certify_all.sh` | Path verification and certification of one cell / of a list of cells |
| `src/motorc.c` | Motor C: rigorous Arb pre-filter of candidate paths (large cells) |
| `lab/` | Windows packages for the n = 9, any order runs |
| `tests/` | Independent checks: t1 (B vs sympy), t2 (B query vs enumeration), t4 (random expressions), motor C soundness, regression |
| `results/certified.jsonl` | One JSON line per certified cell: optimum, error ball (25 digits), formulas, candidates, checks |
| `results/resolution_files.sha256` | Hashes of the resolution files used |
| `src/v1_2026-10-28/` | First versions of motors A and C |
| `src/autotest_motores/` | Tests T1–T5 of the motors, with their expected hashes |
| `c5_final/` | Motor X, final version (Theorem 2, Appendix C): code, the five checker versions, tests, the records of n = 9 in any order and how to reproduce them |
| `c5/` | Motor X as of 11 Nov 2026 (the cells n ≤ 8) and its tests |
| `c11/` | The variant without unary minus (Section 8) |
| `analysis/` | The random model of Section 9.2 and Figure 1 |
| `demostraciones_e/` | The proofs for e (Appendix B, Table of B.7) |
| `lab/` | Windows packages of the laboratory runs |
| `docs/motorx_changes.md` | Supplement to Appendix C.9 |
| `results/` | Certificates and logs (see Results below) |

## Requirements (versions used)

gcc 13 (Linux) or MinGW-w64 gcc 13 (Windows), Python 3.13, python-flint 0.9.0, mpmath 1.3.0,
sympy 1.14.0, FLINT 3.0.1 (C library, only for motor C), wine 9 (only to test the Windows build).

## Build

```
cd src
./build_linux.sh          # motora (optimised version, -mfma) and motorc; checks motora.c = head + body + tail
```
`-ffp-contract=off` is required (no contraction of operations); never use `-ffast-math`. The optimised motor A
needs a CPU with FMA (it exits with code 6 otherwise); `-mfma` only turns the explicit fused multiply-add of
CORE-MATH's `pow` into the hardware instruction, which is correctly rounded like the software one.
`build_windows.sh` gives the Windows binaries used in the laboratory, reproducibly (MinGW-w64 gcc 13.2,
GMP 6.3.0, MPFR 4.2.1, FLINT 3.0.1); their SHA-256 are in `src/autotest_motores/SHA256SUMS_windows_binaries.txt`,
and `src/autotest_motores/ESPERADO.txt` gives the tests T1-T5 that they must pass.
The first versions of motors A and C (28 Oct 2026), with which most cells of Theorem 1 were certified, are in
`src/v1_2026-10-28/`; the optimised versions reproduce their outputs (Appendix A, after Lemma 15).

## Certify a cell

```
./build_res.sh 8 1 res_8_1.txt                     # n = 8, in order: tables + resolution of singular derivations
python3 certify.py 8 1 pi 5.279748e-09 res_8_1.txt out.jsonl
```
Arguments: n, mode (1 in order, 0 any order), constant (pi, e, gamma, delta), a guess of the optimal
error (any value at least the optimum works; the record of a search is the natural choice), the
resolution file. For the large cells, see the laboratory packages `lab/fase2b` (motor C in blocks, then motor B with `FULLMAX=4`) and `lab/fase2c`.

`certify_all.sh "1:2-9 0:2-7"` certifies every cell in those ranges using `cells/guesses.txt`.

## Checks

* Regression: `python3 ../tests/regress.py new.jsonl results/certified.jsonl` must report 0 problems.
* Linux/Windows identity: `src/autotest_motores/ESPERADO.txt` lists the tests T1–T5 with the expected hashes, the same
  on Linux and on Windows.
* `tests/t3_constants.py`: the constants λ₁, λ₂ of Lemma 7 and the digits of δ.
* `tests/t4_random_expr.py n ord nwin seed` must report 0 windows with differences.

## Expressions out of range: motor X (`c5/`, Theorem 2)

Motor X proves that no expression with an intermediate value outside [1e-300, 1e300] falls in the window
of the certified optimum (and of the minus-free optimum, Section 8). It reuses motor A's tables through a
small shared library.

```
cd c5
gcc -O2 -ffp-contract=off -fopenmp -fPIC -shared -o xlib.so xlib.c ../src/coremath/pow.c -lm
python3 c5run.py n ord consts k J out.jsonl     # part k of J of a cell, one line per (constant, split), resumable
```
`c5run.py` queries each top-level split (Theorem X) and checks every value of every hit with `lvl.py`, an
independent interval arithmetic with unbounded exponents, in the same pass. With `E11_<const>` set it also
checks the window of the minus-free optimum. A cell is closed when every split is recorded with 0 undecided
values and 0 undecided operations of motor B. Tests: the `t_*.py` files (one line each in their header);
`results/c5/resultados_tests.txt` lists their results. `docs/motorx_changes.md` describes the 50 changes of the
optimised version and the versions of the checker, each with its proof and test (supplement to Appendix C.9).

### Final version and n = 9 in any order (`c5_final/`)

`c5_final/` is the code that closed n = 9 in any order (and reproduces n = 8: 381 of 381 pairs (constant, split)
with the same hits). `c5_final/README.md` explains how to build it, run its 40 tests, recompute any split and
compare it with its record, and recompute the closure from the records. Its folder names are in Spanish:

| Folder | Content |
|---|---|
| `codigo/c5/` | Motor X (`motorx.py`, `qx.py`, `ivl.py`, `xl.py`, `xlib.c`), `c5run.py`, and the final checker `lvl.py` + `hitcheck.py` (patch 7c) |
| `codigo/cert/` | The certifier files that motor X imports, and CORE-MATH. `motora.c` here is the first version (`src/v1_2026-10-28/`), as on the laboratory machine; motor X only uses it in the fingerprint of its table cache |
| `codigo/versiones_verificador/` | The checker versions that checked the records (original, patches 3, 4, 5, 7, 7c), each with its fingerprint, and `huella.py` |
| `codigo/cierre/` | The closure: `cierre.py` (`closure()`, the same function as for n ≤ 8), `cierre_final.py`, the radii |
| `codigo/datos/` | The `certified.jsonl` and `c11.jsonl` the runs used (see Known limitations) |
| `tests/` | 40 tests and `run_tests.py` |
| `resultados_n9/` | The records of n = 9, the records used by the closure, `tabla_versiones.csv` (checker of every record), the closure output and `estado_c5.csv` |
| `muestra_7c/` | Splits 9, 17, 70 and 191 (one per earlier checker) recomputed with the final checker: same hits, 0 undecided |

## Variant without unary minus (`c11/`)

| File | What it is |
|---|---|
| `nm.py` | Motor B with two flags per magnitude (+u / −u reachable without unary minus) |
| `c11_certify.py` | One cell: motor A window (×4 until a minus-free value appears), flagged verification, Arb minimum, independent check |
| `c11_par.py` | Same, with the verification split across processes and a resumable checkpoint (large cells) |
| `mitm.py` | Upper bound for the minus-free optimum (meet in the middle, checked at 80 digits) |
| `test_small.py` | Brute-force check of the flags on every set of ≤ 4 digits |

`python3 c11_certify.py n ord const E_start ../src/res_n_ord.txt out.jsonl`. Results: `results/c11.jsonl`.

## Random model (`analysis/`)

`dp.py` and `run_w.py` count, in double precision, the distinct values near each constant (heuristic,
not part of any proof); `fig_d1.py` draws Figure 1 of the article.

## Results

| File | Content |
|---|---|
| `results/certified.jsonl` | Theorem 1: one line per certified cell (58 cells: π, γ, δ for n = 2–9 in both modes, e for n ≤ 6) |
| `results/e_7_1.jsonl` | Theorem 1 for e, n = 7 in order |
| `results/c7_n9_any/` | n = 9 in any order: motor B certificates, the paths that survived motor C, motor B checkpoints |
| `results/c5/` | Theorem 2: `estado_c5.csv` (one line per cell and constant: radii, splits, hits, undecided values for C7 and C11; all 59 closed), the session records and the tests. The records of n = 9 in any order are in `c5_final/resultados_n9/` |
| `results/c11.jsonl` | Variant without unary minus: all 59 cells |
| `results/c11_n9_any/`, `results/c11_n8_delta/` | The large cells of C11: surviving paths, verifier checkpoints, log |
| `results/resolution_files.sha256` | Hashes of the resolution files (built with `OMP_NUM_THREADS=1`; with more threads the lines come out in another order, and `sort file | sha256sum` is the invariant) |
| `results/c5_cells_q.jsonl`, `results/c5_hitcheck_9_1.jsonl` | Theorem 2 with the first version of motor X (2 Nov 2026), kept for reference |

## Lab packages (`lab/`)

Windows packages for the n = 9 cells in any order. They run with the Python of the machine (3.13) plus the
pure-Python libraries of the package, and the MinGW builds of motors A and C; nothing is installed.
- `lab/` (root): motor A for C7 (one pass for the three windows, 47 h).
- `lab/fase2b/`: motor C in blocks and motor B in parallel for C7.
- `lab/fase2c/`: C11, one constant at a time (motor A → motor C → join → verifier), with memory guards and
  process identification by PID and creation time.
- The laboratory driver of motor X for n = 8 and 9 in any order is `c5_final/codigo/cierre/cierre.py`.

## Changes since version 2.0 (2 Nov 2026)

* Optimised motors A and C (`src/`), used for n = 9 in any order; the first versions are in `src/v1_2026-10-28/`.
* Optimised motor X (`c5/`), with the check of the hits in the same pass and the window of C11.
* `c11/c11_par.py`: optional recycling of worker processes (`C11_MAXTASKS`); `FULLMAX=3` bounds the memory.
* Results of n = 9 in any order (C7, C11) and of δ with n = 8 in any order (C11); laboratory packages.
* Final motor X (`c5_final/`), with the records of n = 9 in any order: Theorem 2 closed for every cell.

## Known limitations of the code (documented, not changed)

* **Motor C, resumption after a power cut.** In block mode, motor C writes a checkpoint after `fflush` but
  without `fsync`. After a crash of the operating system (not after killing the process, which is tested), a
  checkpoint could survive without the data it refers to. After such a crash, delete the `.ckpt` files of the
  blocks that were running and let them start again. The runs of the article had no interruption of this kind.
* **`MOTORA_MAXITEMS`** is a test option: if it is set, motor A enumerates only part of the work items and
  the output does not say so. It must not be set for a certificate (it was not in any run of the article).

* **Motor A, an upper bound that overflows (latent).** In `finish_real`, if the upper bound of a derivation is
  +∞ and the lower bound is in range, `code_for` returns a class that does not contain the interval (Lemma 9 needs a
  finite bound). Found by the independent review; an instrumented count over every table used by the article
  (about 1.5·10⁹ calls, `tests/revision_C10/probe_inf*.c`) found no such derivation. The guard
  `if (!isfinite(hi)) return D_SING_W;` would remove the issue without changing any table used.
* **Motor X, cache of class members (latent).** `QX.members` is keyed by the position of a class in its table; if
  a class whose logarithm encloses 0 without being the exact 1 were refined, positions would shift. A count over
  every V-table used found no such class (`tests/revision_C10/check_amb.py`).
* **Motor C, very long lines.** `decide_line` copies a candidate line into a buffer of 8192 bytes and splits it into
  at most 16 tokens; real paths have at most 7 tokens and are far shorter.
* **Order of the writings.** When an optimal value has several writings of the same length, `certify.py` may list
  them in a different order from one run to another; the certificate (value, error, uniqueness) is identical.

* **C11, negative bases.** `nm.py` accepts a power with a negative base only when the exponent is an exact rational
  integer. An instrumented rerun of all 59 cells (`tests/audit_c11_negative_bases/`) found no negative base whose
  exponent is a non-rational number whose ball may contain an integer, and reproduced every result.

* **Motor X, the radius depends on `certified.jsonl`.** `c5run.py` takes E* from `certified.jsonl` when the file
  holds the cell, and from `E_<const>` otherwise. The runs of n = 9 in any order used a file without that cell
  (`c5_final/codigo/datos/`), with the radii of `c5_final/codigo/cierre/cierre_config.json`. With the
  `results/certified.jsonl` of this package the window of γ would be narrower and its hits would not match the
  records (the result would still be valid). Use `PANDI_RESULTS=c5_final/codigo/datos`, as `c5_final/README.md` says.
* **Motor X, speed of the final checker.** On some splits the checker of patch 7c is much slower than that of patch 7
  (split 85, piece NU: 66 min against 0.3 min for π), with the same hits and 0 undecided.
* **Motor X, records without fingerprint.** The records made before 13 Nov 2026 (splits 0–16 of n = 9) have no
  fingerprint or time field; they were checked by the original checker (or patch 3, which differs only in printing).
  The closure lets any later record of the same key replace them.

The folder `tests/revision_C10/` holds the independent review of the manuscript and the package (in Spanish), its
scripts and the outputs cited in the article (including motor X for δ with n = 8 and the minus-free window).

## Status (18 Nov 2026)

* Theorem 1: all 59 cells certified.
* Theorem 2: closed for every cell (59 cells, C7 and C11), the last one, n = 9 in any order, on 18 Nov 2026.
* C11: all 59 cells certified; with Theorem 2 and the C11 window, the optima are optimal over all minus-free expressions.

## Licence and citation

Code: MIT licence (`LICENSE`); CORE-MATH's `pow.c` keeps its own MIT licence. Data: CC BY 4.0.
See `CITATION.cff`.

## Contact

Questions, corrections and better records are welcome: please open an issue in the GitHub repository.
