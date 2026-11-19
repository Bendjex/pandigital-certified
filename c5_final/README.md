# Motor X (Theorem 2): how to reproduce

Motor X proves Theorem 2. It shows that no expression with an intermediate value outside [10⁻³⁰⁰, 10³⁰⁰] comes within
the certified optimum of π, γ or δ. The motor runs a query that returns every enclosure (hit) that may lie in the
window, and the independent evaluator `lvl.py` then proves every value of every hit to lie outside the window.

## Requirements

The versions below are the ones the closing runs used.

| | Linux sessions (201 splits of n = 9, tests) | Windows machine (54 splits of n = 9, the closure) |
|---|---|---|
| Python | 3.13.18 | 3.13.7 |
| numpy | 2.5.3 | 2.2.6 |
| python-flint | 0.9.0 | 0.9.0 |
| mpmath | 1.3.0 | 1.3.0 |
| sympy | 1.14.0 | 1.14.0 |
| C compiler | gcc 13.3 (`-fopenmp`) | mingw-w64 gcc |

You also need the `src/` folder of the certifier of Theorem 1 (same repository). It holds the resolution files
`res_<n>_<o>.txt`; motor X uses `res_9_0.txt`, SHA-256 `a04480cf…`, as listed in the certifier's
`results/resolution_files.sha256`.

Set `PANDI_RESULTS` to `codigo/datos`, never to the certifier's `results/`. `codigo/datos` holds the `certified.jsonl`
(SHA-256 `5768bb6f…`) and `c11.jsonl` that the runs used. The reason:

- `c5run.py` takes the radius E* from `certified.jsonl` when the file holds the cell, and from `E_<const>` otherwise.
- The runs' file has no cell for n = 9 in any order, so they used the radii `E_*` given below.
- With a file that holds that cell, the window of γ would shrink from 2.12·10⁻¹⁵ to 1.98·10⁻¹⁵ and its hits would no
  longer match the records. Every record states its radius in its field `Eq`.

`codigo/cert/` is a copy of the certifier files that motor X imports. Each file is byte-identical to the certifier's.
`codigo/cert/coremath/` is CORE-MATH's correctly rounded `pow`, under the MIT license.

## Build

```bash
cd codigo/c5
gcc -O2 -ffp-contract=off -fopenmp -fPIC -shared -o xlib.so xlib.c ../cert/coremath/pow.c -I ../cert -lm
```

On Windows, build `xlib.dll` instead, adding `../cert/roundeven_shim.c`:

```bash
x86_64-w64-mingw32-gcc -O2 -ffp-contract=off -fopenmp -shared -o xlib.dll xlib.c ../cert/coremath/pow.c ../cert/roundeven_shim.c -I ../cert -lm
```

Check the build with `python tests/run_tests.py quick --only t_xlibbits`. Its hashes must be:

```text
exp_lo 5fe9f9e981066623
exp_hi e85f7546a001322c
ln_lo 7098318d80f0ebec
ln_hi 77cf15352f8a7dd7
tables_n5 263217 0b54bfc4dc710db3
```

The Linux and Windows builds give the same bits.

## Tests

```bash
export PANDI_SRC=/path/to/certifier/src PANDI_RESULTS=$PWD/codigo/datos
python tests/run_tests.py quick        # rules of lvl and of the query; about 1-2 h
QX_NCACHE=/path/to/work/ncache9 python tests/run_tests.py n9   # tests on the tables of n = 9 (cache of several GB)
```

Each test passes when it exits with code 0 and reports no problem. The summary is written to
`tests/resultados_tests.txt` and each test's full output to `tests/logs/`.

## Recompute one split of n = 9 (any order) and compare it with its record

The two windows use the radii of the closing, from `codigo/cierre/cierre_config.json`: the window with unary minus
(`E_*`) and the window without unary minus (`E11_*`).

```bash
export PYTHONPATH=$PWD/codigo/c5:$PWD/codigo/cert PANDI_SRC=/path/to/certifier/src \
       PANDI_RESULTS=$PWD/codigo/datos RESDIR=/path/to/certifier/src QX_NCACHE=$PWD/work/ncache9 \
       PYTHONHASHSEED=0 OMP_NUM_THREADS=1 QX_CHUNK=1000000 C5RUN_CHUNKS=1 \
       E_pi=3.23e-14 E_gamma=2.12e-15 E_delta=2.2e-14 \
       E11_pi=3.532867681348361e-13 E11_gamma=1.980173745532652e-15 E11_delta=1.4966108821940592e-12
python codigo/c5/c5run.py 9 0 pi,gamma,delta 70 255 split70.jsonl   # split number 70 (of 0..254)
python codigo/comparar_con_registro.py split70.jsonl resultados_n9/registros_usados_n9.jsonl
```

The last line of the comparison must read `COMPARISON: OK`.

- **Time.** A split takes from minutes to about 4 CPU hours. The first split also builds the NE tables of n = 9 in
  `QX_NCACHE`, which takes about 1 hour and a few GB.
- **Memory.** The splits with |S1| = 2 need up to about 12 GB.
- **Other checkers.** Every record carries the fingerprint of its checker (`resultados_n9/tabla_versiones.csv`). To
  repeat a record with its own checker, copy `lvl.py` and `hitcheck.py` from `codigo/versiones_verificador/<version>/`
  into `codigo/c5/`. Check that the fingerprint matches with
  `python codigo/versiones_verificador/huella.py codigo/c5`.

## Recompute the closure from the records

```bash
mkdir -p codigo/cierre/work/res_c5 && cp resultados_n9/registros/*.jsonl codigo/cierre/work/res_c5/
cd codigo/cierre && python cierre_final.py
```

The output must end with `CERRADA (C7 y C11, las tres constantes)`, and `work/cierre_9.json` must equal
`resultados_n9/cierre_9.json`. The closure only reads the records:

- For each key, it takes the newest record.
- A constant is closed when all 255 splits are recorded with 0 values undecided, 0 values certainly inside the C11
  window and 0 undecided operations of motor B.
- CPU is the sum over the pieces, with a piece made for the three constants together counted once.

## Time for the whole of n = 9

The records add up to 336 CPU hours for the three constants together:

- 89 h on the Windows PC (Intel Core i5-7600, 4 cores, 32 GB), for 54 splits;
- 247 h on 4-core, 16 GB Linux containers, for 201 splits.

On one 4-core machine with 32 GB, running 3 splits at a time because of memory, expect about 4-5 days.
