#!/bin/bash
# Linux builds (Ubuntu 24.04: gcc 13, libflint-dev = FLINT 3.0.1).  -ffp-contract=off is required; never -ffast-math.
set -e
cd "$(dirname "$0")"
cat motora_head.c motora_body.c motora_tail.c | cmp - motora.c        # motora.c is the concatenation
gcc -O2 -ffp-contract=off -fopenmp -mfma -o motora motora.c coremath/pow.c -lm             # TBL = 5 (needs a CPU with FMA)
gcc -O2 -ffp-contract=off -o motorc motorc.c -lflint -lmpfr -lgmp -lm
gcc -O2 -ffp-contract=off -fopenmp -mfma -DTBL=4 -o motora4 motora.c coremath/pow.c -lm   # TBL = 4, used by build_res4.sh
echo "built motora motorc motora4"
