#!/bin/bash
# usage: run_cells.sh FILTER   (lines of cells.txt matching FILTER, run in order, skipping done ones)
cd "$(dirname "$0")"
OUT=out/c11.jsonl
grep -E "$1" cells.txt | while read n o c E; do
  if grep -q "\"n\": $n, \"ord\": $o, \"const\": \"$c\"" $OUT 2>/dev/null; then continue; fi
  res=../src/res_${n}_${o}.txt
  [ -e $res ] || { echo "no $res"; continue; }
  echo "=== $c $n $o $(date +%T)"
  OMP_NUM_THREADS=${THREADS:-1} WORK=run_${n}_${o} python3 -I c11_certify.py $n $o $c $E $res $OUT > out/log_${c}_${n}_${o}.txt 2>&1 || echo "FAILED $c $n $o"
  python3 -c "
import json; r=[json.loads(l) for l in open('$OUT')][-1]
for x in r['optimal_no_minus']: print('   ', r['const'], r['n'], r['ord'], x['error'][:30], x['formula_no_minus'], 'um', x['has_unary_minus'], x['digits_ok'], x['valid_in_range'], x['sign_ok'], 'rounds', len(r['rounds']), 'unres', r['unresolved'], 'warn', r['warn_parity'], r['seconds'])
"
done
