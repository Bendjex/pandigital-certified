#!/bin/bash
# certify every cell listed in cells/guesses.txt (const n ord err value time), smallest first
out=${OUT:-cells/cert_results.jsonl}
order="$1"   # e.g. "1:2-8 0:2-6"
for spec in $order; do
  o=${spec%%:*}; r=${spec#*:}; a=${r%-*}; b=${r#*-}
  for n in $(seq $a $b); do
    for c in pi e gamma delta; do
      line=$(grep "^$c $n $o " cells/guesses.txt) || continue
      E=$(echo $line | awk '{print $4}')
      if grep -q "\"n\": $n, \"ord\": $o, \"const\": \"$c\"" $out 2>/dev/null; then continue; fi
      echo "=== $c n=$n ord=$o E=$E  $(date +%T)"
      OMP_NUM_THREADS=${THREADS:-1} python3 -u certify.py $n $o $c $E res_${n}_${o}.txt $out > cells/log_${c}_${n}_${o}.txt 2>&1 || echo "FAILED $c $n $o"
      python3 -c "
import json; r=[json.loads(l) for l in open('$out')][-1]
g=r['optimal']
print('   ', r['const'], r['n'], r['ord'], 'cand', r['candidates'], 'vals', r['values_in_window'], 'ties', r['ties'], 'unres', r['unresolved'], str(r['seconds'])+'s')
for x in g: print('     ', x['error'], x['writings'], x['checked'][0]['formula'], x['checked'][0]['valid_in_range'])
"
    done
  done
done
