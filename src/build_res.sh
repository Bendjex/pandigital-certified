#!/bin/bash
# staged table construction with motor B resolutions:  build_res.sh n ord resfile
n=$1; o=$2; res=$3; touch $res
pend=pend_${n}_${o}.tmp; log=stage_${n}_${o}.log
for k in 2 3 4 5; do
  [ $k -gt $n ] && break
  for pass in 1 2 3; do
    ./motora stage $n $o $k $res $pend 2> $log || { cat $log; exit 1; }
    np=$(grep -c " | " $pend)
    echo "size $k pass $pass: pending $np  ($(grep "tables size $k" $log | sed 's/.*sets, //'))"
    [ "$np" = "0" ] && break
    python3 resolve.py $n $o $pend $res || exit 1
  done
done
