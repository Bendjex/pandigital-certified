#!/bin/bash
# Verification of the C7 output (n = 9, any order): motor C (Lemma 15) at 1024 bits with 128-bit tables
# upgraded on use, two worker processes sharing the tables; then 2048 bits on the remainder; then motor B.
# usage: c7_verify.sh CANDIDATES [constants]
#   CANDIDATES: c7_cand.txt, or a .zip / .gz holding it (read as a stream, never stored uncompressed)
#   lines "W0 C ..." pi, "W1 C ..." gamma, "W2 C ..." delta; constants default "pi gamma delta"
set -e
cd /tmp/claude-0/-home-claude/38c805de-6fea-5116-a7a4-e9000a897543/scratchpad/cert
IN=$1; WHICH=${2:-"pi gamma delta"}; D=c7work; mkdir -p $D/tmp
names=(pi gamma delta); E=(2.82e-14 1.15e-15 1.39e-14)
stream() { case "$IN" in *.zip) unzip -p "$IN" ;; *.gz) zcat "$IN" ;; *) cat "$IN" ;; esac; }
for k in 0 1 2; do
  c=${names[$k]}
  case " $WHICH " in *" $c "*) ;; *) continue ;; esac
  sed -n "$((k+1))p" w_c7_9_0.txt > $D/w_$c.txt
  stream | grep "^W$k " | sed "s/^W$k /W0 /" | sort -u -S 1500M -T $D/tmp --compress-program=gzip > $D/u_$c.txt
  echo "[$(date +%T)] $c: $(wc -l < $D/u_$c.txt) unique candidates"
  MOTORC_TPREC=128 MOTORC_JOBS=2 ./motorc_new 9 0 $D/w_$c.txt $D/u_$c.txt $D/k1_$c.txt 1024 2> $D/mc1_$c.log
  echo "[$(date +%T)] $c: after 1024 bits $(wc -l < $D/k1_$c.txt)"
  MOTORC_TPREC=128 ./motorc_new 9 0 $D/w_$c.txt $D/k1_$c.txt $D/k2_$c.txt 2048 2> $D/mc2_$c.log
  echo "[$(date +%T)] $c: after 2048 bits $(wc -l < $D/k2_$c.txt)"
  cp $D/k2_$c.txt c_${c}_9_0.txt
  OMP_NUM_THREADS=1 REUSE_CAND=1 FULLMAX=4 python3 -u certify.py 9 0 $c ${E[$k]} res_9_0.txt cells/cert_results.jsonl > cells/logB_${c}_9_0.txt 2>&1
  echo "[$(date +%T)] $c: motor B done: $(grep -m1 '"error"' cells/logB_${c}_9_0.txt)"
done
