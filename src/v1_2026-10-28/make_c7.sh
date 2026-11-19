#!/bin/bash
# build the C7 lab package from the current (fixed) motor A and res_9_0.txt
set -e
cd "$(dirname "$0")"
x86_64-w64-mingw32-gcc -O2 -ffp-contract=off -fopenmp -static -D__USE_MINGW_ANSI_STDIO=1 -o motora.exe motora.c ../c5/coremath/pow.c ../c5/roundeven_shim.c -lm
export WINEDEBUG=-all OMP_NUM_THREADS=1
./motora multi 6 0 res_6_0.txt w_ident6.txt ../l6.txt 2>/dev/null
/usr/lib/wine/wine64 ./motora.exe multi 6 0 res_6_0.txt w_ident6.txt ../w6.txt 2>/dev/null
tr -d '\r' < ../w6.txt | cmp - ../l6.txt && echo "IDENTICAL linux/windows $(wc -l < ../l6.txt) lines"
H=$(sha256sum ../w6.txt | cut -c1-64)
P=../c7pkg; rm -rf $P; mkdir -p $P/c7
cp motora.exe res_6_0.txt res_9_0.txt w_ident6.txt w_c7_9_0.txt $P/c7/
echo $H > $P/c7/esperado_ident6.txt
cp ../c7_templates/*.bat ../c7_templates/LEEME.txt $P/c7/
(cd $P/c7 && sha256sum motora.exe res_6_0.txt res_9_0.txt w_ident6.txt w_c7_9_0.txt esperado_ident6.txt autotest_c7.bat lanzar_c7.bat estado_c7.bat LEEME.txt > HASHES.txt)
rm -f /mnt/user-data/outputs/pandigital/lab_c7.zip
(cd $P && zip -qr9 /mnt/user-data/outputs/pandigital/lab_c7.zip c7)
ls -la /mnt/user-data/outputs/pandigital/lab_c7.zip; sha256sum /mnt/user-data/outputs/pandigital/lab_c7.zip
