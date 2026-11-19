#!/bin/bash
# Windows builds (static), cross-compiled on Linux with mingw-w64 gcc 13 and the prebuilt static GMP/MPFR/FLINT
# of winlibs/ (opt_sesion_parte2.zip): no external code is downloaded or compiled.
#   usage: build_windows.sh PATH/TO/mingw64_static_gmp-6.3.0_mpfr-4.2.1_flint-3.0.1.tar.xz
set -e
cd "$(dirname "$0")"
T=winlibs_tmp; rm -rf $T; mkdir $T; tar -xJf "$1" -C $T; W=$T/prefix     # fixed path: reproducible binaries
x86_64-w64-mingw32-gcc-posix -O2 -ffp-contract=off -fopenmp -mfma -static -s -Wl,--no-insert-timestamp -D__USE_MINGW_ANSI_STDIO=1 \
    -o motora.exe motora.c coremath/pow.c roundeven_shim.c -lm
x86_64-w64-mingw32-gcc -O2 -ffp-contract=off -static -s -Wl,--no-insert-timestamp -D__USE_MINGW_ANSI_STDIO=1 -I$W/include -o motorc.exe motorc.c \
    -L$W/lib -lflint -lmpfr -lgmp -lpthread -lm
rm -rf $T
sha256sum motora.exe motorc.exe
