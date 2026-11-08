#!/bin/sh
# Build against FLINT >= 3 (includes Arb). Set FLINT_ROOT to an install prefix if not system-wide.
cd "$(dirname "$0")"
R=${FLINT_ROOT:-/workspace/rh-research/tools/flintc/root/usr}
I="-I$R/include -I$R/include/x86_64-linux-gnu"; L="-L$R/lib -L$R/lib/x86_64-linux-gnu -Wl,-rpath,$R/lib:$R/lib/x86_64-linux-gnu"
for v in v9; do gcc -O2 -march=native -ffp-contract=off $I fftzi_$v.c -o fftzi_$v $L -lflint -lmpfr -lgmp -lm || exit 1; done
echo built
