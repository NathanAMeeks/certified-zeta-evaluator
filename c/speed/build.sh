#!/bin/sh
# Build against FLINT >= 3 (includes Arb). Set FLINT_ROOT to an install prefix if not system-wide.
cd "$(dirname "$0")"
R=${FLINT_ROOT:-/usr}
I="-I$R/include -I$R/include/x86_64-linux-gnu"; L="-L$R/lib -L$R/lib/x86_64-linux-gnu -Wl,-rpath,$R/lib:$R/lib/x86_64-linux-gnu"
gcc -O2 -std=gnu99 $I fftzsp.c -o fftzsp $L -lflint -lmpfr -lgmp -lm || exit 1
gcc -O2 -std=gnu99 $I dpcheck.c -o dpcheck $L -lflint -lmpfr -lgmp -lm || exit 1
echo built
