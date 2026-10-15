#!/bin/sh
# Build fftz against FLINT >= 3 (includes Arb). Set FLINT_ROOT to an install prefix if not system-wide.
cd "$(dirname "$0")"
R=${FLINT_ROOT:-/usr}
gcc -O2 -std=gnu99 -Wall -I$R/include fftz.c -o fftz -L$R/lib -Wl,-rpath,$R/lib -lflint -lmpfr -lgmp -lm && echo built
