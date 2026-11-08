# Reproducible build and bounded tests
Nathan A. Meeks & Grok. Candidate/unrefereed; not RH.

Requirements: Linux/POSIX, C compiler (GCC tested), FLINT >=3 with Arb headers, GMP, MPFR, libm; Python 3 with `python-flint` (0.9.0 tested), `mpmath`, `numpy`. Standard `xz`, `tar`, `sha256sum`, bash and util-linux/taskset are needed for the historical wrappers and data. Tested C FLINT was 3.6.0, libflint.so.24 (system installation); Python's python-flint supplies a different FLINT build. Version details are recorded in the preparation report. Obtain dependencies from their upstream distributions; no vendored binaries or secret credentials are required.

Set FLINT_ROOT to the installed prefix when calling the original build script:

```sh
cd proof-runtime/v9p
env FLINT_ROOT=/path/to/flint-prefix sh build.sh
```

The original build defaults to a historical `/workspace` prefix; override it. It uses `gcc -O2 -march=native -ffp-contract=off` and no fast-math. For an explicit portable safety-flags build (substitute real include/lib directories):

```sh
gcc -O2 -ffp-contract=off -fno-fast-math -I/PREFIX/include   fftzi_v9.c -L/PREFIX/lib -Wl,-rpath,/PREFIX/lib   -lflint -lmpfr -lgmp -lm -o fftzi_v9
```

Use the same flags/link libraries to compile `t_psi.c` to `t_psi`, and `tests/t_nstep_straddle.c` to `tests/t_nstep_straddle`. Run:

```sh
./tests/t_nstep_straddle
python3 val_psi.py
```

`val_psi.py` is an mpmath-based regression comparison, not an independent rigorous proof; examine its `violations` count because its exit code alone does not report violations. N-step regression checks refusal of a straddling ball and side-index acceptance; it is not a complete continuity theorem test.

One bounded single-core evaluator smoke test, writing only to your current scratch directory (never reuse a live cache/output):

```sh
env AHCACHE="$PWD/smoke-cache.bin" ZDUMP="$PWD/smoke-brackets.txt"   ./fftzi_v9 2999999999920 0.24 131072 80 31040 60 /dev/null dbl
```

Expected: `status=PROVEN`, `found=max_possible=132508`, Na=12362403019661, Nb=12362403152169. On a cache miss the evaluator creates its coefficient cache; no cache is a mandatory publication dependency. Preparation reproduced these values and rebuilt a binary byte-identical to the live evaluator. Run on one pinned allowed CPU via taskset if timing it; this smoke run is not a new benchmark.

The exact historical `v9/v9run.sh` retry wrapper is retained, including old environment/path conventions, for provenance. The direct command above avoids production paths and shell environment changes. Never compile with `-ffast-math`; preserve `-ffp-contract=off`.

Preparation package versions: mpmath 1.3.0, numpy 2.5.4; Python 3.13.5; GCC 14.2.0; python-flint 0.9.0.
