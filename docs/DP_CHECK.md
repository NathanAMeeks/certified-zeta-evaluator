# Referee note: double-precision error bounds in c/speed/fftzsp.c

Date: 2026-10-17. Harness: `c/speed/dpcheck.c`, built by `c/speed/build.sh`, run as `dpcheck T delta P` (fftzsp.c with main replaced; runs the all-Arb path and the double path on the same grid and compares F point by point).

Test: a violation means |F_dbl - mid(F_arb)| + rad(F_arb) > EPH for the real or imaginary part of any grid point.

| T | delta | P | N | EPH (bound) | max observed diff | diff / bound | violations | certified signs Arb / dbl |
|---|---|---|---|---|---|---|---|---|
| 3e12 | 0.02 | 2^18 | 690988 | 3.30e-8 | 2.94e-11 | 8.9e-4 | 0 | 262144 / 262144 |
| 1e13 | 0.02 | 2^18 | 1261566 | 4.69e-8 | 3.88e-11 | 8.3e-4 | 0 | 262144 / 262144 |
| 3e12 | 0.01 | 2^20 | 690988 | 6.54e-8 | 8.71e-11 | 1.3e-3 | 0 | 1048576 / 1048575 |
| 1e13 | 0.02 | 2^20 | 1261566 | 1.84e-7 | 1.88e-10 | 1.0e-3 | 0 | 1048576 / 1048576 |
| 3000000010320 | 0.001 | 4096 | 690988 | 2.93e-7 | 2.52e-14 | 8.6e-8 | 0 | 4096 / 4096 |
| 1e13 | 0.0001 | 2048 | 1261566 | 1.38e-6 | 5.60e-14 | 4.1e-8 | 0 | 2048 / 2048 |

Result: zero violations over about 2.37M grid points (both components). Observed error is at most about 0.13% of the hand bound, so the bound is conservative by roughly 750x or more. This is an empirical consistency check, not a proof of the bounds.

Code-read findings:
1. q = hf*X (line 88) is exact only because P (hence half = P/2) is a power of two. The 12u phase budget does not cover rounding here otherwise. All current callers use powers of two (argv P = 2^20, refinement PRw = 256*2^k), but nothing enforces it. Suggest: abort unless P is a power of two (and RPR likewise).
2. X~ = frac(delta*a_hi) drops delta*a_lo; together with the product rounding this is <= 2u*delta*a_max, inside the stated 2.5u*delta*a_max.
3. Phase reduction rounding count (fma two-product exact, frac steps exact, three rounded adds with |values| < 2.01) gives about 4u, inside 12u.
