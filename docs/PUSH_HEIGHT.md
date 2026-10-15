# Pushing verified height beyond Platt–Trudgian 3e12 (fftarb.py + Turing/Trudgian) [NUM windows / HEUR roadmap]
## (1) Measured windows: fully Arb, 131,072 grid points, one Python core
| window | pts/unit | pts/zero | Arb time | certified | zeros proven (Turing) | N(start) proven | time/zero (grid) |
|---|---|---|---|---|---|---|---|
| 3e12+[60,460] | 250 | 58 | 19.8 s | 100% | 1,712 | 12,362,403,019,918 | 8.8 ms |
| 3e12+[300,1800] | 60 | 14 | 21.9 s | 100% | 6,420 | 12,362,403,020,945 | 2.3 ms |
| 1e13+[60,460] | 250 | 56 | 26.2 s | 100% | 1,789 | 43,124,192,297,371 | 11 ms |
In all three runs, float-vs-Arb referee: 0 violations, with slack ≥ 1.2e4×. Gabcke + Trudgian constants are cited. These are the first windows above 3e12, so they are tiny (~10^-9 of the range). They do not extend the verified height in the sense of Platt–Trudgian (all zeros from 0 to T).
Points per zero: 14 sufficed at 3e12 (60/unit). 20/unit failed earlier at 1e6, so ~10–15 pts/zero is the working floor.
Cost model, fitted from 1e6, 2.68e11, 3e12, 1e13: time ≈ 5.6e-7·N·L + 2.55e-6·M·L s, where N=√(t/2π), L≈19+log10(t/1e13) and M=2P. With P=2^20 per window, the cost is ≈1.7–2.1 ms/zero in Python. It is nearly flat in t, as Odlyzko–Schönhage predicts.
## (2) Roadmap: verify ALL zeros from 3e12 to T (HEUR extrapolation, $0.03/core-h)
| T | zeros to check | core-h (Python) | $ (Python) | $ (C/FLINT, guess ÷20) |
|---|---|---|---|---|
| 1e13 | 3.1e13 | 1.5e7 | $434k | ~$22k |
| 3e13 | 1.2e14 | 6.0e7 | $1.8M | ~$89k |
| 1e14 | 4.6e14 | 2.3e8 | $7.0M | ~$350k |
| 1e15 | 5.0e15 | 2.9e9 | $87M | ~$4.4M |
## Limits / what needs C+FLINT
- Python-Arb objects cost ~100–200 B each, so P ≳ 2^21 points per window exhausts RAM. Production needs P ≫ N to amortise the N·L term (N=1.3e7 at 1e15), so it needs flat C arrays (FLINT acb_dft / arb vectors, or double-double FFT with a proven error bound).
- The ÷20 C speedup is a guess, not measured. A C port and benchmark is the next real step.
- Precision: 96 bits is ample to 1e15 (radii ~3e-10 vs |Z|≫). T·ln n mod 2π needs ~log2(T)+60 bits; grow prec ~4 bits/decade.
- Parallel: windows are independent, so the work is embarrassingly parallel. The bottleneck is the scatter (N·L) and memory bandwidth.
- Rigor gaps for a publishable extension: an independent second implementation; Gabcke and Trudgian constants rechecked; a complete, gapless chain of Turing windows with every window proven; and a check for multiple/close zeros (Turing fails on them, which is detectable, but they then need a local refinement).
Verdict: our method proves zero counts in windows just above 3e12 and at 1e13 (6,420 and 1,789 zeros, all on the line). Extending the full verified height is a $20k-scale (C, guessed) to $400k-scale (Python) job for 1e13. That would move U from ~0.163 to ~0.157 (HEUR). Doing it needs a C/FLINT port first.
