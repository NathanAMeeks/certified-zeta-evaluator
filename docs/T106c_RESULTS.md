# t106c results (Nathan A. Meeks & Grok; candidate, unrefereed) — 2026-11-02 15:00 MST
All runs scratch-only, nice 19, 8 vCPU Xeon @ 2.4 GHz, sharing the box with the live loop (load avg 7-18). Wall times are load-noisy; CPU times are the steadier comparison. Raw: bench/res.txt, bench/sweep_res.txt, bench/winlen_curve.csv.
zeros/ms/GHz/core: wall basis = net_zeros/(wall_ms*2.4*8); CPU basis = net_zeros/((user+sys)_ms*2.4).

## Per-step (start idx 74898)
| step | N | load | wall s | cpu s | s/win wall | z/ms/GHz/core wall | cpu |
|---|---|---|---|---|---|---|---|
| t103 | 100 | 15.7-18.3 | 85.2 | 122.0 | 0.852 | 2.7 | 15.0 |
| s0 t106b | 100 | 18.3-15.5 | 35.8 | 115.7 | 0.358 | 6.4 | 15.8 |
| s1 XZLEVEL=1 | 100 | 15.5-14.2 | 60.6 | 58.0 | 0.606 | 3.8 | 31.6 |
| s3 GATE_PROCS=6/7/8 | 100 | ~15 | 46.4/84.6/54.9 | 56.0/57.9/56.2 | — | noise-dominated | ~32 |
| s4 V9DIRECT=1 | 100 | 12.3-13.6 | 134.9 (26.0 rc5 crash run before _R.args fix) | 74.9 | — | load spike | 24.5 |
| final400 (xz1,direct,gp7; GRS_TO_MIN=5) | 400 | 13.6-14.4 | 159.5 | 287.2 | 0.399 | 4.8 | 25.5 |
| **final400_to30** (same + GRS_TO_MIN=30) | 400 | 7.8-7.4 | **46.1** | 269.2 | **0.115** | **19.9** | **27.2** |
| **S1 PC_FAST=1** (+ above) | 400 | 8.0-7.8 | **42.5** | 276.1 | **0.106** | **21.6** | 26.6 |
Net zeros: 4,399,837 (N=100), 17,599,347 (N=400). HYB_TMESH=1 was already on in live and all benches (no change). xz -1 halves CPU (biggest real win). GATE_PROCS 6-8 indistinguishable under load; 7 kept. V9DIRECT saves ~0.1 cpu-s/win (wrapper fork/grep/awk/sed). S1 gives ~8% wall, ~0 CPU — marginal; keep off by default unless Audit wants it.
t103 vs final400_to30 were at different loads (t103 at ~17); live t103 reference is 0.52 s/win, so t106c is ~4.5x on wall vs live t103, ~1.8x on CPU vs the bench t103.

## Byte identity vs t103 (t106b/bench/t103_v2, 400 windows)
S1_pcfast_to30 and final400_to30: verdicts 400/400 (cmp_verdicts.py, 19 IND-gated), ledger cols 1-12 sha 43cc5aed = t103, uncompressed brackets (concatenated, sorted) sha 701ccaa5 = t103, 400 PROVEN, sum zeros_turing 17,667,837 both.

## WIN_LEN sweep (bench/winlen_curve.csv; same t-range ≈ 4.112e6 units, N=(4112000+L-121)/(L-120))
| WL | N | wall s | cpu s | ms/zero cpu | retries | fails | GRS timeouts | maxrss MB | redo est cpu s |
|---|---|---|---|---|---|---|---|---|---|
| 10400 | 400 | 143.8 | 298.9 | 0.0170 | 0 | 0 | 6 | 434 | 12.8 |
| 15600 | 266 | 347.3 | 517.9 | 0.0294 | 0 | 0 | 4 | 421 | 14.2 |
| 20800 | 199 | 642.3 | 612.6 | 0.0348 | 0 | 0 | 9 | 481 | 23.9 |
| 31000 | 134 | 204.0 | 380.6 | 0.0215 | 0 | 0 | 2 | 465 | 20.4 |
| 41600 | 100 | 535.3 | 604.3 | 0.0340 | 0 | 0 | 8 | 624 | 46.2 |
| 62400 | 67 | 539.4 | 672.9 | 0.0377 | 0 | 0 | 8 | 929 | 79.0 |
All lengths: every window PROVEN, ledger N chain continuous, first Na = t103, and N(T) equals t103 at every shared endpoint (400/400 at 10400; 0-2 shared endpoints at other lengths since grids differ). No retries or fails at any length, so the retry term never enters; cost per zero rises with length (bigger P/FFT and G-RS point count per window) and memory rises above 31000 (929 MB at 62400).
**Best: WL=10400 (current default)** — lowest CPU/zero, lowest redo. 31000 is second, but no length beats 10400; no change recommended.

### Kill/resume walls (1229 s at 10400, 935 s at 15600): load artifact, not real redo
Resume only re-verifies committed rows (7.6 s at 10400 for 40 rows) and redoes ≤7 in-flight windows (~13 cpu-s). The long walls come from GRS_TIMEOUT under load: with GRS_TO_MIN=5 the fast G-RS tiers timed out (8 per resume), were dropped, and points went to Arb (CPU 787 s for 355 windows vs 299 s for 400 in the full run). Same cause made final400 159.5 s; with GRS_TO_MIN=30 the identical run is 46.1 s. Recommend live GRS_TO_MIN=30 (bound still caps a hung child; verdicts unchanged). Real redo cost after a kill = verify + in-flight windows, column redo_cost_est_cpu_s.

## Deferred (design work, not staged)
- S2 local repair: gates.py close_pair is already local (±20 units around the bracket); the remaining part, local repair inside the v9 retry ladder, needs a certified splice of a partial re-run into the window's Turing bracket — a design + proof task, and with 0 retries at all lengths the payoff is ~0 now.
- S3 fftzi v10 Z'' skip: needs new C plus an audited Z'' bound from the interpolation and a proof that a skipped zero always forces a Turing mismatch; not safely stageable in a bench session.
- S5 persistent fftzi: fftzi has function-level statics; a pipe-fed persistent mode needs a statics/reset audit first so no state leaks between windows. V9DIRECT already removed the shell-wrapper part of the per-window overhead.
So the 16-combo test reduced to S1 on/off (both identical to t103).
