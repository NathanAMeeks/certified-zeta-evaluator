# v9p(h) vs T106b: where the 8x goes (Audit, 2026-11-02 ~11:55 MST, read-only)
Data: live ledger batch idx 107318-108617 (11:09:5x-11:14:47 MST, 1300 windows, 57.2M zeros, ~282 s wall, 8 cores = 2256 core-s);
bench40.txt = 40 fresh v9p windows at live height T=3.00117e12 (core 7, private ahcache copy, box load ~13); xz timing on a live bracket copy.
Units: Nathan's zeros/ms/GHz/core = 0.4165 / (core-ms per zero)  (0.0050 -> 83.3).

## (a) v9p(h) real speed, certification
- 40 windows: 39/40 PROVEN on attempt 1, 1 not (needs ladder). Median 0.314 s/window = 0.0071 core-ms/zero (58.7 units); mean incl. the miss 0.0087 (48 units). 83.3 holds only on an idle box (Oct bench 0.33 s/win incl. spawn).
- Certified? Signs yes (every bracket endpoint is a certified Arb sign, Lemma CP; Turing count). But NO ledger, NO independent referee (G-Turing/G-RS), and it cannot finish the ~6% close-pair windows by itself. Its 0.0050 excludes the windows that cost the most.

## (b) T106b budget for that batch (2256 core-s)
| # | component | core-s | share | fix | est. saving | keeps rigor |
|---|---|---|---|---|---|---|
| 1 | xz -6 bracket compression (XZPAR=6, 6 procs at 100% seen live) | ~740 (0.57 s/win measured) | 33% | xz -1 (0.038 s/win, +28% size) or zstd; brackets_sha is on the raw file, so no ledger change | ~690 s | y |
| 2 | gated windows (80 = 6.2%: close_pair/pair_cert + ind_gate referee G-Turing+G-RS), mean 10.7 s, max 63 s | 858 | 38% | turn on audited HYB_VC+HYB_TMESH (referee 2.07x, refine/hybrid); refine pair_cert to cover more slow cases | ~430 s | y (VC breaks G-Turing/G-RS independence: use TMESH only, or VC within one gate) |
| 3 | v9p evaluator, clean windows (0.358 s live vs 0.22 idle: contention from 14 gate procs + 6 xz on 8 cores) | 466 | 21% | after #1, GATE_PROCS 14->7 | ~136 s | y |
| 4 | per-window driver overhead (wall - v9: Python, bash v9run.sh spawning grep/awk/sed, ZDUMP, sha, atomic write+fsync) | 115 | 5% | call fftzi_v9 directly / persistent worker; keep sha+fsync | ~85 s | y |
| 5 | rest: orchestrator Python (~25% CPU), in-order commit stalls (largest 6 commit gaps = 56 s wall), checkpoints | ~77 | 3% | none needed (reorder buffer already exists) | small | - |
Useful evaluator work = 466/2256 = 21% of the 8 cores. Safety checks (hash pin, sha, fsync, lock, ledger guards) are inside #4/#5: <2 ms each, not a bottleneck.

## (c) In T106b, not in v9p
- Required for the certified ledger/referee: ledger+brackets_sha+fsync/atomic write (cheap, keep); close-pair cert + ind_gate referee on gated windows (required, can be made cheaper: #2); evaluator hash pin (once per process, cheap); Turing check (inside v9p already).
- Required but cheaper possible: compression (#1, xz -6 is pure size optimisation); v9run.sh bash ladder (#4); gate pool size (#3).
- Unnecessary/duplicate: xz -6 vs -1 level; oversubscribed pool (14 procs/8 cores); grep/awk/sed spawns per window.

## (d) Path and projection (units)
now 10.6 -> #1 15.2 -> +#2 20.9 -> +#3 23.8 -> +#4 ~26. Ceiling with the referee kept ~30-35; full 83.3 is not reachable while 6% of windows need certified repair, and 83.3 itself is an idle-box number (measured 48-59 under current load).
Biggest single lever is #1 (a config change, xz level), then #2.

## Virtual-core idea status
Tested 2026-11-01 as "One giant core" (meeks_test/refine/orchestra/RESULTS.md §A): WIN for latency (heavy window ~2x faster wall), KILL for throughput (+2-5% Turing CPU, +28% G-RS CPU); sim (§C) showed no throughput gain. Related "M windows in one fftzi run" (refine/hybrid RESULTS.md, Grid amortization): KILL (0.56x at M=8). What shipped instead was t102 (WHITEPAPER G.5.t102): core 0 orchestrator + shared pinned gate pool across all 8 cores. t103 "per-core worker processes" option was never built. Parked, not lost; not worth reviving for throughput.
