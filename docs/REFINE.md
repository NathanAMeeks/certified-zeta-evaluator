# Close-pair refinement prototype (3e12, 2^20 window) [NUM]
Date: 2026-10-15. Not pushed to the public repo.
1. Coarse run: `./fftz 3000000000000 0.01 1048576 60 10380 50 sg0.bin dbl` -> all 1,048,576 signs certified, found 44,167, Turing max_possible 44,169, NOT proven (16.3 s).
2. Localize (loc.py, heuristic only, used just to choose where to look): found-count minus (theta/pi) steps down by 2.00 at u~5477.5; the widest gap there is 5476.855-5477.865 (4.3x mean spacing).
3. Refine: `./fftz 3000000005476 0.001 2048 ... sg2.bin dbl` (certified Arb signs on a 1e-3 grid over 2 units) -> close pair at u in (5477.492,5477.493) and (5477.499,5477.500), separation ~0.007 (~3% of mean spacing).
4. Merge (merge2.py): union of certified points from both grids; 205 coincident points agree in sign (0 disagreements). Found 44,169 = max_possible -> PROVEN, all 44,169 zeros in [3e12+60, 3e12+10380].
Rigor: localization is heuristic but only picks where to add points; every added sign is a certified Arb ball, and the proof is still certified sign changes = Turing/Trudgian count. Refinement was interior, so Turing end-window bounds are unchanged.
Cost: extra points 0.2% of the window, but each fftz call has ~2.4 s fixed setup (O(N) phases + NUFFT), i.e. +15% for one flagged gap per window. Next: refinement in-process reusing the window's setup (point cost then ~0), and an automated loop (flag -> refine -> recheck) for multiple deficits / 2^21 windows.

## Window 2 (T=3e12+20000, same 2^20 / 0.01 grid) [NUM]
Date: 2026-10-16. Automation advance of "test more windows".
1. Coarse: `./fftz 3000000020000 0.01 1048576 60 10380 50 sg_w2.bin dbl` -> all signs certified, found 44,167, max_possible 44,169, NOT proven (15.0 s).
2. Localize (same heuristic): step of 2.00 at u~2717.5.
3. Refine: `./fftz 3000000022716 0.001 4096 ... sg_w2_r.bin dbl` over 4 units -> close pair at u in (2718.640,2718.641) and (2718.649,2718.650), sep ~0.009.
4. Merge: 400 coincident points, 0 disagreements. Found 44,169 = max_possible -> PROVEN.
Confirms the first-window pattern: one close pair per deficit-2 window; separate refine pass is enough to prove, but still not in-process.

## In-process refinement, one pass (inproc/fftzr.c) [NUM]
Date: 2026-10-16. Not pushed to the public repo. fftz.c's grid code moved into zgrid() (per-call FFT buffers, no statics); main runs the coarse grid, Turing/Trudgian bounds, then, if found < max_possible, refines the widest interior gaps (>1.5x mean spacing) on a delta/10 certified Arb grid inside the same process, merges certified points, and recounts. Up to 4 rounds, (deficit/2+4) gaps per round. Heuristic only picks where to look; proof = certified sign changes in [a,b] equal the Turing max, with 0 disagreements at coincident grid points. Coarse signs are byte-identical to fftz for windows 1 and 2.
| window (2^20 pts, delta 0.01, [60,10380]) | coarse found / max | refine windows | extra pts | found | proven | total s |
|---|---|---|---|---|---|---|
| 3e12 | 44,167 / 44,169 | 3 | 6,144 | 44,169 | yes | 21.5 |
| 3e12+20000 | 44,167 / 44,169 | 10 | 18,432 | 44,169 | yes | 37.8 |
| 3e12+40000 (new) | 44,164 / 44,170 | 19 | 33,792 | 44,170 | yes | 59.2 |
| 3e12+60000 (new) | 44,169 / 44,169 | 0 | 0 | 44,169 | yes | 14.5 |
| 1e13 (new) | 46,145 / 46,145 | 0 | 0 | 46,145 | yes | 16.6 |
Total 222,822 zeros proven, all on the line, each window in one command. Extra points <=3.3% of the window, but each sub-window pays ~2.3 s of O(N) setup, so a 3-pair window costs ~4x a clean one. Gap-width ranking picks wrong gaps first (w2: 10 sub-windows for 1 pair). Next: one batched sub-grid call for all candidate gaps (or deficit-step localization) to cut refine time to roughly one setup; windows are still not contiguous, and this is not an independent implementation.

## Contiguous stretch 3e12+60 to 3e12+51660 (5 chained windows) [NUM]
Date: 2026-10-16 (afternoon). Not pushed to the public repo. Windows at T=3e12+10320k, k=0..4, each [60,10380], so each window's b is the next one's a; Turing/Trudgian endpoint bounds chain exactly (Nb_upper of k = Na_lower of k+1). Runs in inproc/contig/.
| k | coarse found / max | refine windows | found | proven | total s |
|---|---|---|---|---|---|
| 0 (run_w1) | 44,167 / 44,169 | 3 | 44,169 | yes | 21.5 |
| 1 | 44,163 / 44,169 | 22 (4 rounds) | 44,167 | NO with fftzr (4-round cap) | 63.8 |
| 1 (fftzr8, 8 rounds) | 44,163 / 44,169 | 27 | 44,169 | yes | 71.2 |
| 2 | 44,167 / 44,169 | 10 | 44,169 | yes | 37.8 |
| 3 | 44,168 / 44,170 | 5 | 44,170 | yes | 26.9 |
| 4 | 44,166 / 44,170 | 18 | 44,170 | yes | 55.3 |
Total 220,847 zeros = Nb_upper(k=4) - Na_lower(k=0) = 12362403240765 - 12362403019918, all on the critical line, 0 disagreements. fftzr8 = fftzr with rounds<8 (only change); coarse signs byte-identical to fftzr. Finding: the 4-round, widest-gap-first loop is not robust (k=1 had 3 close pairs and missed one). Fix before publishing: rank by the deficit-step localizer (loc.py step at u~9372 for k=1) or loop until proven with a point budget, and batch sub-grids to cut the ~2.3 s per sub-window setup.
