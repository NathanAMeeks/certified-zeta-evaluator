# Close-pair finder: deficit-targeted certified refinement [NUM][LEM] — CANDIDATE, UNREFEREED
Date: 2026-10-16. No RH claims: statements are about zeros of ζ in specific finite windows only.
Code: `fftzcp.c` (build: `build.sh`), driver `runall.sh`, outputs `w_*.json`.
Usage: `./fftzcp T delta P a b Lw signfile dbl [point_budget]` (default budget 4e6 extra points, 200 rounds).

## Procedure
1. Coarse grid: certified signs s(t_p) of Z at t_p = T + p·δ (Arb balls; batched Taylor-NUFFT; optional double FFT with rigorous Higham bound).
2. Turing/Trudgian bracket: integers Na ≤ N(T+a) and Nb ≥ N(T+b) from certified sign data in [a−Lw,a] and [b,b+Lw]; max = Nb − Na.
3. Count: found = number of sign changes between consecutive certified points of the merged set in [a,b].
4. While found < max: (a) for each gap between consecutive located zeros, score = step in D(u) = (#changes ≤ u) − (u−a)·log(T/2π)/2π (mean D over 40 changes before minus 40 after; a missed pair is a step of ≈ 2) + 0.25·(expected zeros in the gap); (b) take the top (max − found + 4) gaps not already sampled at the current step; (c) evaluate all of them in one pass with a shared T-independent setup cache (log n, n^{-1/2}, NUFFT bin/offsets), on integer-anchored grids of step dr = δ/10; (d) merge and recount. If no candidate gap is left, dr ← dr/10 (floor δ·1e-4). Stop at found = max (PROVEN), or budget/round cap/dr floor (UNPROVEN_*), or overcount found > max (ERROR), or any two points within 1e-9 with opposite certified signs (UNPROVEN_disagree).
The targeting (steps a–b) is heuristic. It only decides where to evaluate and plays no part in the proof.

## Lemma CP (soundness)
Let S be a finite set of reals in [T+a, T+b], each carrying a sign σ(t) ∈ {±1} with Z(t)σ(t) > 0 rigorously (assumption A1). Order S as t_1<…<t_m, and let k = #{i : σ(t_i) ≠ σ(t_{i+1})}. Let Na, Nb be integers with N(T+a) ≥ Na and N(T+b) ≤ Nb (A2). Suppose ζ has no zeros with ordinate exactly T+a or T+b (A3; implied by Z(T+a), Z(T+b) ≠ 0 being certified, plus A2 applied to the strip). If k = Nb − Na, then every zero of ζ in the strip T+a < Im s < T+b is simple and lies on Re s = 1/2, and their number is exactly k.
Proof. Z is real and continuous on R. Each index i with a sign change gives an open interval (t_i, t_{i+1}) holding a zero of Z of odd multiplicity. These intervals are pairwise disjoint, so Z has at least k distinct real zeros in (T+a, T+b), i.e. at least k critical-line zeros of ζ counted without multiplicity. The total count with multiplicity is N(T+b) − N(T+a) ≤ Nb − Na = k. So the critical-line zeros found already exhaust the count: no zero lies off the line, each interval holds exactly one zero, and it has multiplicity one. ∎
Corollary (never a false PROVEN). The program returns PROVEN only when found = max with no disagreement. Refinement only adds points to S, so k can only grow, and the lemma gives k ≤ N(T+b) − N(T+a) ≤ max whenever A1–A2 hold. So found > max means a broken assumption, and the program reports ERROR. Targeting quality affects run time and whether the stop is PROVEN or UNPROVEN, never correctness.
Implementation obligations behind A1:
- Positions. Refinement grids start at integer offsets o, so T+o is exact in binary64, and t = (T+o) + q·dr is formed exactly in Arb (dr is the exact binary64 value). The double u used for sorting has error ≲ 1e-12. Any two points within 1e-9 are treated as coincident: same sign means the duplicate is skipped, opposite signs abort as a disagreement. So sorting by u gives the true order of all retained points.
- Z(t) = 2Re(e^{iθ(t)} Σ_{n≤N} n^{-1/2} n^{-it}) + C0 + R. The Σ is from the batched Taylor-NUFFT: Arb phases, Taylor truncation enclosed by the explicit bound S0(πr)^L e^{πr}/L!, and the DFT done in Arb (acb_dft) or in binary64 with the Higham (2002, Thm 24.2) radix-2 error bound plus input/output rounding, added as ball radii. θ comes from acb_poly_lgamma_series with a certified cubic Taylor remainder. C0 is evaluated in Arb, and a point is certified only if |cos 2πp| > 1e-3 holds certifiably. R: |R₀| ≤ 0.127 t^{-3/4} (Gabcke 1979; see docs/GABCKE_SOURCES.md).
Behind A2: Turing's method in Trudgian's form, |∫_{t1}^{t2} S(t)dt| ≤ 2.067 + 0.059 log t2 (Trudgian 2011, t ≥ 168π), with the ∫θ integrals enclosed by Arb Taylor bounds, and the gap integrals ∫ taken from certified sign changes only.

## Dependencies and status
- Arb/FLINT 3.6 ball arithmetic: assumed correct. It is mature but not formally verified.
- Higham FFT bound: published theorem. The constants as coded need a referee check, including that W[k] is rounded once.
- Gabcke 1979, Riemann–Siegel remainder: |R_0(t)| ≤ 0.127 t^{-3/4}, t ≥ 200. Second sources: Kuznetsov, arXiv:2312.03261; Gourdon–Sebah 2003; Platt–Trudgian arXiv:1507.01261; Arias de Reyna Math. Comp. 80 (2011). Earlier code used the looser 0.127(t/2π)^{-3/4} (about 3.97× larger); the code now uses 0.127 t^{-3/4}. Citation check only, not a re-proof; see `docs/GABCKE_SOURCES.md`, `docs/GABCKE_CHECK.md`.
- Trudgian 2011 S(t)-integral bound: published. Constants are as quoted in the literature.
- No independent second implementation of fftzcp yet.

## Results (δ=0.01, P=2^20, window [T+60, T+10380], Lw=50, binary64 FFT with rigorous bound; 4 jobs in parallel on 8 cores)
| window T | coarse found / max | refine sub-grids (rounds) | extra pts | found | status | refine s | total s |
|---|---|---|---|---|---|---|---|
| 3e12 (old fail) | 44167/44169 | 6 (1) | 24,576 | 44169 | PROVEN | 11.2 | 26.8 |
| 3e12+10320 (old 3-pair, failed at cap 4) | 44163/44169 | 10 (1) | 40,960 | 44169 | PROVEN | 17.7 | 34.0 |
| 3e12+20000 (old fail) | 44167/44169 | 6 (1) | 24,576 | 44169 | PROVEN | 11.0 | 27.0 |
| 3e12+40000 | 44164/44170 | 10 (1) | 40,960 | 44170 | PROVEN | 18.6 | 34.7 |
| 3e12+51600 (new) | 44166/44168 | 6 (1) | 24,576 | 44168 | PROVEN | 10.9 | 27.3 |
| 3e12+61920 (new) | 44169/44169 | 0 | 0 | 44169 | PROVEN | 0 | 15.0 |
| 3e12+72240 (new) | 44167/44169 | 6 (1) | 24,576 | 44169 | PROVEN | 11.4 | 26.5 |
| 1e13 (new) | 46145/46145 | 0 | 0 | 46145 | PROVEN | 0 | 17.4 |
| 1e13+10320 (new) | 46144/46148 | 8 (1) | 32,768 | 46148 | PROVEN | 28.1 | 45.5 |
| 1e13+20640 (new) | 46140/46146 | 10 (1) | 40,960 | 46146 | PROVEN | 33.9 | 50.5 |
- Disagreements: 0 in all runs. Every deficit window was resolved in one round.
- The old widest-gap loop needed up to 27 sub-grids and 8 rounds on the 3-pair window; this one needs 10 sub-grids and 1 round.
- Contiguous chains, checked by matching endpoint Turing counts (Nb of window k = Na of window k+1):
  - 3e12+51,660 → +82,620: Na 12362403240765 → Nb 12362403373271, so 132,506 zeros. This joins the earlier fftzr/fftzr8 chain at 12362403240765, but that chain's windows k=2..4 were not re-run with fftzcp.
  - 1e13+60 → +31,020: 43124192297371 → 43124192435810, so 138,439 zeros.
- UNPROVEN path test: at 3e12 with budget 0 the program returned found 44167/44169, UNPROVEN_budget, proven=0.
- Cost: clean windows take 15–17 s. A deficit window adds ≈1.8 s per sub-grid at 3e12 and ≈3.4 s at 1e13, so +70–190% here. The shared cache removes only the log/rsqrt part. The per-call O(N) phases and the L-term NUFFT still dominate.

## Remaining gaps
1. Referee check of the Higham constants as coded (empirical check in `docs/DP_CHECK.md`, not a proof); an independent re-implementation. (Gabcke second sources: done, `docs/GABCKE_SOURCES.md`, `docs/GABCKE_CHECK.md`.)
2. Cost: each extra sub-grid costs O(N·L). A true batched multi-interval NUFFT (all suspect gaps in one transform), or a narrower first pick (want = deficit instead of deficit+4), would cut the overhead.
3. Targeting is not tested on pairs closer than dr/10 = 1e-4 or on deficits near the window edges. Those can end UNPROVEN (sound, but needing a rerun), not wrong.
4. Points within 1e-9 of each other are not resolved (they abort as a disagreement). A real pair closer than that would need exact-rational positions.
5. The A3 boundary condition holds because Z(T+a) and Z(T+b) are evaluated as certified points of the coarse grid (60 and 10380 are multiples of δ); the code does not check this explicitly.

## Extension: Lemma TZ zero-free gap test [LEM][NUM] (code `c/closepair/tz/fftztz.c`, runs `results/closepair/tz/w_*.json`; CANDIDATE, UNREFEREED)
Notation. For t in the window, Z(t) = f(t) + C0(t) + R(t), f(t) = 2 Re(e^{iθ(t)} Σ_{n≤N} n^{-1/2} n^{-it}) = 2Σ n^{-1/2} cos φ_n(t), φ_n = θ − t log n, N = ⌊√(t/2π)⌋.

**Lemma TZ.** Let g ∈ C²[a,b], h = b−a, |g''| ≤ M on [a,b], and ℓ the linear interpolant of g at a, b. Then for t ∈ [a,b]: |g(t) − ℓ(t)| ≤ M(t−a)(b−t)/2 ≤ M h²/8. Hence if σ ∈ {±1} and min(σg(a), σg(b)) − M h²/8 > E ≥ 0, then σg(t) > E on [a,b].
Proof. e = g − ℓ satisfies e(a) = e(b) = 0 and e'' = g''. For fixed t ∈ (a,b), set k(s) = e(s) − K(s−a)(s−b) with K chosen so that k(t) = 0. Then k vanishes at a, t, b, so by Rolle twice there is ξ with 0 = k''(ξ) = g''(ξ) − 2K. So e(t) = g''(ξ)(t−a)(t−b)/2, and |e(t)| ≤ M(t−a)(b−t)/2 ≤ Mh²/8. Since ℓ is affine, σℓ(t) ≥ min(σg(a), σg(b)). ∎

**Corollary TZ (certified zero-free gap).** Take consecutive coarse points a, b = a+δ in one window. Suppose (i) N is constant on [T+A, T+B] (checked: ⌊√(t/2π)⌋ over the window ball is one exact integer), (ii) |C0(t)| ≤ C̄ and |R(t)| ≤ ḡ on the window, (iii) certified σf(a), σf(b) ≥ F_a, F_b, and (iv) |f''(a)|, |f''(b)| ≤ D_a, D_b. Let M = max(D_a, D_b) + M4 δ²/8, where M4 = 2Σ_{n≤N} n^{-1/2}(w_n⁴ + 6w_n²e2 + 4w_n e3 + 3e2² + e4), with w_n ≥ |θ' − log n| and e_k ≥ |θ^{(k)}| on the window. If min(F_a, F_b) > Mδ²/8 + C̄ + ḡ, then σZ > 0 on [a,b], so the gap holds no zero of Z, and in particular no missed pair.
Proof. f is analytic, and by Faà di Bruno |d⁴cos φ_n| ≤ φ'⁴ + 6φ'²|φ''| + 4|φ'||φ'''| + 3φ''² + |φ''''|, with φ_n^{(k)} = θ^{(k)} for k ≥ 2. So |f''''| ≤ M4. Lemma TZ applied to g = f'' (two-sided) gives |f''| ≤ M on [a,b]. Lemma TZ applied to g = f with E = C̄ + ḡ gives σf > C̄ + ḡ ≥ |C0 + R|, so σZ > 0. ∎
The formula used is f'' = 2Re(e^{iθ}[(−θ'² + iθ'')S0 + 2θ'S1 − S2]), with S_k = Σ n^{-1/2}(log n)^k n^{-it}.
Implementation:
- S1 and S2 come from two extra batched NUFFT passes, with weights n^{-1/2}(log n)^k and their own truncation bound Σw_n (πr)^L e^{πr}/L!.
- θ', θ'' are taken as Arb balls over each 2-unit block. θ^{(k)}, k ≤ 5, are Arb balls over the whole window (lgamma series of length 6).
- C̄ comes from c0term on the window ball, and ḡ from Gabcke at T+A.
- Final comparisons are on outward-rounded doubles, with a 1e-9 relative margin covering ≤ 6 positive binary64 operations.
**Dependencies:** the same as Lemma CP (Arb, Higham FFT bound, Gabcke 1979 constant, second source in `docs/GABCKE_CHECK.md`). It additionally needs the RS main-term structure Z = f + C0 + R with |R| ≤ Gabcke for real t. No bound on R' or R'' is needed, because only f is differentiated.
**Not done (rigorously): uniqueness inside sign-change gaps.** Proving that Z has exactly one zero in a sign-change gap needs control of Z' = f' + C0' + R', and no certified bound on R' is in hand. Uniqueness is still supplied by Lemma CP (Turing count). Missed pairs are counted only as pairs of certified sign changes located by refinement.
**Role in the procedure:** gaps failing the test are "suspects". Refinement now considers only inter-zero gaps containing a suspect, ranked by severity thr/max(F_a, F_b), and picks ⌈deficit/2⌉+1 per round. Targeting remains heuristic, and soundness rests on Lemma CP. The consistency check `tz_violations` counts certified-zero-free gaps in which refinement later found a sign change, and must be 0.

### Results (same 10 windows, δ=0.01, P=2^20, [T+60,T+10380], 4 jobs/8 cores)
| window | same-sign gaps | zero-free | suspects | missed-pair gaps (all in suspects?) | sub-grids old→new | rounds | refine s old→new | total s old→new |
|---|---|---|---|---|---|---|---|---|
| 3e12 | 987,833 | 986,036 (99.818%) | 1,797 | 1 (yes) | 6→4 | 2 | 11.2→7.4 | 26.8→60.3 |
| 3e12+10320 | 987,837 | 986,028 | 1,809 | 3 (yes) | 10→4 | 1 | 17.7→7.3 | 34.0→59.2 |
| 3e12+20000 | 987,833 | 986,028 | 1,805 | 1 (yes) | 6→2 | 1 | 11.0→4.0 | 27.0→56.8 |
| 3e12+40000 | 987,836 | 986,010 | 1,826 | 3 (yes) | 10→6 | 2 | 18.6→10.5 | 34.7→47.6 |
| 3e12+51600 | 987,834 | 986,062 | 1,772 | 1 (yes) | 6→6 | 3 | 10.9→10.7 | 27.3→64.0 |
| 3e12+61920 | 987,831 | 986,037 | 1,794 | 0 | 0→0 | 0 | 0→0 | 15.0→37.6 |
| 3e12+72240 | 987,833 | 986,082 | 1,751 | 1 (yes) | 6→6 | 3 | 11.4→10.5 | 26.5→48.6 |
| 1e13 | 985,855 | 983,947 | 1,908 | 0 | 0→0 | 0 | 0→0 | 17.4→45.1 |
| 1e13+10320 | 985,856 | 984,037 | 1,819 | 2 (yes) | 8→3 | 1 | 28.1→10.4 | 45.5→54.8 |
| 1e13+20640 | 985,860 | 984,018 | 1,842 | 3 (yes) | 10→6 | 2 | 33.9→19.9 | 50.5→63.7 |
- All 10 windows PROVEN. Refinement disagreements: 0. tz_violations: 0. Na/Nb are identical to the previous runs, so the chains 3e12+51,660→+82,620 (132,506 zeros) and 1e13+60→+31,020 (138,439 zeros) re-verify.
- Zero-free share: 99.81–99.82% of same-sign gaps, about 1,800 suspects per window. The 15 actual missed pairs (one coarse gap each) were all inside suspects, so suspects outnumber real pairs about 1,200:1.
- Cost: refinement fell from 62 to 37 sub-grids and from 143 s to 81 s total. The two extra NUFFT passes cost +22–38 s per window, though, so the whole run is about 1.3–2.5× SLOWER than the current finder. Net: TZ is a soundness and audit gain (missed pairs provably confined to 0.18% of gaps), not a speed gain at δ=0.01.
- Possible speed paths, untested:
  - Use the crude S1/S2-free bound M ≤ 2Σn^{-1/2}((θ'−log n)² + |θ''|) ≈ 32√N. That gives a threshold ≈ 0.33 at 3e12, likely too weak for most gaps.
  - Evaluate S1 and S2 only on blocks near low-|Z| points.
  - Compute S1 and S2 at lower L or precision, since the bound only needs to be crude.
