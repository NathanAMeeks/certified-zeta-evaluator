# SPEED.md — speeding up the certified Z evaluator (candidate, unrefereed; no RH claims)
Code: c/speed/fftzsp.c (from c/closepair/fftzcp.c; run with trailing `dbl`). Data: results/speed/t_*.json, s_*.json. Box CPU: Intel Xeon (model 207), 2.4 GHz, 8 vCPU; runs were 2–4 at a time.

## Validation
Chains rerun as single windows: 3e12+51,660→+82,620 and 1e13+60→+31,020, P chosen per δ.
- δ=0.02 (≈5 evals/zero) and δ=0.03 (≈7.5): both PROVEN, 0 refinement disagreements, Turing endpoints identical to before
  (12362403240765→12362403373271 = 132,506; 43124192297371→43124192435810 = 138,439). Also PROVEN at δ=0.0075, 0.015, 0.04 (3e12).
- 1e13 window m1, δ=0.01: final code's sign file is byte-identical to fftzcp's (sg_m1.bin).
- δ≥0.06 FAILS (UNPROVEN): misses grow fast, and the Turing zones get undersampled, so Nb_upper loosens (…373273, …435814) and
  refinement (which only targets [A,B]) can't make up the difference. Nothing false was certified. Per-window chains were not redone.
- Lemma TZ (tz/fftztz.c) was not ported: it was already slower. Idea 1 (Z'' only near small |Z|) isn't needed in this path.

## Changes and timings (ms/zero/core)
Window m1 (1e13, δ=0.01, 46,145 zeros): baseline 17.44 s, 0.378 ms/zero
 (3) double phases (double-double reduction plus Arb-certified 1024-entry cis table with a degree-7 Taylor series), double scatter: 12.1 s, 0.263
 (3) double per-point F = Σ_k fac_k·Ĝ_k: 6.67 s, 0.145
 (4) per-point Arb at 64 bits after reducing θ0 mod 2π per block, plus Gabcke bound once per block (taken at the block's lower t): 5.65 s, 0.122
Refinement (3e12 chain at δ=0.03, 182 refinement calls): 304 s → 81 (double phase) → 41 (double-double instead of Arb) → 33 (no unused Arb vectors)
 → 30 (M=16P, so Lt 18→11) → 17.7 (refinement windows 256 points, not 4096) → 4.8 s (fused one-pass scatter over all k, cached twiddles). 26 ms per refinement, not 1.7 s.
 (5) True single-transform batching wasn't needed once each refinement call cost about 26 ms.
 (2) Density sweep, total seconds per 3-window chain (3e12 / 1e13): δ=.0075 30.7/33.8; .015 22.8/23.4; .02 11.9/14.8; .03 11.9/18.8; .04 19.6/–; .06 and .12 FAIL.
  Best: δ≈0.02. Getting to Gourdon's ~1.2 evals/zero would need rigorous band-limited interpolation (Odlyzko–Schönhage/Platt), which isn't done.

**Before → after:** 3e12 chain 68.8 s → 11.9 s (0.519 → 0.090 ms/zero, 5.8×); 1e13 chain 113.5 s → 14.8 s (0.820 → 0.107, 7.7×).

## Rigor of the double paths (constants coded in fftzsp.c; need a referee)
u=2^-53. All bounds are absolute, added to F as Arb error (EPH, ≈1–2e-7 at P=2^20–2^22, which is far above Gabcke; it is dominated by P·δX):
- Phase ψ_n=frac(T a_n + (P/2)X_n), a_n=log n/2π, held as double-double (Arb-derived, radius ARAD). T must be an integer double <9e15.
  fma two-product, exact frac steps: δψ ≤ T·ARAD + 12u. X̃_n=frac(δ·a_hi): δX ≤ 2.5u·δ·a_max. Phase-model error ≤ 2π·S0·(δψ + P·δX).
- Terms: c_n = n^{-1/2}·TB[j]·cis(−2πr), |r|≤2^-11. Taylor truncation is about 1e-25. Relative error budget ρ0=16u, and after k multiplications by e_n, ρ_{k+1}=ρ_k+u(1+ρ_k).
- Bucket sums: γ_cmax·Σ|c|, where Σ|c| ≤ S0·2^-k(1+ρ)/(1−ρ) and |e_n|≤1/2. FFT: unchanged Higham Thm 24.2 bound.
- F accumulation: |fac_k| ≤ (π/2)^k/k!. Per term: F_k(√2·err_k + (S_k+√2·err_k)·4(k+2)u), plus γ_{Lt+2} times the sum of the terms.
- Assumes IEEE-754 binary64 with round-to-nearest and correctly rounded sqrt/div/fma; there is no libm sin/cos in the certified path. These are Higham-style hand bounds, not formally verified.
Still open from before: second source for Gabcke's constant; Trudgian's; no independent reimplementation.

## Cost roadmap, 3e12→1e13 (3.076e13 zeros, $0.03/core-hr)
Before: 0.31 ms/zero (C_PORT) → about $79k by this count (C_PORT quoted $94k). Now, 0.090–0.107 ms/zero → 7.7–9.1e5 core-hr → **about $23–27k**.
(Short-window extrapolation; real runs add Turing-zone overlap and I/O, so budget about $30k.)

## Platt–Trudgian check (arXiv:2004.09765, p.3, verified from the PDF)
"some 7.5 million core hours on 3.6GHz Intel Xeon processors, so each GHz-hour processed a piece of the half line of length about 110 000";
Gourdon "525 days on a single 2.4GHz Pentium 4, so about 80 000 000 … per GHz-hour"; "Gourdon's computation was 725 times quicker than ours". Confirmed.
Zeros to H=3.0001753e12: ≈1.236e13.
- Platt–Trudgian: 7.5e6 h / 1.236e13 = 2.18 ms/zero/core → **7.9 ms·GHz per zero**.
- Gourdon: 8e7 half-line per GHz-hour × ~4.1 zeros/unit → **≈0.011 ms·GHz per zero** (not interval-certified).
- Ours: 0.090–0.107 ms × 2.4 GHz → **0.22–0.26 ms·GHz per zero**. That is about 30× fewer GHz-cycles per zero than P–T and about 20× more than Gourdon.
Caveats: theirs averages over all heights and ours is top-height short windows; 2020 vs 2026 CPUs (instructions per clock differ, so GHz-normalizing is crude);
we use double plus hand error bounds where P–T used mostly Arb; and their total includes production overheads.
