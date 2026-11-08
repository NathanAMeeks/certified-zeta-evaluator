# Riemann–Siegel N-step fix in v9p — [LEM] candidate, unrefereed
Grok for Nathan A. Meeks, Oct 28 2026. Patched copy c/interp/fftzi_v9.c (`c0term`); /workspace/cze_repo untouched. Not RH.

**Problem.** fftzi_v9 sets N = floor(√(T/2π)) once per window and builds F_N(t) = Σ_{n≤N} n^{-1/2} e^{-it log n} by NUFFT/BL interpolation.
Every assembly site computes Z ≈ 2Re(e^{iθ(t)}F_N(t)) + c0term(t,N) ± Gabcke. Past t_N = 2π(N+1)², p = √(t/2π) − N ≥ 1, so c0term refused
those points (ok=0). A window crossing t_N could never be PROVEN (window 519: found 38839 vs a bogus Turing upper count of 44297).

**Lemma NSTEP.** For t with k = ⌊√(t/2π)⌋ − N ≥ 1, set N' = N+k and p' = p − k ∈ [0,1). The standard RS formula with the correct index is
Z(t) = 2Σ_{n≤N'} n^{-1/2} cos(θ(t) − t log n) + (−1)^{N'−1}(t/2π)^{-1/4}Ψ(p') + R, with |R| ≤ 0.127 t^{-3/4} (Gabcke, t ≥ 200; unchanged).
Since 2Re(e^{iθ}F_N) = 2Σ_{n≤N} n^{-1/2} cos(θ − t log n), this equals the existing assembly once c0term returns
  C(t) = (−1)^{N'−1}(t/2π)^{-1/4}Ψ(p') + 2Σ_{n=N+1}^{N'} n^{-1/2} cos(θ(t) − t log n).
The patched c0term does exactly that, in Arb, on the caller's t-ball. Here θ(t) = Im logΓ(1/4 + it/2) − (t/2) log π (acb_lgamma), and Ψ(p') comes from the C0 cut (C0_FIX.md).
It needs ⌊p⌋ constant on the ball; a ball that straddles t_N returns ok=0, and the caller's existing per-point fallback handles those few points.
The error bounds keep their form: C is an Arb enclosure, used exactly where the old C0 ball was (C0m/C0r in the binary64 paths, an arb add elsewhere).
For t < t_N nothing changes except one extra comparison, which is why normal windows run at the same speed. NSTEPOFF=1 restores the old refusal, for A/B tests only.

**Validation [NUM].**
- Window 519 (T=3000006078280, crosses t_N ≈ 3000006087434.6): PROVEN on the first attempt in 0.53 s, all 65536 grid points certified, 44170 zeros = Turing count.
- Future crossing window T=3000014765664 (t_N ≈ 3000014770664.8, N 690989→690990), tested in scratch: PROVEN on the first attempt, 44171 zeros.
- Z(t) spot check vs mpmath.siegelz (25 digits), 10 grid points from t−t_N = −1155 to +1145, including ±0.3–5 around t_N: max |diff| 4.1e-9, within the evaluator's per-point error budget (EPH ≈ 1e-7).
- Benchmarks (10 bench windows, 2 chains, 2e13+20800, 2 reps each): 26/26 PROVEN, counts identical to unpatched v9; wall 12.90 s vs 12.94 s.
- Driver: cheap pre-check flags windows crossing t_N (logs/nstep.log) and skips the slow MAXR retry rung for them.
