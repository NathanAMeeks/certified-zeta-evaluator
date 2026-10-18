# Certified Taylor-reuse Z(t) [NUM] (cert.py, cert_results.json)
The certified bound has three parts:
- Taylor remainder: S0·(R|v|)^K/K!·e^{R|v|}, where S0=Σn^{-1/2}, R=max|ln n − c| and c is the center of ln n. The order K is chosen so this is ≤ 0.1×Gabcke.
- Gabcke's RS remainder bound after the C0 term: 0.127(t/2π)^{-3/4}.
- Higham-style float rounding: S0·e^{R|v|}·(δφ+(N+2K+8)u), with u=2^-53. The phases T·ln n mod 2π come from mpmath at 40 digits, and θ(t) from mpmath at each point, so θ has no Taylor remainder.
A sign is certified when |Z̃| > bound. Points near the removable singularity of C0 (|cos 2πp|<1e-3) are flagged. Fallback is the direct RS sum with the same Gabcke + rounding bound. Arb's own zeta is infeasible at 2.68e11, because its cost grows like t.
10,000 grid points per height. Best H (Taylor radius) per height:
| t | H | K | certified | uncertified after fallback | max bound | speedup (incl. θ) | speedup (sum only) | sign changes vs table |
|---|---|---|---|---|---|---|---|---|
| 1e4 | 4 | 34 | 99.93% | 7 (also uncertified by direct) | 5e-4 | 0.99× | 0.7× | 945 = 945 |
| 1e6 | 4 | 54 | 100% | 0 | 3e-5 | 1.15× | 5.3× | 953 (no table) |
| 2.68e11 | 2 | 63 | 99.99% | 0 | 8.7e-3 | 64× | 163× | 148 = 148 table zeros (+7 below table start) |
At 2.68e11, H=4 fails: the e^{R|v|} rounding growth leaves only 70% certified. H=1 gives tighter bounds (2e-5) at 58×.
There were 0 sign disagreements with direct RS anywhere. Per-point θ via mpmath is the bottleneck at low t; computing θ by its own certified Taylor series would remove it.
Verdict [NUM]: the rigorous version keeps nearly all of the heuristic's gain at large t. It is 64× faster than certified direct RS at 2.68e11 (163× on the sum), with every sign certified after fallback and exact agreement with Odlyzko's table. There is no gain at 1e4, where N is only 39. Caveats: the float rounding bound is a hand-derived Higham bound, not Arb, and Gabcke's constant is taken from the literature (t ≥ 200). Both need a referee check.

## v2: all-Arb (cert2.py, cert2_{a,b,c}.json) [NUM]
Everything is in Arb (python-flint 0.9, 110-bit precision):
- The phases n^{-iT}, the per-center Taylor coefficients P_k=Σ b_n w_n^k, and the per-point series evaluation all use ball arithmetic, so there is no hand-derived rounding bound.
- θ(c+v) uses Arb series coefficients at the center up to v², plus a Lagrange remainder bounded by |θ'''/6| enclosed over the ball [c±H] (acb_series lgamma). There is no per-point mpmath.
- The Taylor remainder S0(R|v|)^K/K!·e^{R|v|} is added to the ball, and C0 is evaluated in Arb.
- Baseline: the certified direct RS sum, fully in Arb, plus Gabcke. Its per-point time is measured on 3 points.
| t | H | K | certified | fallbacks / uncertified | Taylor time | Arb direct /pt | speedup | sign changes vs table |
|---|---|---|---|---|---|---|---|---|
| 1e4 | 2 | 20 | 99.93% | 7 / 7 | 0.54 s | 0.12 ms | 2.2× | 945 = 945 |
| 1e6 | 2 | 32 | 100% | 0 / 0 | 0.87 s | 0.50 ms | 5.7× | 953 (no table) |
| 2.68e11 | 2 | 63 | 100% | 0 / 0 | 32 s | 0.237 s | 73× | 155 = 148 table + 7 below table start |
At 1e4, the 7 uncertified points are |Z| below Gabcke's 5e-4 margin, right at zeros. The direct method can't certify them either. None of them hides a sign change.
Remaining caveat: only Gabcke's constant, |Z(t) − (RS main sum + C0 term)| < 0.127·t^{-3/4} for t ≥ 200 (the code uses the looser 0.127·(t/2π)^{-3/4}). Source: W. Gabcke, "Neue Herleitung und explizite Restabschätzung der Riemann-Siegel-Formel", Dissertation, Göttingen, 1979, doi:10.53846/goediss-5113. Restated in A. Kuznetsov, "Series expansions for the Riemann zeta function", arXiv:2312.03261. (An earlier version of this note wrongly cited Platt–Trudgian; their paper, Bull. London Math. Soc. 53 (2021), does not state it.) See `docs/GABCKE_CHECK.md`.
Verdict [NUM]: certified sign evaluation is 73× faster than all-Arb direct RS at 2.68e11. Every grid sign is certified, and the sign changes match Odlyzko's table exactly. The speedup grows with N, roughly like N/K.

## v3: auto-tuner (auto_HK in cert2.py; tune.py, tune_results*.json) [NUM]
Rule: choose H to minimise the modelled cost per point, (N·K + 5500)/(2H·ρ) + 10K, where ρ is grid points per unit t. K is the smallest order with Taylor remainder ≤ 0.1×Gabcke, which gives K ≈ e·R·H + log(S0/tol), with R = ½ln N. Precision is prec = 64 + log2 T + log2 N bits, which is 82–119 bits here. The constants were calibrated from timings of this python-flint code, so the rule is implementation-specific.
| t | tuned H | K | prec | certified | speedup (tuned) | speedup (H=2) | Taylor time tuned vs H=2 |
|---|---|---|---|---|---|---|---|
| 1e4 | 1.86 | 19 | 82 | 99.93% (7 at zeros) | 1.8× | 2.6× | 0.53 vs 0.50 s |
| 1e5 | 1.49 | 21 | 87 | 99.91% (9 at zeros) | 3.3× | 3.5× | 0.60 vs 0.59 s |
| 1e6 | 1.34 | 24 | 92 | 100% | 4.9× | 10.7× | 0.94 vs 0.83 s |
| 1e8 | 2.07 | 44 | 102 | 100% | 19.0× | 16.5× | 2.51 vs 2.51 s |
| 2.68e11 | 2.9–6.2 | 82–156 | 119 | 100% | 76–86× | 72× | 27–29 vs 30 s |
Speedups use a 3-point timing of the direct method, so they are noisy (±2× at low t). Taylor time is the reliable comparison. Fitted optimum: H ≈ 1.3–2 for t ≤ 1e8, rising to ~3–6 at 2.68e11. The cost minimum is flat, so H=2 is within ~10% of optimal everywhere. A first tuner without the per-center overhead term chose H=0.1–0.4 at low t and was 2–10× slower, which is why that term is in the model.
Sign changes match the table at 1e4 (945 = 945) and 2.68e11 (148 + 7 below the table start). Lower precision did not reduce the certified fraction.
Verdict: the tuner works but gains only ~10%, at large t. Fixed H=2 is already near-optimal. Further speedup needs batching (FFT over centers, Odlyzko–Schönhage), not tuning.
