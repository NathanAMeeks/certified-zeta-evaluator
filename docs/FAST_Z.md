# Fast heuristic Z(t) evaluators [HEUR] (fastz.py, fastz2.py, fastz_results.json)
The reference is direct Riemann–Siegel with the C0 correction in float64; phases T·log n mod 2π were precomputed in mpmath. 10,000 grid points per height.
- (1) Truncated RS sum (first N/10 or N/3 terms): FAILS. Sign-error rate is 18–23%, and it misses 6–21% of the sign changes.
- (3) Taylor expansion of the RS sum about centers 2H apart, reused across nearby t (K ≈ e·H·range(ln n)+25 terms):
| t | RS terms | speedup | max abs err | median rel err | sign errors | sign changes (ref) |
|---|---|---|---|---|---|---|
| 1e4 | 39 | 0.4× | 9.5e-13 | 6e-14 | 0 | 921 (921) |
| 1e6 | 398 | 1.6–2.3× | 4.6e-12 | 1.3e-13 | 0 | 953 (953) |
| 2.68e11 | 206,393 | 122× | 1.1e-6 | 3.5e-12 | 0 | 155 (155)* |
- (2) FFT/chirp batch (Odlyzko–Schönhage-lite): not built. It is the same idea as (3) at scale.
Chudnovsky 1e-14 relative criterion, with fallback to exact evaluation when (10·max err)/|Z| > 1e-14: 100% of points are flagged at every height, so speedup including fallback is ≤1×. float64 cannot reach 1e-14 relative error on Z ~ O(1) sums of 2e5 terms, and the float64 reference can't verify it either. Meeting it would need double-double or mpmath arithmetic in both evaluators.
Sign-safe criterion (flag if |Z| < 10·max err): 0% flagged at 1e4/1e6 and 0.04% at 2.68e11. Speedup is 116× at 2.68e11 including fallback.
Verdict: Taylor-reuse gives 0 sign errors on 3×10^4 points (this sample can only show a rate below ~3e-5, not 1e-6) and large speedups at big t. The 14-digit criterion is not met in double precision. HEUR only: it ranks and locates zeros; certification still needs error-bounded (Arb/Gabcke) evaluation.
*Odlyzko table lists 148 zeros in [T, T+40], but both evaluators see 155 sign changes. That points to a defect in my float RS reference at this height (phase precision or the θ Taylor step) or in how I parsed the table. Unresolved, so the 2.68e11 accuracy numbers are only self-consistency vs the reference, not vs truth.

## Follow-up: reference fix + truncated-sum tail predictor [HEUR] (fix.py, tail.py, tail_results.json)
(0) The 155 vs 148 mismatch is RESOLVED, and the reference was not at fault. The table begins at zero #10^12+1, at T+1.8475. The 7 extra sign changes all sit in [T, T+1.63], below the table's start, so they are zeros #≤10^12. Every one of the table's 148 zeros matched a sign change. mpmath.siegelz (dps 25) agrees with the float RS reference at the extra points to ~1e-10. The earlier 2.68e11 accuracy numbers stand.
(1) Tail predictor:
- Analytic: no shortcut exists. The tail Σ_{m<n≤N} n^{-1/2}e^{-it ln n} has stationary points at n=t/(2πk). By stationary phase it maps to a dual sum of length ~N(N/m−1), which is longer than the tail itself for m ≤ N/2.
- Fitted: |error| from a random-phase RMS model (2·sqrt(Σ1/n)), or max |error| on a training half, times safety s=1,2,4. Tested on a held-out half (5,000 points per case).
(2) Results, best cases:
| t | keep | best zero-sign-error setting | fallback | net speedup |
|---|---|---|---|---|
| 1e4 | N/3 | fit, s=1 | 97% | 0.75× |
| 1e6 | N/10 | stat, s=4 | 99.7% | 0.90× |
| 2.68e11 | N/10 | fit, s=1 | 91.5% | 1.01× |
| 2.68e11 | N/3 | fit, s=1 | 76% | 0.93× |
Any setting with net speedup >1.02× had sign-error rates of 0.4–1.4% (stat, s=1).
Verdict: FAIL. The dropped tail is O(1), the same size as Z. Every point the method can trust is one where |Z| exceeds an O(1) window, which is rare, so 63–100% of points fall back and there is no net speedup. Taylor reuse (above) stays the only working idea, at 116× at 2.68e11 with 0.04% fallback.
