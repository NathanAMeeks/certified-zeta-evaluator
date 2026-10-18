# Batched certified Z(t), Odlyzko–Schönhage style [NUM*] (fftb.py, fftb_run.py, fftb_results.json)
Method: I batched over the grid points directly, not over Taylor centers, since one transform covers every point. For u_p=pδ, the main sum F(u_p)=Σ a_n e^{-i u_p ln n} is a nonuniform DFT in x_n=δ ln n/2π mod 1.
- Each x_n is snapped to the nearest node of an M-point grid, M = 2^⌈log2 2P⌉, and the offset is handled by a Taylor-NUFFT: L numpy FFTs. The shift p=p'+P/2 keeps |2πp'd_n| ≤ π/4.
- The order L is the smallest with S0(π/4)^L/L!·e^{π/4} ≤ 1e-4×Gabcke, which gives L = 11–17.
- Phases T·ln n mod 2π use Arb. The θ(t) Taylor coefficients per 2-unit block use Arb, with a Lagrange remainder enclosed on the block.
- The bound is NUFFT truncation + float rounding + θ remainder + C0 rounding + Gabcke. A sign counts only when |Z| > bound; otherwise the point falls back to all-Arb direct RS.
| t | points | span | certified | uncertified after fallback | batch time | cert2 auto_HK (est.) | speedup vs cert2 | vs Arb direct | sign changes | table check |
|---|---|---|---|---|---|---|---|---|---|---|
| 1e6 | 10,000 | 500 | 100% | 0 | 0.04 s | 0.94 s | 24× | 128× | 953 | 953 = cert2 |
| 1e6 | 131,072 | 6,554 | 99.95% | 65 (at zeros) | 0.29 s | 12.3 s | 43× | 236× | 12,491 | — |
| 1e8 | 10,000 | 200 | 100% | 0 | 0.016 s | 2.5 s | 160× | 3,200× | 527 | 527 = cert2 |
| 1e8 | 131,072 | 2,621 | 100% | 0 | 0.21 s | 33 s | 157× | 2,800× | 6,917 | — |
| 2.68e11 | 10,000 | 40 | 100% | 0 | 0.27 s | 29 s | 107× | 8,300× | 155 | 148 table + 7 below its start |
| 2.68e11 | 131,072 | 524 | 100% | 0 | 0.49 s | 385 s | 790× | 61,000× | 2,042 | 2,035 table + 7 below its start |
The cert2 and Arb-direct times are scaled linearly from 10k-point and 100-point timings. Spot checks: 100 random points per run against all-Arb direct RS gave 0 disagreements.
Remaining rigor caveats:
1. The NUFFT/FFT path is float64 with a hand-derived error bound, not Arb balls. It uses a Higham-type summation bound (N+L+4)u·S0 and an FFT bound of 5u·log2M·√M·‖G‖₁. It also adds phase-rounding terms for the P/2 shift and for the x_n mod 1 reduction. The bounds are conservative but need a referee to check them.
2. The θ Taylor series and C0 are evaluated in float with explicit slack, and C0 is flagged when |cos2πp|<1e-3.
3. Gabcke's constant 0.127 t^{-3/4}, t≥200 (Gabcke 1979; see docs/GABCKE_SOURCES.md).
4. A grid of certified signs gives only a LOWER bound on the zero count. Proving that no zeros were missed between grid points still needs Turing's method / an S(t) bound, which was not run here. The table agreement at 2.68e11 is empirical confirmation only.
Verdict: batching works. It is 24–790× faster than cert2 auto_HK and up to 6×10^4× faster than certified direct RS, with a ≥99.95% certified fraction and exact agreement with Odlyzko's table. The gain grows with N and with the number of points per transform. It is NUM*, not NUM, until caveat 1 is replaced by Arb (or a referee accepts the float bound).

## v2: referee, fully Arb path, Turing completeness (fftarb.py, fftarb_run.py, fftarb_{a,a2,c}.json, REFEREE_FFT.md)
(1) Referee: there are 5 flaws in the float bound (see REFEREE_FFT.md). The worst are a dead `*0` term, a float-vs-exact point mismatch not budgeted (≈ Gabcke size at 2.68e11), and an unproven pocketfft twiddle constant. Against Arb, observed errors sat 294–31,000× below the bound, with 0 violations. Even so, the float bound is not a proof.
(2) Fully Arb (python-flint acb.dft, 96-bit precision): this covers the phases, the x_n mod 1 reduction, the offsets e_n, the scatter, the DFTs, the Taylor-NUFFT combination with Arb truncation tail, the θ Taylor series + Lagrange remainder, C0, and Gabcke. The final enclosure is rigorous, except for the Gabcke and Trudgian constants.
| t | points | Arb time | vs float | certified | vs Arb direct | vs cert2 auto_HK |
|---|---|---|---|---|---|---|
| 1e6 (20/unit) | 131,072 | 9.1 s | 32× slower | 99.95% | 8× | 1.4× |
| 1e6 (80/unit) | 131,072 | 8.3 s | 49× slower | 99.90% | 7× | — |
| 2.68e11 | 131,072 | 13.6 s | 29× slower | 100% | 2,300× | ~28× (385 s) |
(3) Completeness: Turing's method in integral form uses Trudgian (Math. Comp. 80, 2011): |∫_{t1}^{t2}S| ≤ 2.067+0.059 log t2 for t2>t1>168π. Certified sign changes give n(t). θ is integrated in Arb with a remainder. The argument then bounds N(T+b) from above via [b,b+L] and N(T+a) from below via [a−L,a].
| window | L | sign changes found | N(a) ≥ | N(b) ≤ | max possible | proven |
|---|---|---|---|---|---|---|
| 1e6+[500,6000], 20/unit | 250 | 10,482 | 1,748,099 | 1,758,587 | 10,488 | NO |
| 1e6+[300,1300], 80/unit | 250 | 1,906 | 1,747,718 | 1,749,624 | 1,906 | YES |
| 2.68e11+[60,460] | 50 | 1,558 | 10^12+226 | 10^12+1,784 | 1,558 | YES (= Odlyzko table 1,558; N(a)=10^12+226 also matches the table) |
At 20/unit, up to 6 zeros are unaccounted for: they may be close pairs missed between grid points or at uncertified points, or simply bound slack. Turing correctly refused to certify the count. At 80/unit the count is proven.
Verdict [NUM]: the fully Arb batched evaluator plus Turing proves exact zero counts: 1,906 zeros at 1e6 and 1,558 at 2.68e11. At 2.68e11 it is ~2,300× faster than Arb direct RS and ~28× faster than cert2. Arb costs ~30–50× over the float path. Remaining caveats are the cited constants only: Gabcke 1979 (RS remainder) and Trudgian 2011 (S(t) integral). The grid must also be dense enough, and Turing detects when it is not.
