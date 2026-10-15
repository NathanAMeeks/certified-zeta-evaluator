# C/FLINT port of fftarb (c/fftz.c) [NUM]
Build: gcc against the Debian FLINT 3.6 build (same flags as c/build.sh).
Usage: `./fftz T delta P a b Lw [signfile] [dbl]`.
Pipeline:
- Arb phases, x_n, e_n, and the Arb scatter, with L Taylor-NUFFT passes.
- The DFT either uses acb_dft (Arb), or, with the `dbl` flag, a radix-2 double FFT with a rigorous enclosure. That enclosure covers input rounding + rad(G) as ‖ΔG‖₁, plus Higham's FFT forward-error bound, Thm 24.2 in *Accuracy and Stability of Numerical Algorithms* (2nd ed., 2002): log2M·η/(1−log2M·η)·√M·‖G‖₂ with η=μ+γ4(√2+μ). Twiddles come from Arb rounded to nearest, μ=2u. The result is added as a ball radius.
- The combine, θ Taylor + remainder, C0 and Gabcke are all in Arb. Turing/Trudgian is computed in Arb inside C.
Not ported: the per-point Arb-direct fallback. It was never triggered in these runs, with 100% certified.
## Validation vs Python (fftarb.py)
| window | Python proven | C (acb_dft) | C (dbl FFT) | sign arrays (131,072 pts) |
|---|---|---|---|---|
| 2.68e11+[60,460] | 1,558, N(a)=10^12+226 | same | — | — |
| 3e12+[60,460] | 1,712, N(a)=12,362,403,019,918 | same | same | 0 disagreements (both C modes) |
| 1e13+[60,460] | 1,789, N(a)=43,124,192,297,371 | same | same | 0 disagreements |
## Benchmarks (single core per process)
| config | time | sign changes | ms/zero | RSS |
|---|---|---|---|---|
| Python fftarb, 3e12, 131k pts | 19.8 s | 2,245 | 8.8 (250/unit) | — |
| C acb_dft, 3e12, 131k pts | 11.6 s | 2,245 | 5.2 | — |
| C dbl, 3e12, 131k pts | 3.65 s | 2,245 | 1.6 (250/unit) | 196 MB |
| C dbl, 3e12, 2^20 pts, 100/unit | 13.8 s | 44,877 | 0.31 | 826 MB |
| C dbl, 3e12, 2^21 pts, 100/unit | 26.0 s | 89,749 | 0.29 | 1.55 GB |
| C dbl, 1e13, 2^20 pts, 100/unit | 16.5 s | 46,888 | 0.35 | ~830 MB; **proven 44,201 zeros** |
Memory is ≈ 740 B/point, so a 2^22 window is ≈3 GB. That is the max feasible window on this 16 GB box, which has ~3.8 GB free. The algorithm itself has no limit.
Multi-core: 1, 2 and 4 concurrent processes (2^18 pts) each took 5.3 s wall, so throughput scales linearly to 4 cores. The limit is memory, not CPU.
Completeness: windows of 2^20–2^21 points at 3e12 (60 and 100/unit) were NOT proven, with 2–18 zeros unaccounted for. Larger windows contain rare close pairs that slip between grid points. Production needs a local refinement step (extra certified evaluations in gaps flagged by Turing), which is not built. It should cost about +1–10%. The 1e13 2^20 window was proven without it.
## Roadmap with measured C rate (replaces ÷20 guess)
Rate fit: ms/zero ≈ 0.31 + 0.077·log10(t/3e12), from 0.31 at 3e12 and 0.35 at 1e13. Costs include +10% for refinement/overlap, at $0.03/core-h. [HEUR extrapolation]
| verify all zeros 3e12→T | zeros | ms/zero at T | core-h | $ |
|---|---|---|---|---|
| 1e13 | 3.1e13 | 0.35 | 3.1e6 | ~$94k |
| 3e13 | 1.2e14 | 0.39 | 1.35e7 | ~$406k |
| 1e14 | 4.6e14 | 0.43 | 5.5e7 | ~$1.7M |
| 1e15 | 5.0e15 | 0.50 | 7.3e8 | ~$22M |
The measured C gain is ~5–6× over Python, not 20×, so the ÷20 guess was ~4× too optimistic.
Next speedups:
1. Move the per-point combine (P·L Arb addmul) and the scatter to double, with rigorous a-posteriori bounds. These dominate now; est. 3–10×.
2. Larger windows with refinement.
3. Lower precision (64 bits) where the bounds allow.
Verdict: the C port is validated against Python, with identical proven counts and identical certified signs. It is 5.4× faster on the same window, and ~7× per zero in production-size windows (0.31 ms vs Python’s 2.3 ms at similar density). It proved 44,201 consecutive zeros at 1e13 in 16.5 s on one core. A full verification to 1e13 would cost ~$94k at today's rate, so about $10–30k needs the double-precision combine plus refinement first.
