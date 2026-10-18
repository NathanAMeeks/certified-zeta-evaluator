# SPEED2.md: certified band-limited interpolation for the Z evaluator (candidate, unrefereed; no RH/Λ claims)
Code: c/interp/fftzi_v7.c (developed as refine/speed2/fftzi.c: fftzsp.c with all referee patches, plus Lemma BL, plus the v7 midpoint probe). Build: c/interp/build.sh (set FLINT_ROOT).
Run: `./fftzi T delta P A B Lw /dev/null dbl` (P must be a power of 2, otherwise EPH=1e300 and nothing is certified. That is a safe failure.)
Env: NOINTERP=1 gives the old behaviour. DIPTAU (default 3), NDIV (8), BLTGT (92), MAXR (20 interp rounds).
Gabcke: uses the published 0.127·t^{-3/4}, t ≥ 200 (tightened 2026-10-28 from the cautious 0.127·(t/2π)^{-3/4}, about 3.97× looser; see docs/GABCKE_SOURCES.md). v7/v8 timings below were measured with the looser form.
Status: candidate, unrefereed. On GitHub main since 2026-10-28.

## 1. Ranked options (gain = wall-time factor on top of the current build)
| # | Idea | Gain | Certification path | Effort | Status |
|---|---|---|---|---|---|
| 1 | Certified Gaussian-windowed sinc (Whittaker–Shannon) interpolation of the RS main sum from a coarse FFT grid. Replaces O(N)-cost refinement evaluations and lets δ grow from 0.015–0.02 to 0.06–0.08. This is the Odlyzko–Schönhage/Gourdon/Platt idea. | **2.5–3× [NUM]** | Lemma BL below (Poisson summation; explicit erfc aliasing bound and Gaussian tail), plus grid ball and hand double-rounding bound | done (prototype) | §2 |
| 2 | **(done in v2, §5)** Extremum-targeted refinement: use the interpolated Z' (differentiate the kernel) to put a point at each local extremum of \|Z\| between grid points, as Gourdon/Platt do. Then δ≈0.12–0.2 (≈1–2 grid pts/zero) should prove windows. At δ=0.12 today the parabola/deficit targeting misses 3–5 pairs per window. | ~1.5–2× [HEUR] | Targeting is heuristic. Soundness unchanged (certified signs + Turing) | 1–2 days | next |
| 3 | Per-point θ, C0, Gabcke in double with hand-proven bounds (that phase is now ~30% of run time and still uses Arb), plus a tighter EPH (measured slack ~10^5–10^6×, §2) | ~1.3–1.4× [HEUR] | Same Higham-style analysis as SPEED.md; needs a referee | 1–2 days | open |
| 4 | Lemma "one zero / no zero" via certified Z' and a Z'' bound from the BL representation (avoids fine sign sampling near close pairs) | small now (refinement is <3% of time) [HEUR] | provable with the same kernel bounds | 1 day | low priority |
| 5 | GPU: FFT and per-point work with directed-rounding intrinsics, or ball arithmetic (CUDA __dadd_rd etc.; CAMPARY, Joldes–Muller–Popescu–Tucker 2016) | hardware ~10×/$ possible [HEUR] | rounding-mode-correct kernels; reuse the existing hand bounds | weeks | later |
| 6 | NUFFT with ES kernel (FINUFFT; Barnett–Magland–af Klinteberg, SIAM J. Sci. Comput. 41 (2019) C479–C504; aliasing analysis Barnett, ACHA 51 (2021)) instead of Taylor-NUFFT (L≈18 FFTs) | ≤1.5× on the NUFFT phase only [HEUR] | needs a fully rigorous kernel-error bound (Barnett 2021's bound: check whether it is fully explicit; unverified) | 2–4 days | maybe |
| 7 | Error-free transformations / double-double (Ogita–Rump–Oishi, SIAM J. Sci. Comput. 26 (2005)) | accuracy, not speed; useful for item 3 bounds | standard | — | supporting |
| 8 | Fast multipole with rigorous bounds | ~none: the 1-D exponential sum is already handled by FFT/Taylor | — | — | skip |
| 9 | Hiary's sublinear methods (Annals of Math. 174 (2011) 891–946, O(t^{1/3+ε}) per value) | impractical for bulk verification at 10^12–10^13: they target single values at much larger t, and our amortised cost is already O(1) FFT work per point [HEUR] | — | — | skip |
Gourdon's "Taylor-in-shift" multi-evaluation is essentially what our Taylor-NUFFT already does (Odlyzko–Schönhage, Trans. AMS 309 (1988) 797–809), so there is no separate gain there.
Sources: Platt, Math. Comp. 86 (2017) 2449–2467 (interpolation + Turing; exact form of his interpolation lemma not re-checked here [unverified]); Platt–Trudgian, Bull. LMS 53 (2021), arXiv:2004.09765; Gourdon (2004 preprint, "The 10^13 first zeros…").

## 2. Lemma BL [LEM, candidate; own derivation, needs referee]
F(u)=Σ_{n≤N} n^{-1/2} e^{-i(T+u)log n}. g(u)=e^{iβu}F(u), β=½log N, has all frequencies in [-β,β]. Grid u_k=kδ, A=π/δ>β,
K(x)=sinc(Ax)·e^{-x²/2h²}. Poisson summation for one exponential e^{iωs}, |ω|≤β, gives Σ_k e^{iωkδ}K(u−kδ)=Σ_m e^{i(ω−2Am)u}p(ω−2Am),
where p(ξ)=Pr(|ξ+Y|≤A), Y~N(0,h^{-2}). Hence the aliasing+window error is ≤ E_al=erfc(z1)(1+1/(1−e^{−2z1D})), z1=(A−β)h/√2, D=√2Ah
(using erfc(x+y)≤erfc(x)e^{−2xy}). Truncation to |j|≤K (X=(K−½)δ): E_tr=2e^{−X²/2h²}/(AX(1−e^{−Xδ/h²})). Then
|F(u) − Σ_{|j|≤K} F̃_k e^{−iβx_j}K(x_j)| ≤ S0(E_al+E_tr) + √2·r_grid·Σ|K(x_j)| + fp, with S0=Σn^{-1/2} (Arb upper bound) and fp a hand bound
2(ε_t+γ_{2K+2})Σ|F̃_k||K_j|, ε_t=((K+3)²+20(K+1))u (recurrences for the Gaussian, cis and sinc; sin(A x_j)=(−1)^j sin(Ax0) exactly).
The transcendental seeds (sin, sinc, exp, cis) come from Arb at 128 bits. θ, C0 and Gabcke are added in Arb as in fftzsp.
Parameters: h=√(X/(A−β)), X(A−β)≈92, so K=36–41 and E_al+E_tr ≈1e-22. The error is dominated by the grid ball (≈2e-7) times Σ|K|≈4.5.
Check [NUM]: 20,000 interpolated Z values at 3e12 (δ=0.06 and δ=0.12 grids) vs a direct δ=0.01 grid: max |diff| 3.2e-12,
certified radius ≈1.2e-6 (bound/actual ≈ 4e5), all signs agree.

## 3. Results [NUM] (2.4 GHz Xeon, shared box, load 6–10)
Chains (single windows A=60, B=31020, Lw=50), run in the same batch:
| chain | old fftzsp δ=0.015 P=2^21 | new fftzi δ=0.06 P=2^19 |
|---|---|---|
| 3e12 (132,506 zeros) | 13.30 s, 0.100 ms/zero, 16.0 evals/zero, PROVEN | 5.38 s, **0.041 ms/zero**, **4.06 evals/zero**, PROVEN, 0 disagreements |
| 1e13 (138,439 zeros) | 14.54 s, 0.105 ms/zero, 15.3 evals/zero, PROVEN | 5.83 s, **0.042 ms/zero**, **3.90 evals/zero**, PROVEN, 0 disagreements |
On a lighter-loaded run: 4.45 s / 4.51 s (0.034 / 0.033 ms/zero). Turing endpoints are identical to before (12362403240765→…373271; 43124192297371→…435810).
SPEED.md's δ=0.02 best (0.090/0.107 ms/zero) is ≈11–12 grid evals/zero. Its old "≈5" label was wrong and is now corrected in SPEED.md.
Per GHz: 0.08–0.10 ms·GHz/zero, vs Platt–Trudgian 7.9 (~80–100× fewer cycles, crude) and Gourdon ~0.011 (still ~7–9× behind; Gourdon was not interval-certified).

10 test windows (A=60, B=10380, Lw=50; δ=0.08, P=2^17): 0, 20000, 40000, n1, n2, n3, m1, m2, m3 are all PROVEN in 2.0–3.2 s each
(0.045–0.069 ms/zero, 3.2–3.25 evals/zero; old closepair runs took 15–51 s at δ=0.01). The Turing bounds are identical to the old runs, with 0 disagreements.
c1: δ=0.08 is UNPROVEN (2 short, one close pair missed, nothing false). The retry at δ=0.04 (P=2^18) is PROVEN with 44,169 zeros in 2.0 s (6.9 s wasted on the failed try).
Policy: on UNPROVEN, retry at δ/2. Cost with retries is still well under the old figure.
δ=0.12: Turing bounds tighten to the exact values, but 3–6 zeros are still missed per window, which leaves it UNPROVEN. Deficit-loop overhead is O(want·gaps), so it ran slowly. This is item 2.
Cost roadmap (3e12→1e13), extrapolated [HEUR]: 0.034–0.042 ms/zero → about $9–11k at $0.03/core-hr (was $23–27k). Budget about $12–13k.

## 4. Open referee items for this build
Lemma BL proof check (Poisson-summation applicability and the erfc tail lemma); the fp bound ε_t (hand derivation); whether Platt 2017's published lemma matches (if so, cite it instead).
Existing items carry over (Gabcke wording, independent reimplementation).

## 5. v2: extremum-targeted refinement (idea 2) [NUM, candidate]
Targeting only (heuristic). Soundness is unchanged: only certified Lemma BL signs count, together with Turing and Lemma CP.
1. Z'(u_p) at every grid point from the BL kernel derivative. At u=u_p the kernel collapses to F'(u_p) ≈ −iβF_p + Σ_{j≠0} F_{p−j}e^{−iβjδ}(−1)^j G(jδ)/(jδ),
   so Z' = 2Re(e^{iθ}(iθ'F + F')). This is double precision and uncertified, with θ' saved from the per-point step.
2. Cubic Hermite (Z, Z') on each grid interval, sampled at 33 points. The interval is a candidate if it predicts more sign changes than the grid shows,
   has an interior same-sign local min of |H| below KTAU·(endpoint scale) (KTAU=1), or touches an uncertified grid point.
   That is about 3% of intervals (4.2k–5k per chain).
3. For each candidate: 65 approximate BL samples (precomputed weight table, so 84 complex MACs each), corrected linearly to the certified endpoint values (this absorbs C0).
   Then one certified Lemma BL evaluation at the midpoint between each pair of predicted roots. This is Platt 2017's "estimate the location, then check rigorously" pattern.
   Two bugs found during the build: Hermite-argmin probing alone missed close pairs, and omitting C0 (≈1e-3) from the approximate samples gave false root counts next to grid points with |Z|≈1e-4.
4. Leftover deficit goes to the old deficit loop (MAXR=20), then the δ/2 retry policy. Neither was needed in any v2 run below.
NOEXT=1 gives v1 (parabola dip pass). The old fftzi_v1.c is kept.

Same batch (6 jobs on 8 vCPU, load ≈8). All runs PROVEN, 0 disagreements, Turing endpoints identical:
| chain | old fftzsp δ=0.015 | v1 δ=0.06 | **v2 δ=0.12 (P=2^18)** |
|---|---|---|---|
| 3e12, 132,506 zeros | 14.02 s, 0.106 ms/zero, 16.1 ev/zero | 6.17 s, 0.047, 4.06 | **3.40 s, 0.0256 ms/zero, 2.01 ev/zero** |
| 1e13, 138,439 zeros | 14.86 s, 0.107 ms/zero, 15.4 ev/zero | 6.64 s, 0.048, 3.90 | **3.49 s, 0.0252 ms/zero, 1.93 ev/zero** |
Lighter load: v2 2.42 s / 3.01 s (0.018 / 0.022 ms/zero). v2 time split (3e12): setup 1.17 s (O(N) Arb tables), NUFFT 0.95, per-point Arb 0.98, targeting+probes 0.26.
In cycles: 0.04–0.06 ms·GHz/zero. That is ~130–190× fewer cycles than Platt–Trudgian's 7.9 (crude) and still ~4–5× behind Gourdon's ~0.011, which was not interval-certified.
10 test windows at **δ=0.16 (P=2^16)**: all 10 PROVEN (0, 20000, 40000, c1, n1, n2, n3, m1, m2, m3), 0 disagreements, Turing bounds identical to the old closepair runs.
They took 1.2–2.3 s each (0.028–0.051 ms/zero; higher than the chains because the fixed setup cost is spread over 3× fewer zeros), at 1.50–1.60 evals/zero. 4 windows used 1 deficit round, and none needed a retry.
δ=0.24 chains: UNPROVEN (28–36 short; Hermite screening at 0.24 is too coarse; nothing false was certified).
Cost (3e12→1e13), extrapolated [HEUR]: 0.018–0.026 ms/zero → about $5–7k at $0.03/core-hr. Budget about $8k.

## 6. Platt 2017 vs Lemma BL (checked against the author-accepted manuscript, Bristol repository; Math. Comp. 86 (2017) 2449–2467, doi:10.1090/mcom/3198)
Platt §4.2 does rigorous Whittaker–Shannon up-sampling (his Thm 4.3), with Weiss's theorem (Thm 4.4: error ≤ 4∫_B^∞|F(x)|dx for nearly
band-limited f). He applies it to the Gaussian-windowed W(t)=Λ(t)·exp(πt/4 − (t−t0)²/2H²), with explicit bounds deferred to his Appendix C (not re-derived here).
Zeros are located non-rigorously (Newton–Raphson) and then confirmed by two rigorous evaluations, one either side.
**Same family, different lemma.** Platt windows the *function* (Λ·Gaussian is only nearly band-limited, so a Weiss tail bound is needed).
Lemma BL uses the *exact* band-limitedness of the finite RS main sum and puts the Gaussian in the *kernel* (Poisson-summation aliasing bound via erfc).
So BL is not a restatement of a published lemma. It is a close analogue, and it is still unrefereed. Platt 2017 is the right prior-art citation for the method.
His implementation parameters (H=2089/16384, 70 samples each side using every fifth sample) are confirmed from the AAM PDF (REFEREE_BL_DP.md). Platt 2017 does not use Riemann–Siegel or Gabcke.

## 7. Next
(a) Setup is now the largest single cost (O(N) Arb log/rsqrt tables per window). Use longer windows (P=2^20 at δ=0.12 covers ~126k units) or cache them across windows. Expect ~1.3–1.5× [HEUR].
(b) Idea 3 (double-precision per-point θ/C0/Gabcke with hand bounds) is still open. It needs a careful rounding analysis and wasn't cheap, so it was not done this round.
(c) Make δ=0.2 work with a better screen (Z'' or approximate samples on more intervals).
(d) Referee: Lemma BL proof, the ε_t bound, and comparison with Platt's Appendix C.

## 8. v3 (2026-10-19): setup cache, double-precision per-point step, δ=0.2–0.24 screen, compiler [candidate, unrefereed]
Build: `gcc -O3 -march=native -ffp-contract=off` (no -ffast-math/-Ofast). Env: AHCACHE=path, LTF (default 100), PH (block half-width, default 8), SCREEN, NZT, OLDPTS=1, DBLTEST=1.
1. **Setup cache [NUM]:** a_n=log n/2π as double-double (AH, AL) with Arb radius ARAD, kept in a file (`ahcache.bin`, 20 MB, FNV-1a checksum, refused if N is too small).
   With the cache, the per-window O(N) Arb loop is skipped, and S0 uses the bound Σn^{-1/2} ≤ 2√N+1. Setup 1.2–1.6 s → 0.06–0.15 s.
2. **Lemma DP [LEM, candidate]: per-point θ, C0, Gabcke in binary64.** Per block (half-width H=8), Arb gives co0 (reduced mod 2π), co1, co2 → doubles with radii cor_k (|co_k−co_kd| upper bounds),
   th3 ≥ |θ'''/6| over the block, C0 as one Arb ball over the t-ball [c−H, c+H] (mid C0m, radius C0r; the block falls back to Arb if c0term refuses), Gabcke at the block's lower end, and the trunc bound tr.
   Per point: v_d=fl(fl(pδ)−cu), |v_d−v| ≤ ev=1.01u(|pδ|+|v_d|), av=|v_d|+ev; θ̃=co0+v(co1+v·co2) (Horner). The error vs the true θ is
   eθ = Σcor_k av^k + av³th3 + |co1|ev + 2av|co2|ev·1.01 + 3.03u(|co0|+av|co1|+av²|co2|).
   e^{iθ̃} uses the existing Arb-certified 1024-entry table plus degree-6/7 Taylor; the reduction f=θ̃·fl(1/2π) adds ≤2.02u|θ̃|, so εcs = 16u + 2.02u|θ̃|.
   zr = 2(Fr·ec − Fi·es) + C0m. Bound E = [2(|Fr|+|Fi|)(eθ+εcs) + 2.02u(|a|+|b|+|a−b|) + 1.01u|zr| + C0r] + 2√2·tr(1+eθ+εcs) + gabcke, times (1+1e-12).
   The sign is certified iff |zr| > E. Every term is a sum of per-operation errors ≤ u·|result|, which is valid only without contraction, hence -ffp-contract=off.
   (GCC's GNU-C default is fp-contract=fast, which would fuse a·b±c once -march exposes FMA. The explicit fma() in the phase code is intended and unaffected.)
   **Adversarial test [NUM]:** 256-bit Arb truth 2Re(e^{iθ(t)}F̂)+C0(t) at every 13th grid point (θ taken directly at t, not from the Taylor series), on the -O3 -march=native build:
   3e12 δ=0.12: 20,165 points, max |err|/bound 0.044; 1e13: 0.020; δ=0.24 runs: 0.037 and 0.020; PH=32: 0.022 and 0.008. **0 violations.**
   Per-point phase 0.8–1.0 s → 0.05 s.
3. **Screen for δ≥0.15 [HEUR targeting]:** 9 approximate BL samples per interval (table weights) in addition to the Hermite test.
   Also flag intervals with a zero within 0.1δ of an end (|Z| < 0.1|δZ'|): a pair hiding right next to a grid point broke window 40000 at δ=0.2.
   Plus certified points exactly at the four Turing-zone boundaries (A−Lw, A, B, B+Lw). At δ=0.24 the boundaries are no longer grid points, and a zero straddling B loosened Nb_upper by 1.
   Plus certified bisection of zone intervals to width 0.02. The grid origin is shifted 20 units down (T−20, A+20, B+20) so the lower Turing zone has K≈62 grid points of interpolation support. The absolute windows are the same.
4. **Compiler/SIMD:** -O3 -march=native, plus an explicit textbook complex product cmul (bit-identical on finite values, but no __muldc3 NaN fallback calls; there were 16 before) in the FFT butterfly, the combine loop, and the BL and screen MAC loops.
   LTF=100 (Taylor-NUFFT truncation ≤100×Gabcke ≈4e-8, still below EPH≈2–3e-7; rigorous for any value) cut L from 18–19 to 14.
   Measured: about 1.1× total (NUFFT 0.55–0.68 → 0.46–0.52 s), within noise on this loaded box. The radix-2 FFT is the remaining hot spot. A hand-vectorised radix-4 AVX2 FFT is needed for more, and its Higham bound would need re-checking.

### Results (same batch, 6 jobs on 8 vCPU, load ≈7.5; all PROVEN, 0 disagreements, Turing endpoints identical)
| chain | old fftzsp δ=0.015 | v3 δ=0.12 (P=2^18) | v3 δ=0.24 (P=2^17) |
|---|---|---|---|
| 3e12 132,506 | 14.46 s, 0.109 ms/zero, 16.1 ev/zero | **1.34 s, 0.0101, 2.02** | 1.53 s, 0.0115, **1.18** |
| 1e13 138,439 | 16.20 s, 0.117 ms/zero, 15.4 ev/zero | **1.21 s, 0.0088, 1.94** | 1.73 s, 0.0125, **1.16** |
Lighter load (v3 before the screen was added, δ=0.12): 0.93 s / 0.90 s (0.0070 / 0.0065 ms/zero).
10 windows at δ=0.2 (P=2^16, origin −20): **all 10 PROVEN**, 0 deficit rounds, 0 disagreements, Turing identical to the closepair runs, 0.49–0.69 s each (0.011–0.016 ms/zero), 1.59–1.64 evals/zero.
At δ=0.24 the evals/zero (1.16–1.18) match Gourdon's ~1.2. Wall time there is dominated by ~25k certified probes (still Arb seeds, ~15–20 µs each), so δ=0.12 is currently fastest.
**ms·GHz:** 0.016–0.024 ms·GHz/zero (light to loaded), vs Gourdon ~0.011 (~1.5–2.2× behind, crude, and Gourdon wasn't interval-certified) and Platt–Trudgian 7.9 (~330–490× fewer cycles, crude).
**Cost (3e12→1e13, 3.076e13 zeros, $0.03/core-hr) [HEUR]:** 0.0065–0.010 ms/zero → 5.6–8.5e4 core-hr → **about $1.7–2.6k**. Budget about $3k (Turing overlap, I/O, retries).

### Next
(a) A double-precision certified probe (BL seeds via the TB table plus rigorous short Taylor for exp/sinc, θ from the block data as in Lemma DP). That would make δ=0.24 the fastest and remove most of the remaining 0.3–1.0 s.
(b) AVX2 radix-4 FFT with a re-derived Higham bound. (c) Longer windows to amortise the Turing zones. (d) Referee: Lemma BL, Lemma DP, ε_t, Platt Appendix C.

## 9. v4 (2026-10-19): double-precision certified extra points (Lemma BL-DP) [candidate, unrefereed]
- `zinterp_d`: the Lemma BL sum and the Lemma DP final step entirely in binary64. Arb runs only once per 8-unit block (θ Taylor data, C0 ball, Gabcke) and at init.
  The interpolant is evaluated at u'=k0δ+x0d, x0d=fma(−k0,δ,u). |u'−u| ≤ ex=1.01u|x0d| costs |F(u)−F(u')| ≤ S0·logN·ex (termwise |F'| ≤ Σn^{-1/2}log n).
  Seeds: sinc(Ax0) by degree-24 Taylor (|Ax0| ≤ π/2, remainder <1e-22), Horner bound γ_{2n+3}·e^{|y|}. exp(−w) seeds (w ≤ 0.5) by degree-12 Taylor, same bound.
  e^{−iβx0} via the Arb-certified TB table (16u + 4.1u|φ|). Per-term relative budget ε_t = ((K+3)² + 60(K+1) + 24 + 4.1|φ|)u, plus the absolute seed errors (s, sinc) times Σ|F||G/(Ax)|.
  If the DP bound doesn't certify, the probe falls back to the Arb probe (337–409 of ~25k).
- **Bug caught by the pipeline:** the first DP build summed exp with signed coefficients in (−w), which gives e^{+w}. Lemma CP flagged it as ERROR_overcount / disagreement and nothing was reported as proven. Fixed.
- **Adversarial test [NUM]:** 256-bit Arb truth of the same interpolant (same K, same grid values), with θ at t directly and C0(t). Over the δ=0.24 chains (3e12: 24,968 probes; 1e13: 28,413), max |err|/bound was 0.021 and 0.006, **0 violations.**
- Targeting tweaks [HEUR]: two-level fine sampling (17 coarse samples, then full resolution only near sign changes or |z| minima) and 4-lane SoA dot products. These are heuristic only; reassociation is harmless there.
  Tried and reverted: sampling only the near-end eighth for "zero near a grid end" candidates. It lost pairs (UNPROVEN, and a 57 s deficit loop).
- **AH cache reproducibility:** regenerated from scratch, `ahcache.bin` is byte-identical (cmp; sha256 prefix c13d3b248f52a6ea).
- **Bernstein route (CREATIVE_SPEEDUPS #1), checked:** for g(u)=Σc_n e^{iω_n u}, |ω_n| ≤ β, the bound |g^{(k)}| ≤ β^k Σ|c_n| = β^k S0 holds termwise and is elementary.
  (Bernstein's inequality proper gives β^k·sup|g|, which is sharper but needs a sup bound.) The termwise version is what the DP probe already uses for the u' shift.
  As a probe method it would need derivatives g^{(k)}(u_p) up to k≈25 (βδ/2≈0.84 at δ=0.24), i.e. more transforms or interpolations. It doesn't beat the BL-DP probe (≈5.5 µs/probe, 0.13–0.15 s per chain), so it wasn't pursued.

### Results (load ≈8–9 on 8 vCPU; all PROVEN, 0 disagreements, Turing endpoints identical)
| chain | δ=0.12 (P=2^18) | δ=0.24 (P=2^17, origin −20) |
|---|---|---|
| 3e12 | 0.0071–0.0074 ms/zero, 2.02 ev/zero | 0.0073–0.0077 ms/zero, **1.18** ev/zero |
| 1e13 | 0.0071–0.0078 ms/zero, 1.94 ev/zero | 0.0080–0.0081 ms/zero, **1.16** ev/zero |
δ=0.24 now **ties** δ=0.12 but doesn't beat it. Its NUFFT saving (~0.3 s) is eaten by targeting: fine sampling ≈0.25 s for ~57k candidate intervals, the screen ≈0.13 s, probes ≈0.14 s.
10 windows: all PROVEN at δ=0.24 and at δ=0.2 (P=2^16, origin −20), 0 deficit rounds, Turing identical to closepair, 0.41–0.61 s each (≈0.010–0.013 ms/zero).
ms·GHz ≈0.017–0.019 (loaded) vs Gourdon ~0.011 (not certified; ~1.6× behind, crude). Old fftzsp in the same session: 0.109–0.117 ms/zero, so ≈15×.
Cost 3e12→1e13 [HEUR]: ≈0.0072 ms/zero → 6.2e4 core-hr → **≈$1.9k**. Budget about $2.5k.
Next: cut the candidate count at δ=0.24 (the "zero near a grid end" rule flags ~43% of intervals; a Z''-aware or derivative-sign test could replace it). Then a vectorised FFT if it can be done cheaply with a re-derived bound.

## 10. v5 (2026-10-19): sharper screen for δ>0.15 [HEUR targeting only; soundness unchanged: certified signs + Turing + Lemma CP]
- Z'' at every grid point from the BL kernel's second derivative at x=jδ: K''(jδ)=(−1)^j G(−2/x²−2/h²) for j≠0, K''(0)=−A²/3−1/h²; Z''≈2Re(e^{iθ}(−θ'²F+2iθ'F'+F'')).
- A quintic Hermite (Z, Z', Z'') per interval replaces the cubic. The interval is flagged if it predicts extra sign changes, has an interior same-sign |H| min <0.5·end scale,
  or |H| stays below 0.1× the neighbouring amplitude (a flat near-zero stretch hid a pair at 3e12 u≈22897).
- The "zero within 0.1δ of a grid end" rule is removed. The 9-sample approximate screen now runs only when the quintic shows a small non-crossing value (<0.25× local amplitude; ~65% of intervals).
- The quintic alone missed one triple-zero cluster at 1e13 (u≈28654.9; quintic error ~δ⁶ω⁶/46080·amp ≈3% of amplitude at δ=0.24), so the gated sample screen stays in.
- **Flagged intervals:** 43% → **18.5–21.2%** (chains), 18.4–21.3% (windows δ=0.24), 10.4–12.3% (windows δ=0.2).
### Acceptance (load ≈7–9; all PROVEN, 0 disagreements, 0 deficit rounds, Turing identical)
| chain | δ=0.12 | δ=0.24 |
|---|---|---|
| 3e12 | 0.0074 ms/zero, 2.02 ev/zero | **0.0066–0.0069**, 1.18 |
| 1e13 | 0.0073–0.0075, 1.94 | **0.0070**, 1.16 |
10 windows: all PROVEN at δ=0.24 and δ=0.2 (0.45–0.67 s each), Turing equal to closepair.
δ=0.24 is now the fastest. That's ≈0.016–0.017 ms·GHz/zero vs Gourdon ~0.011 (≈1.5× behind, crude, not certified); ≈16× the old fftzsp.
Cost 3e12→1e13 [HEUR]: ≈0.0068 ms/zero → 5.8e4 core-hr → **≈$1.7k** (budget about $2.3k).

## 11. Rounding proofs for Lemma BL-DP (zinterp_d) and Lemma DP (closes REFEREE_BL_DP.md gaps) [LEM, candidate; the referee is an AI, so a human referee is still needed]
Model: IEEE binary64, round-to-nearest, u=2^-53, -ffp-contract=off, so each +,−,×,/ rounds once, |fl(x∘y)−x∘y| ≤ u|x∘y|. The explicit fma calls round once. γ_k=ku/(1−ku).
First-order counts below are absorbed by the stated multipliers (≥1.01) and the final (1+1e-12).

(a) Perturbed abscissa. x0d=fma(−k0,δ,u) gives |x0d−x0| ≤ u|x0d|(1+u) ≤ ex. The BL interpolant is evaluated exactly at u*=k0δ+x0d.
Lemma BL holds at any real point, and sin(A(x0d+jδ))=(−1)^j sin(A·x0d) is exact in the reals (A=π/δ). |F(u)−F(u*)| ≤ ex·Σn^{-1/2}log n ≤ ex·S0·logN.

(b) Seeds (Taylor + Horner). With rounded coefficients and z=fl(y²), Horner gives |p̂−p(z)| ≤ γ_{2n}Σ|c_k||z|^k (Higham 2002, eq. 5.3).
Coefficient rounding adds u·Σ, and the error in z adds ≤ n·u·Σ. Total ≤ γ_{2n+3}·p̃(|z|), p̃ ≤ e^{|y|} (coded with n=12, γ_29).
- sinc: degree 24, |y| ≤ 1.6 (checked), truncation <1e-22 (added). y=fl(fl(π/δ)·x0d) has |y−A·x0d| ≤ 2.02u|y|=ey, and |sinc'| ≤ ½.
  So esc = γ_29·e^{|y|} + 0.5ey + 1e-22. s=fl(y·ŝc) gives es = |y|·esc + 1.01u|s| + ey.
- exp(−w): degree 12, signed coefficients, 0 ≤ w ≤ 0.5 (checked), truncation ≪ u. Relative error ≤ γ_29·e^{w}/e^{−w} ≤ 2.2e-15 ≈ 20u.
  The argument w carries ≤3u relative error, which moves exp by ≤3u·w (relative). Budgeted as 40(K+1)u in ε_t, because R seeds are reused up to K times.
- cis: Arb-certified table entries rounded (≤u each), degree-6/7 Taylor for |x| ≤ 2π·2^-11 (truncation <1e-25, Horner ≤4u), complex product ≤√5u, x=fl(r·fl(2π)) ≤2u|x|.
  Reduction f=fl(φ·fl(1/2π)): ≤2.02u|φ| rad. f−floor(f) is exact for f ≥ 0 (Sterbenz). **For −1<f<0 it equals f+1 with one rounding (≤u/2 turn = πu rad).**
  Total ≤ (2+4+√5+π+0.01)u + 2.02u|φ| ≤ 16u + 3.15u + 2.02u|φ|. Coded: final cis (16+3.15)u + 2.02u|θ̃| (both paths); E0 28u + 4.1u|φ| (includes β_d rounding).

(c) ε_t, per-term relative error of term j (|j| ≤ K). Term_j = F_k·E_j·K_j, E_j = E0·w^j, K_j = (−1)^j s·G_j/(A_d·x_j), G_j = G0·Π_{i<j}R_i, R_i = R0·Q^i.
- R_i carries ≤ (2i+1)u, so G_j ≤ u + Σ_{i<j}(2i+2)u = (j²+j+1)u, plus seed error 40(j+1)u. That is ≤ (j+3)² + 40(j+1).
- E_j: the w seed ≤√2u per use, plus j complex products ≤√5u each: 3.7j·u, plus E0's 28u + 4.1u|φ|.
- K_j: x_j=fl(x0d±fl(jδ)) relative ≤ 6u (|x_j| ≥ (|j|−½)δ); A_d, product and division: 3u. Total 9u. The absolute error of s goes into es·SB (SB=Σ|F_k|G_j/(A|x_j|)), and j=0 into esc·|F_k0|G0.
- f·E·K_j: √5u + u ≤ 3.3u.
Sum ≤ (j+3)² + 40(j+1) + 3.7j + 12.3 + 28 + 4.1|φ| ≤ ((K+3)² + 60(K+1) + 28 + 4.1|φ|)u = ε_t (coded). The sum of 2K+1 terms adds γ_{2K+2}·Σ|terms|.
**Efp = 2(ε_t + γ_{2K+2}(1+ε_t))·SA·1.01 + es·SB·1.02 + esc·|F_k0|·G0·1.02 + S0·logN·ex** (the factor 2 is a deliberate margin over the first-order count).
Arb fallback zinterp: Arb seeds (≤u + radius), ε_t = ((K+3)² + 20(K+1))u. That covers the referee count (j²+j+1) + 3.6j + 5.

(d) |F(u)−Ŝ| ≤ Eint = Efp + S0(E_al+E_tr) + √2·r_grid·Σ|K_j|·1.01.

(e) Final step (shared with Lemma DP). Horner θ̂ = co0 + v(co1 + v·co2) has 4 roundings, so the error is ≤ γ_4(|co0| + av|co1| + av²|co2|) (textbook; replaces 3.03u).
v_d: grid |v_d−v| ≤ 1.01u(|pδ|+|v_d|); probe 1.01u|v_d|. Mean-value term ev·((|co1|+cor1) + 2av(|co2|+cor2)) (cor1/cor2 now included).
Plus Σcor_k·av^k and Taylor av³·th3, th3 ≥ sup|θ'''|/6 on the t-ball. **eθ ×1.25 margin** (the referee found the Taylor term attained to 0.999).
|zr−Z| ≤ 2(|Fr|+|Fi|)(eθ+εcs)·1.000001 + 2.02u(|a|+|b|+|a−b|)·1.01 + 1.01u|zr| + C0r + 2·Eint(1+eθ+εcs) + Gabcke (grid path: 2√2·tr in place of 2·Eint).
The sign is certified iff |zr| exceeds that bound × (1+1e-12).
**Rerun after these bound changes [NUM]:** both δ=0.24 chains (132,506 / 138,439; 0.0064 / 0.0069 ms/zero) and all 10 windows at δ=0.24 are PROVEN, 0 disagreements, Turing identical. Speed is unchanged within noise.

## Context reference (not used in any proof here)
OpenAI's quasi-RH release proves that ζ(s) and every Dirichlet L-function has no zeros with Re s > 7/8 (by the functional equation, all nontrivial zeros lie in 1/8 ≤ Re s ≤ 7/8). This is far from the critical line, so it changes none of the bounds or lemmas this evaluator uses. Status: the statement is machine-checked in Lean against their proof; human review is ongoing. [LEM, external]

```bibtex
@misc{OAI:The-Quasi-Riemann-Hypothesis-September-30-2026,
  author = {{OpenAI}},
  title = {{The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane $\mathrm{Re}(s)>7/8$}},
  howpublished = {OpenAI Math Release preprint, 30 September 2026,
    \url{https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf}},
  note = {Repository commit 1581d4c; Lean statement lean/ComparatorChallenges/QuasiRiemannHypothesis.lean. Companion (11/12 version): OAI:The-Quasi-Riemann-Hypothesis-October-5-2026. Lean-checked; human review ongoing.},
  year = {2026}
}
```

## 12. v6 (2026-10-21): jitter and hard-window test [NUM, candidate, unrefereed; soundness unchanged]
Code: fftzi_v6.c (v5 kept as fftzi.c). Build: v6/build_v6.sh. Raw data and scripts: v6/ (bench.sh, summ.py, base_*.txt, v6_pin.txt, delta_scan.txt).
Recovered run args (they reproduce the stored Turing endpoints exactly): windows `T 0.24 65536 80 10400 60 /dev/null dbl`; chains `T 0.24 131072 80 31040 60 /dev/null dbl`; AHCACHE=ahcache.bin.
Box: 8 vCPU Xeon, load ≈0.8–1.4 (lighter than §9–10's 7–9, so absolute numbers are lower than before). 3 reps per run.

Changes (search/targeting only; the certified probes, Lemma BL/BL-DP/DP and Turing are untouched):
- Quintic/Hermite basis polynomials precomputed once (same expressions, so H is bit-identical).
- The 9-sample screen's raw values are reused by the candidate fine sampler at i≡0 (mod 8). Same formula and same p, so the values are bit-identical, and 9 of the 17 coarse samples are skipped for screened candidates.
- **Check:** found, Na_lower, Nb_upper, evals/zero and status are identical v5↔v6 on all 36 runs (diff empty). All PROVEN, so Lemma CP holds (certified sign changes = Turing count) for all 10 windows and both chains.

| window | v5 unpinned min/med/max | v5 pinned | v6 pinned | evals/zero (all) | PROVEN |
|---|---|---|---|---|---|
| 0 | 0.0087/0.0089/0.0098 | 0.0075/0.0078/0.0084 | 0.0075/0.0075/0.0076 | 1.708 | 3/3 each |
| 20000 | 0.0076/0.0076/0.0079 | 0.0075/0.0079/0.0081 | 0.0074/0.0074/0.0076 | 1.706 | 3/3 each |
| 40000 | 0.0076/0.0076/0.0076 | 0.0076/0.0076/0.0079 | 0.0073/0.0073/0.0085 | 1.706 | 3/3 each |
| c1 | 0.0075/0.0082/0.0087 | 0.0075/0.0076/0.0077 | 0.0073/0.0077/0.0079 | 1.704 | 3/3 each |
| m1 | 0.0082/0.0085/0.0089 | 0.0081/0.0083/0.0085 | 0.0080/0.0081/0.0084 | 1.659 | 3/3 each |
| m2 | 0.0083/0.0086/0.0088 | 0.0080/0.0082/0.0085 | 0.0080/0.0082/0.0087 | 1.66 | 3/3 each |
| m3 | 0.0082/0.0090/0.0091 | 0.0083/0.0085/0.0089 | 0.0081/0.0081/0.0082 | 1.661 | 3/3 each |
| n1 | 0.0076/0.0078/0.0079 | 0.0076/0.0076/0.0084 | 0.0074/0.0074/0.0075 | 1.706 | 3/3 each |
| n2 | 0.0077/0.0085/0.0088 | 0.0076/0.0077/0.0080 | 0.0072/0.0074/0.0075 | 1.706 | 3/3 each |
| n3 | 0.0078/0.0078/0.0080 | 0.0075/0.0077/0.0078 | 0.0073/0.0074/0.0086 | 1.708 | 3/3 each |
| chain3e12 | 0.0058/0.0063/0.0064 | 0.0057/0.0059/0.0059 | 0.0055/0.0057/0.0059 | 1.183 | 3/3 each |
| chain1e13 | 0.0060/0.0062/0.0064 | 0.0059/0.0060/0.0061 | 0.0059/0.0059/0.0061 | 1.157 | 3/3 each |

Summary:
- **Pinning** (`taskset -c 7`) cut the mean per-window spread from 8.0% to 6.0%, and the worst chain from 0.0064 to 0.0061.
- **v6 code** gives about 1–3% on medians. Chains (best/worst): 0.0055/0.0061 (v6 pinned) vs 0.0057/0.0061 (v5 pinned) vs 0.0058/0.0064 (v5 unpinned).
- **Occasional outliers** (40000, n3 at 0.0085) are still OS noise even when pinned.
- **The structural gap** between 3e12 and 1e13 is about 4–8%. It comes from N∝√t (setup and NUFFT 0.25 vs 0.27 s; probes 23.8k vs 27.2k) and is not reducible by targeting.
- **Windows vs chains:** windows sit ~25% above chains because of fixed per-run cost; that is amortization, not hard windows.

Not done or no gain:
- **Per-height δ from zero density:** δ scan on the chains (v6/delta_scan.txt):
  - δ=0.25 needed 1 deficit round (0.0066–0.0068).
  - δ=0.26 gives 0.0059–0.0061, the same as 0.24.
  - δ=0.28 fails at 1e13 (UNPROVEN_budget, nothing false certified).

  Keep δ=0.24 at both heights.
- **Batching close-pair refinement:** at δ=0.24 refine_points=0 and rounds=0, and probes are already one pass, so there is nothing to batch.
- **madvise/hugepages and ahcache preload:** not tried, because setup is only 0.05 s (~6%) and the file is already page-cached on repeat runs.
- **Biggest remaining cost:** the screen loop (t_loop ≈0.36–0.39 s of 0.77–0.82 s). The next real lever would be cheaper Zp/Zpp stencils or vectorizing the 9-sample mac4. The stencil-shortening lemma was not attempted, per instruction.

## 13. v7 (2026-10-21): quintic fast path ahead of the fine sampler [NUM, candidate, unrefereed; soundness unchanged]
Code: fftzi_v7.c (v5 fftzi.c and v6 kept). Build: v6/build_v7.sh. Data: v6/v6_pin_b.txt, v6/v7_pin.txt. Runs pinned with `taskset -c 7`, 3 reps each, same args as §12.

**Profile (which search parts change the outcome):**
- Almost every probe comes from the quintic Z/Z'/Z'' fit.
- The 9-sample screen adds only 2 candidates per 1e13 chain. Dropping it (SCREEN=0) costs a deficit round at 1e13, so it stays.
- KGATE, KT5 and KFLAT did not change the probe set within noise, so they were left at v5 defaults.

**Change (targeting only):** when the grid shows no sign change across an interval but the quintic shows exactly 2, v7 first places one certified probe at the midpoint of the two quintic roots.
- It is accepted only if the certified sign is opposite to both grid endpoints. That certified sign change brackets the pair.
- Otherwise the interval falls back to the v6 fine sampler.
- Result: 26,234 / 23,021 tries per chain, 100% hit, so fine-sampler candidates drop from 27.4k/24.0k to 1.2k/1.0k.
- A missed fast-path probe can only cost a fallback or a deficit round, never correctness. Certification is still the certified sign changes plus Turing (Lemma CP). Disable with NOFASTH=1.

**Check:** found, Na_lower, Nb_upper, evals/zero and status are identical between v6 and v7 on all 36 runs. All are PROVEN.

| run | v6 pinned min/med/max ms/zero | v7 pinned min/med/max | evals/zero |
|---|---|---|---|
| 0 | 0.0078/0.0080/0.0081 | 0.0069/0.0071/0.0071 | 1.708 |
| 20000 | 0.0075/0.0076/0.0077 | 0.0071/0.0072/0.0073 | 1.706 |
| 40000 | 0.0076/0.0079/0.0079 | 0.0068/0.0069/0.0071 | 1.706 |
| c1 | 0.0075/0.0078/0.0079 | 0.0069/0.0070/0.0070 | 1.704 |
| m1 | 0.0081/0.0084/0.0085 | 0.0076/0.0078/0.0079 | 1.659 |
| m2 | 0.0085/0.0086/0.0088 | 0.0075/0.0077/0.0078 | 1.66 |
| m3 | 0.0084/0.0085/0.0087 | 0.0076/0.0077/0.0078 | 1.661 |
| n1 | 0.0077/0.0079/0.0080 | 0.0070/0.0071/0.0073 | 1.706 |
| n2 | 0.0078/0.0079/0.0082 | 0.0069/0.0070/0.0075 | 1.706 |
| n3 | 0.0076/0.0079/0.0080 | 0.0069/0.0070/0.0074 | 1.708 |
| chain3e12 | 0.0057/0.0057/0.0060 | 0.0050/0.0050/0.0054 | 1.183 |
| chain1e13 | 0.0060/0.0061/0.0066 | 0.0052/0.0052/0.0055 | 1.157 |

**Summary:**
- Chains (best/worst): v6 0.0057/0.0066 → v7 0.0050/0.0055 ms/zero, about 12–15% faster.
- Windows: about 8–10% faster.
- Spread: per-window spread is unchanged (about 4.7%).
- Evals/zero: unchanged, because the probes are the same certified evaluations at different (equivalent) points.

**Tried and not kept:**
- Early exit at found==Turing: not applicable, because the Turing count needs the zone search first, and found coarse is about 76%.
- Lazy Z'': not applicable, because the quintic needs Z'' on every interval.
- Priority ordering (smallest-|Z|-first, sorter idea): not implemented [HEUR: expected no gain], because with zero deficit rounds every candidate is visited anyway.

**Remaining cost:** probes (≈0.12 s, certified, required), the 9-sample screen, and setup/NUFFT.

**Credit (§13):** The certified midpoint sign-flip probe used in v7 is Nathan A. Meeks' idea, implemented here.
Prior record: N. A. Meeks, "Phase-Coherent Holographic Compression," X Article, Feb 18, 2026, https://x.com/NathanAMeeks/status/2024058294244180131

## §14 v8: delta cache extender + gated self-correction [NUM, candidate, unrefereed]

Directed by Nathan A. Meeks (delta cache extender, gated self-correcting loop); implemented here. Proof logic unchanged (Lemma CP, Gabcke/Trudgian, DP error bounds, Lemma BL untouched); all changes are search/setup only.

**Diagnosis (confirmed in code).** v7 `ah_load` rejects the AH cache when stored N < floor(sqrt(T/2π)) (1e13 cache: N=1,261,566; 1.5e13 needs 1,545,096; 2e13 needs 1,784,124). zgrid then rebuilds all N Arb logs from scratch (~1.7–1.9 s) and overwrites the file.

**Delta cache extender (`fftzi_v8.c`).** Load a short cache as a prefix, compute only n = N_old+1..N_target with the same Arb ops at PREC=96, take ARAD = max(old, new), and save. Target = N(T(1+AHXS)) rounded up to AHXC. Checks: extending the 1e13 cache to 2e13 gives a file byte-identical (`cmp`) to a fresh 2e13 build. Found counts, Na and Nb match v7 at 1.5e13 and 2e13. Setup on a cache miss drops from 1.68 s to 0.06 s + 0.10 s extension (1.5e13) and from 1.86 s to 0.06 s + 0.19 s (2e13). `NOAHX=1` restores v7 behaviour.
Squeezed-delta sweep (v8/sweep.txt): 8 sliding windows across 1.5e13, 2e13, 2.5e13, with AHXS ∈ {0, 1e-3, 1e-2, 0.1} and AHXC ∈ {1, 65536}. All settings give one extension per height jump (0.07–0.14 s). Slack makes no measurable difference, because N grows by about 0.0005 per 10^4 of height. Chosen: AHXS=0, AHXC=65536.

**Self-correcting loop (`v8/v8run.sh`).** If a window comes back not PROVEN, the driver re-runs it from scratch at δ = 0.2, then 0.16, then 0.12, doubling P when needed to keep the span. Every attempt is a full independent certification, and only a PROVEN result is accepted. This catches a hard window found during the sweep, T=20000000020800: v7 gives UNPROVEN_budget (47284 of 47286, also with NOFASTH and MAXR=60). v8 proves it at δ=0.2: 47286 = Turing, PROVEN. The failed first attempt costs about 1.6 s of wall time, so the retry is not free.

**Gated self-tuner (`v8/tune.py`, v8/tune.txt).** Knobs: KGATE, KT5, KFLAT, NDIV, DIPTAU, MAXR, δ. Each was run on the 12 bench windows/chains × 3 pinned repeats. A knob is kept only if every run is PROVEN with identical found/Na/Nb AND the median drops by more than the noise (spread of the per-repeat baseline medians, 0.0005 ms/zero).

| Knob | Result |
|---|---|
| δ=0.22, δ=0.23 | 6 of 36 runs fail their proof → rejected |
| KGATE=0.2 (0.00775) | Closest miss, still inside noise |
| All other knobs | Inside noise |

KEPT = none, so v8 ships v7's search defaults.

**Before/after (pinned, 3 reps, v8/final.txt):**

| set | ver | min | med | max ms/zero | setup med s | evals/zero | PROVEN |
|---|---|---|---|---|---|---|---|
| 12 bench | v7 | .0056 | .0082 | .0106 | .038 | 1.705 | 36/36 |
| 12 bench | v8 | .0078 | .0085 | .0095 | .036 | 1.706 | 36/36 |
| 6 new (4e12…2e13 + hard 2e13+20800) | v7 | .0079 | .0089 | .0353 | .054 | 1.673 | 15/18 |
| 6 new | v8 | .0079 | .0087 | .0124 | .054 | 1.657 | 18/18 |

Notes:
- Bench timings are the same within noise. v8's win is robustness: cache misses (about 5x faster setup on the first window at a new height) and self-corrected proofs.
- In the v8 rows of the table, ms/zero for the retried window counts only the attempt that succeeded.
- The final table was run on the already-extended cache, so it does not show v7's cold-miss penalty. That penalty is shown above (1.7–1.9 s vs 0.16–0.25 s).
