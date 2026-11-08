# Hybrid verifier: results [NUM]. Credit: Nathan A. Meeks & Grok
These are candidate results. They are unrefereed and make no RH claims. Everything runs on a copy of /workspace/rh_extend, and the original is untouched. The box was shared with Math Researcher's live run (load ~7–9 on 8 cores), so expect about ±15% timing noise.
Contents: gates.py and spotcheck/indgates.py are patched behind env flags (all off by default, so a plain run reproduces the baseline). bench.py is the harness. Logs are B.log, Hy.log, ref.log, refh.log, pred.log, and per-window brackets are in out/.
Run: `V9BIN=$PWD/v9p/fftzi_v9 V9RUN=$PWD/v9p/v9/v9run.sh HYB_AHC=1 HYB_VC=1 HYB_TMESH=1 [HYB_REUSE=1] python3 bench.py TAG idx...`

## Profile (window 10202)
- Each close-pair probe or chunk run of fftzi_v9 took 0.72 s, of which sec_setup was 0.68 s. gates._run never passed AHCACHE, so every probe recomputed the coefficient table from scratch. A window runs 7–9 probes.
- The referee (ind_gate on the spliced chunk) took 66 s: G-Turing 35.7 s (748 Arb evaluations) and G-RS 30.7 s (563 endpoints). Some points were evaluated twice.

## Per-idea table
| Idea | Measured | Correctness | Verdict |
|---|---|---|---|
| DPA#4 widest-gap / widest-bracket prediction (HYB_PREDICT) | +1.5–2 s per window; 0/9 hits | PROVEN and counts unchanged | **KILL** (slower) |
| DPA#2 reuse the probe chunk certificate instead of a final rerun (HYB_REUSE) | fires on ~3/9 windows; saves ~0.2–0.7 s there | same counts and N; bracket endpoints differ inside the chunk (delta-0.24 certificate vs delta-0.02), both certified by Lemma CP | **small WIN, optional** |
| DPA#1 interval Newton on Z', #3 Taylor model, #5 block-margin check | not implemented (C-level Z''/remainder proofs) | — | **NOT TESTED** |
| Error-correction: pass AHCACHE to gate probes (HYB_AHC) | probe 0.72→0.17 s; close-pair windows 6.5–12.8 s → 2.1–6.2 s | same binary path the main run already uses; bracket files byte-identical | **WIN** |
| Referee: shared Arb value cache plus a meshless G-Turing attempt at bracket endpoints first (HYB_VC, HYB_TMESH) | 66.4 → 32.1 s (G-Turing 35.7→15.5, G-RS 30.7→16.6) | identical N(lo), N(hi), 282/282 certified, RESOLVED. The mesh is search-only; only coinciding Trudgian bounds are accepted, otherwise it falls back to the old configuration | **WIN** |

## Hybrid (AHC + REUSE + VC + TMESH) vs baseline, 14 windows, GATE_IND=0
| window | type | base s | hybrid s | PROVEN both | N/Turing same | brackets byte-same |
|---|---|---|---|---|---|---|
| 10202 | cp | 9.68 | 4.15 | ✓ | ✓ | ✓ |
| 8758 | cp | 12.81 | 6.15 | ✓ | ✓ | ✓ |
| 35009 | cp | 6.75 | 2.23 | ✓ | ✓ | endpoints differ (REUSE) |
| 35540 | cp | 9.29 | 4.92 | ✓ | ✓ | ✓ |
| 36208 | cp | 6.99 | 2.42 | ✓ | ✓ | ✓ |
| 16098 | cp | 6.70 | 2.34 | ✓ | ✓ | ✓ |
| 21377 | cp (slow, not fast path) | 108.5 | 100.3 | ✓ | ✓ | ✓ |
| 35423 | cp | 6.45 | 2.07 | ✓ | ✓ | endpoints differ (REUSE) |
| 36497 | cp | 7.13 | 2.36 | ✓ | ✓ | ✓ |
| 20000 | cp | 7.18 | 2.16 | ✓ | ✓ | endpoints differ (REUSE) |
| 35010, 35541, 36209, 30000 | clean | 0.33–0.36 | 0.34–0.37 | ✓ | ✓ | ✓ |
Fast-path close-pair windows: **2.6× mean (2.1–3.3×)**. Clean windows: no change. 14-window total: 182.9 → 130.5 s (1.4×; dominated by the outlier 21377, which was not investigated).
Referee on 10202: **2.07×**. Estimated whole-run effect, assuming the close-pair repair share of 38% gets 2.6× faster: about 1.3× overall [HEUR], more if the referee is counted.

## Caveats
- With REUSE on, bracket endpoints are not byte-identical (counts and N are). Use HYB_REUSE=0 for byte-identical output. Its gain is small.
- The referee comparison is n=1 window. Run more before adopting it.
- Upstream adoption (Math Researcher's call): in gates._run, add AHCACHE to the env, plus the two indgates.py flags.

## DPA#1 / #3 / #5 follow-up (2026-10-31, Nathan A. Meeks & Grok) [NUM]
Probes: /workspace/meeks_test/refine/dpa/{probe.py,jet.py}, at t≈3.0001e12 (window 10202). Nothing in rh_extend or cze_repo was touched. None of the three ideas passed the cost gate, so none was wired into gates.py. The 14-window counts are therefore unchanged from the hybrid run (14/14 match Turing), and no harness rerun was needed.
| Idea | Measured | Correctness | Verdict |
|---|---|---|---|
| DPA#1 interval Newton on Z' (certified Z'') | Arb Z jet via `acb_dirichlet_hardy_z` takes 0.6 s for len=1 (Z only) and **>200 s (timed out) for len=3 (Z, Z', Z'')**. pair_cert v2 already certifies a pair with 3 Arb signs (~0.2–0.6 s each), so the Z'' route is ≥100× slower per pair. A cheap RS-series bound gives \|Z''\| ≤ 2·Σn^-½·(½log(t/2π))² ≈ 6.0e5 plus R'', and there is no rigorous R'' bound, so it would be search-only. | not adopted; counts unchanged (14/14 match Turing) | **KILL** |
| DPA#3 Taylor model (one Arb ball over [t0±h], Cauchy bound) | Arb Z on a complex ball of radius R gives a real-part radius of ≈9e6·R (R=1e-6→9.0, 1e-4→9e2, 1e-2→1.4e5) in 0.4–0.55 s, i.e. an effective Lipschitz bound of ≈9e6 against a true \|Z'\| of ~1e2. So with \|Z\|~2 it certifies a sign only on h ≲ 2e-7, versus a grid step of 0.02: 0 blocks certified | not adopted; counts unchanged (14/14) | **KILL** |
| DPA#5 block-margin (sign on a block if \|Z(t_i)\| > L·h + err) | The ceiling is the point-check share: 0.018 of 0.271 s (7%) on a clean window. The certified RS Lipschitz bound is L ≈ 2·Σn^-½·½log(t/2π) ≈ 4.5e4 (R' not rigorously bounded), so the margin needs h < \|Z\|/4.5e4 ≈ 1e-4, which is finer than the grid. It cannot merge grid points, so the gain is 0 (≤5% even ideally) | not adopted; counts unchanged (14/14) | **KILL** |

## speed2: profile + prototypes (2026-10-31, Nathan A. Meeks & Grok) [NUM]
Copies are in /workspace/meeks_test/refine/speed2/. Nothing in rh_extend or cze_repo was touched, and nothing was pushed. Box load was 11–17 on 8 cores, so timings are noisy. G-Turing and G-RS stay independent: HYB_VC (shared values) is OFF in every referee run.
**Profile.**
- Clean window 35010 (CPU 7): total 0.60 s = setup 0.06 + NUFFT/FFT 0.19 + point checks 0.03 + dip/BL-interp search 0.27 (fine-sample 0.16, probes 0.16, 260 Arb fallbacks). Interpolation/dip search is the largest share.
- Referee (TMESH on, 3 procs): 10202 took 75 s (G-Turing 31 s for 330 Arb evals; G-RS 44 s for 563 endpoints). 21377 took 45 s (G-Turing 30 s for 238 evals; G-RS 15 s for 188 endpoints). It is roughly 100% Arb Z point evaluations at about 0.2–0.3 s each (python-flint acb.zeta at 128 bits). The Platt fast path never fires because the blocks have 163–282 zeros, under PLATT_MINZ=300.
- 21377: the hybrid close_pair needs a δ=0.02 rerun of chunk [4871,4919] for the narrow pair 4877.430/4877.456/4877.490 (gaps 0.026/0.034). That took 157 s this run (92 s before). pair_cert proves the pair with 3 Arb signs in 4.5 s, which explains Math Researcher's ~31 s live replay.
| Idea | Measured vs baseline | 14/14 match Turing | Verdict |
|---|---|---|---|
| S1 Arb precision ladder starting at 64 bits (sgn/_sv: 64→128→256→512; ball arithmetic, so escalation keeps it certified) | Serial microbench on 80 real points: 20.9→18.0 s (**1.16×**), 0 uncertain. Referee end-to-end: 10202 75→86 s, 21377 45→32 s, which is within load noise | yes (referee RESOLVED on both blocks; finder untouched) | **KILL** (marginal; the gain is below noise) |
| S2 CP-S + Arb pair certificate (unified v2 config) on all 14 windows | 211.8→54.8 s (**3.9×**); 21377 157→4.5 s (35×); 13/14 resolved by pair_cert, and 10202 fell back to close_pair (9.9 s) | yes | **WIN** (Math Researcher already adopted it) |
| S3 multi-window amortization: 4 independent windows in parallel on CPUs 4–7 (no shared state) | 14-window wall 54.8→~18 s (**~3.0× throughput**, from timestamps) | yes (14/14) | **WIN** (throughput only; per-window certification unchanged) |
| S4 referee pool 7 vs 3 procs (with S1) | 21377 32→23 s, 10202 86→84 s; it oversubscribes at load 15 | yes | **KILL** (load-dependent, no reliable gain) |
Next lever (not prototyped) [HEUR]: the referee is dominated by Arb point evals in a short block, so a certified multipoint evaluator would be the big win. FLINT's Platt multieval would do that if OPEN-P is closed; otherwise an Arb-certified batched RS sum over the block. G-RS must use a different evaluator than G-Turing to stay independent.

## mpeval: certified batch Riemann–Siegel evaluator for G-RS (2026-10-31, Nathan A. Meeks & Grok) [NUM], candidate, unrefereed; no RH/Λ claims
Code: refine/mpeval/c/rsbatch.c, hook MPEVAL=1 in mpeval/spot2/indgates.py (copy). Method: Z=2Re(e^{iθ}Σ_{n≤N} c_n e^{-iv log n}) + (−1)^{N−1}τ^{−1/4}C0(p) ± 0.127τ^{−3/4} (Gabcke), c_n=n^{−1/2−it0} in Arb once per chunk, inner sum in double with explicit a-priori bound (coeff/phase/sincos 4 ulp/blocked-sum γ), θ, C0 and the final combination in Arb. Uncertified points fall back to Arb sgn (0 needed). Distinct from G-Turing's acb_zeta; reads no library/finder values. Assumes glibc sin/cos ≤1 ulp (budgeted 4); needs a referee.
| window (chunk, trimmed to 10400) | endpoints | G-RS base (Arb, 3 procs) | G-RS mpeval (3 thr) | speedup | result | per-endpoint signs vs Arb |
|---|---|---|---|---|---|---|
| 10202 [7511,7577] | 563 | 48.6 s | 2.42 s | 20× | RESOLVED, Turing 282 = certified 282 | identical (0 diff) |
| 8758 [8692,8736] | 217 | 15.7 s | 1.43 s | 11× | RESOLVED 188 = 188 | identical |
| 27053 [7831,7873] | 211 | 15.5 s | 1.18 s | 13× | RESOLVED 180 = 180 | identical |
| 38498 [10368,10400] | 160 | 21.2 s | 0.81 s | 26× | RESOLVED 137 = 137 | identical |
Radius ≤3.7e-9 (actual error vs Arb ≈3.4e-12); min |Z| at endpoints 1.25e-6, so 0 fallbacks. Box load 11–18 on 8 cores (the live run was going). G-RS total 101.0 s → 5.84 s (17×); the referee is now about 95% G-Turing (41–54 s).
Verdict: **WIN** [NUM] on these 4 windows. Soundness rests on the stated error bound, the libm ulp assumption and Gabcke's bound, so it should be refereed before adoption.

## mpeval audit: independent rigor check of rsbatch (2026-11-01, Nathan A. Meeks & Grok) [NUM], candidate, unrefereed; no RH/Λ claims
Files: refine/mpeval_audit/ (rsc.h, rsbatch_fix.c, sc.c, sct.c, az.c, cmp.py, replay/).
| check | finding | status |
|---|---|---|
| 1 inner-sum bound | t0·log n is done in Arb (192 bits) inside c_n=n^{-1/2-it0}; double only sees v·log n with \|v\|≤33 in use (≤1000 allowed), so the phase is ≤1.4e4 rad. No 1e13-rad catastrophic reduction. Phase error ln·(ev+3\|v\|u) covers ln rounding and the product. Coefficient error comes from Arb mags. The blocked-sum γ_{B+nb+4} uses Ab ≥ Σ\|terms\|. | correct |
| 2 Gabcke | Gabcke 1979, Satz 4.2.3 (2015 reprint, eq. 8): \|R0(t)\| < 0.127 t^{-3/4} for t ≥ 200, with only the C0 term, written (−1)^{N−1}(t/2π)^{-1/4}C0. The code uses 0.127 τ^{-3/4}, which is (2π)^{3/4} ≈ 3.97× looser. C0=Ψ(p) and the sign are right. | correct but loose (could be tightened ~4×) |
| 3 sin/cos | glibc 2.41: max 0.515 ulp (≤0.504·2^-53 abs) on 3×1e6 arguments (\|x\|≤3.2, 500, 14000) vs Arb. The 4-ulp budget holds empirically, but it is an assumption, not a proof. | **fix applied** |
| fix | rsc.h: Cody–Waite reduction with fma plus Taylor deg 19/20, a priori error ≤64·2^-53 (derivation in header). It refuses \|x\|>2^20 or \|r\|>0.8, which falls back to Arb. Measured max 1.40·2^-53 on 2e6 arguments. Compiled with -ffp-contract=off. | rigorous |
| adversarial | Batch vs Arb hardy_z at 2000 random points each in 10202 [7511,7577] and 38498 [10368,10400], plus 200 points with \|v\|≤990: 0 violations, 0 sign mismatches, 0 refused. Worst err/rad is 0.0133 (orig) and 0.0129 (fix); max err 2.7e-11, rad ≤3.8e-9. | pass |
| 4-window replay (fix) | 10202 282=282, 8758 188=188, 27053 180=180, 38498 137=137, all RESOLVED, 0 fallbacks. G-RS 2.23+1.25+0.97+1.04 = 5.49 s (was 5.84 s), load ~13. | pass |
Verdict: **SOUND-WITH-FIX** [NUM]. The only gap was the unproved libm ulp assumption, now replaced by the rigorous rsc() at no speed cost. Rigor still assumes IEEE round-to-nearest double and a correct Arb/FLINT. Optional: tighten Gabcke to 0.127 t^{-3/4}.

## Grid amortization across M windows (2026-11-01) — Nathan A. Meeks & Grok [NUM]
Setup: one fftzi run of the unified-v2 copy over M consecutive ledger windows (δ=0.24, P = next power of 2 covering the span, its own ahcache copy). 3 heights (ledger rows 1081, 35151 and 59094; 64 consecutive windows each), single core (taskset 7), box load 11–14. Code and data: refine/amort/ (amort.py, amort.jsonl, stats.txt). The referee is untouched and gets no grid values.
| M | grid s/window | grid µs/zero | certify s/window | total s/window | peak RSS MB | speedup vs M=1 | spans PROVEN | per-window counts = ledger |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.181 | 4.11 | 0.211 | 0.392 | 224 | 1.00 | 192/192 | 192/192 |
| 2 | 0.199 | 4.51 | 0.269 | 0.468 | 224 | 0.84 | 93/96 | 127/192 |
| 4 | 0.205 | 4.65 | 0.329 | 0.533 | 245 | 0.74 | 45/48 | 114/192 |
| 8 | 0.202 | 4.58 | 0.499 | 0.700 | 375 | 0.56 | 14/17* | 78/136 |
| 16 | 0.293 | 6.65 | 0.777 | 1.069 | 712 | 0.37 | 6/8* | 75/128 |
(*the run was still finishing when tabulated.) Every PROVEN span had found = Turing = ledger. The per-window mismatches are not wrong counts: a certified bracket straddles a ledger window boundary (the windows overlap by 40 units), so the count is ambiguous there. Non-PROVEN spans would need the escalation ladder.
**Verdict: KILL.** The grid cost per window does not fall, because the NUFFT costs ∝ span and the per-n setup is not the dominant term. The certify stages cost more per window as the span grows. Longer spans fail attempt 1 more often, and memory grows about 3× by M=16. Keep M=1.
