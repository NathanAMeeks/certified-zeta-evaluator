# certified-zeta-evaluator

> **⚠ Candidate, unrefereed research code.** It makes **no claim about the Riemann Hypothesis** and no claim about any de Bruijn–Newman-type bound. All results are labeled [NUM] (numerical, interval-certified modulo the stated published constants) or [HEUR] (heuristic extrapolation). Independent review is welcome.

Author: Nathan A. Meeks. MIT License.

## What it is
This is an Arb-certified, FFT-batched evaluator for the Riemann–Siegel Z function. It also includes a Turing/Trudgian completeness check, which counts the sign changes in a window and proves that no zeros were missed.

- **Batching:** the main sum Σ n^{-1/2} e^{i(t ln n)} is evaluated on a grid of P points through L Taylor-NUFFT passes plus an FFT (`python/fftb.py`, float plus a hand bound; `python/fftarb.py`, full Arb).
- **Error terms:** the θ phase uses its own certified Taylor expansion and remainder. The Riemann–Siegel C0 term and Gabcke's remainder bound are evaluated in Arb.
- **Signs:** a point's sign counts only if its Arb ball excludes 0.
- **Completeness:** Turing's method with Trudgian's bound on ∫S(t) brackets N(a) and N(b). The window is *proven* when the number of certified sign changes equals the maximum allowed.
- **C/FLINT port:** `c/fftz.c` offers two DFT modes. One is `acb_dft`. The other is a double radix-2 FFT with a rigorous Higham (Thm 24.2) forward-error enclosure, added as a ball radius.
- **Earlier version:** `python/cert2.py` is the Taylor-reuse evaluator (all Arb). `python/tune.py` is its auto step-size tuner.

## Unproven-constant dependencies
Rigor depends on these published bounds, which are used as stated and not re-proven here:
- Gabcke (1979), the Riemann–Siegel remainder after C0: |R_0(t)| ≤ 0.127 t^{-3/4} for t ≥ 200 (restated by Kuznetsov, arXiv:2312.03261; same form in Gourdon–Sebah 2003 and Platt–Trudgian arXiv:1507.01261). The code uses exactly this bound (form corrected 2026-10-18 from the looser 0.127 (t/2π)^{-3/4}); see `docs/GABCKE_SOURCES.md` and `docs/GABCKE_CHECK.md`.
- Trudgian (2011), the S(t) integral: |∫S| ≤ 2.067 + 0.059 log t.

## Build and run
Python needs `python-flint`, `numpy` and `mpmath`.
```
cd python
python3 fftarb_run.py e      # 1e13 window [60,460], 131,072 points -> JSON summary
```
For C, install FLINT 3.x:
```
c/build.sh                   # FLINT_ROOT=/prefix c/build.sh for a custom install
./c/fftz T delta P a b Lw [signfile] [dbl]
python3 c/rt.py ./c/fftz ... # wall time + peak RSS
```

## Measured results [NUM]
| window | proven zeros | N(a) |
|---|---|---|
| 2.68e11 + [60,460] | 1,558 | 10^12 + 226 |
| 3e12 + [60,460] | 1,712 | 12,362,403,019,918 |
| 1e13 + [60,460] | 1,789 | 43,124,192,297,371 |
| 1e13, 2^20 points (C, dbl) | 44,201 | — |

The C and Python versions give identical proven counts. Across 131,072 points there are 0 sign disagreements in both C modes.

Speed:
- C dbl takes about 0.31 ms per zero at 3e12 and about 0.35 ms per zero at 1e13, on a single core.
- Memory use is about 740 B per point.
- Throughput scales linearly up to 4 cores.
- Full Arb FFT batching is about 2,300–7,600× faster than per-point Arb direct, and about 28× faster than cert2 at 2.68e11.

The float-vs-Arb referee check found 0 bound violations. Details are in `docs/C_PORT.md`, `docs/FFT_BATCH.md` and `results/`. Cost extrapolations in `docs/PUSH_HEIGHT.md` and `results/*push_table.txt` are [HEUR].

## Close-pair finder [NUM][LEM] (`c/closepair/fftzcp.c`, `docs/CLOSE_PAIR.md`)
When the certified sign-change count falls short of Turing's maximum, the finder adds certified evaluations in targeted gaps (deficit-targeted refinement, in the same pass) and loops until the window is proven or its budget runs out.
- **Lemma CP (soundness):** a window is reported *proven* only when the certified sign changes equal Turing's count exactly. Where refinement looks is heuristic and affects only speed and whether a window finishes proven; it cannot produce a false "proven". With zero refinement budget, windows are correctly reported unproven.
- **Lemma TZ (zero-free gap test, `c/closepair/tz/fftztz.c`):** a proven Taylor/linear-interpolation bound on Z'' certifies ~99.8% of grid gaps as zero-free, so refinement only looks at suspect gaps. Across 10 windows all 15 missed pairs were in suspect gaps; 0 violations. It does not prove uniqueness inside a sign-change gap (no certified bound on R'); Turing's count still supplies that.
- All 10 test windows proven, each needing at most one refinement round, 0 sign disagreements between coarse and fine grids.
- Earlier in-process refinement prototypes: `c/inproc/` (`docs/REFINE.md`).

## Speed results [NUM] (`c/speed/fftzsp.c`, `docs/SPEED.md`)
Main sums in double precision with rigorous (hand-derived) error bounds, lower Arb precision, and ~60× cheaper refinement:
- **0.090 ms/zero/core at 3e12** (was 0.519) and **0.107 ms/zero/core at 1e13** (was 0.820), 2.4 GHz Xeon.
- Contiguous proven chains: **132,506 zeros near 3e12** and **138,439 zeros near 1e13**, with Turing endpoint counts matching across windows.
- Best spacing ~5 evaluations per zero; wider spacing makes windows fail as unproven, never wrongly proven.

**Crude per-GHz comparison [HEUR]:** Platt–Trudgian (arXiv:2004.09765) ≈ 7.9 ms·GHz per zero; Gourdon (2004, not interval-certified) ≈ 0.011; ours ≈ 0.22–0.26. That suggests ~30× fewer cycles per zero than Platt–Trudgian and ~20× more than Gourdon. Caveats: their figures average a full production run over all heights, ours are short top-height windows; different CPU generations make GHz normalization crude; ours is not a production run.

## Band-limited interpolation, v7 [NUM][LEM] (`c/interp/fftzi_v7.c`, `docs/SPEED2.md`)
- Certified Gaussian-windowed sinc interpolation of the main sum (Lemma BL, own derivation, needs a referee) plus extremum targeting and a certified midpoint sign-flip probe for close pairs.
- Measured: two contiguous chains above 3x10^12 and 10^13 (132,506 and 138,439 zeros) at 0.0050-0.0055 ms/zero; all 10 test windows PROVEN at 0.0069-0.0079 ms/zero, with 0 sign disagreements. Certification is unchanged: certified sign changes must equal Turing's count.
- Credit: the certified midpoint sign-flip probe is Nathan A. Meeks' idea, first recorded in "Phase-Coherent Holographic Compression" (X Article, Feb 18, 2026, https://x.com/NathanAMeeks/status/2024058294244180131).
- Build with `FLINT_ROOT=/prefix c/interp/build.sh`. Usage is in `docs/SPEED2.md`.

## Prior art
FFT-batched Riemann–Siegel evaluation is the Odlyzko–Schönhage method (1988). Gourdon (2004) and Platt (with Trudgian, 2021, fully interval-certified) used versions of it for record verifications. **This repository is an independent certified implementation, not a new method of proof.** Its contributions are a fresh reproducible implementation, the close-pair finder, and the measured timings.

## Open items / known limits
- **Referee needed:** the Higham FFT bound application and the hand-derived double-precision error bounds coded in `c/speed/fftzsp.c`. An empirical Arb-vs-double check found 0 violations over about 2.37M grid points, with observed error at most about 0.13% of the bound (`docs/DP_CHECK.md`); that is a consistency check, not a proof.
- **Gabcke (1979) constant:** second sources confirmed (Kuznetsov arXiv:2312.03261; Gourdon–Sebah 2003; Platt–Trudgian arXiv:1507.01261; Arias de Reyna 2011); a citation check, not a re-proof. Form corrected to 0.127 t^{-3/4} on 2026-10-18 in every evaluator, including v7/v8 (`docs/GABCKE_SOURCES.md`, `docs/GABCKE_CHECK.md`).
- **Trudgian (2011) S(t) integral:** constants 2.067 + 0.059 log t₂ cross-checked against Math. Comp. 80 Thm 2.2 on 2026-10-18 (`docs/TRUDGIAN_SOURCES.md`).
- **No independent reimplementation** has checked these results.
- The C code has no per-point Arb-direct fallback (never needed in these runs).
- `docs/REFEREE_FFT.md` lists 5 flaws in the earlier float bound (`python/fftb.py`), fixed in `fftarb.py` and the C code; treat `fftb.py` as superseded.

Build the new programs with `FLINT_ROOT=/prefix c/speed/build.sh` (likewise `c/closepair/build.sh`, `c/closepair/tz/build.sh`, `c/inproc/build.sh`). Run scripts are in `scripts/` (they expect the binary in the working directory).

## Layout
`python/` evaluators · `c/` C/FLINT port · `docs/` method write-ups and referee notes · `results/` JSON outputs and cost tables · `scripts/` run/analysis scripts.
