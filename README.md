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
- Gabcke (1979), the Riemann–Siegel remainder: |R| ≤ 0.127 (t/2π)^{-3/4}.
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

## Known limits
- **No close-pair refinement step.** In 2^20–2^21-point windows at 3e12, rare close pairs of zeros fell between grid points, so Turing's check (correctly) refused to certify those windows: 2–18 zeros were unaccounted for. Production use needs extra certified evaluations in the gaps Turing flags. That step is not built.
- The C port has no per-point Arb-direct fallback. These runs never needed it.
- The combine and scatter are still done in Arb. They are the current bottleneck.
- `docs/REFEREE_FFT.md` lists 5 flaws in the earlier float bound (`fftb.py`). `fftarb.py` and the C port fix them. Treat `fftb.py` as superseded.

## Layout
`python/` evaluators · `c/` C/FLINT port · `docs/` method write-ups and referee notes · `results/` JSON outputs and cost tables.
