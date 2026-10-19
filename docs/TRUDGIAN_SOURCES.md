# Trudgian S(t)-integral constants: source check [NUM]

**Status (2026-10-18):** constants in code match the published theorem.

## Bound used in code

For `t2 > t1 > 168π`,

\[
\Bigl|\int_{t1}^{t2} S(t)\,dt\Bigr| \le 2.067 + 0.059\log t_2
\]

with `log` the natural logarithm (as in Math. Comp. usage; coded via `arb_log` / Python `arb.log`).

Sites: `python/fftarb.py` (`C=lambda t2: arb(2.067)+arb(0.059)*(Ta+arb(t2)).log()`), and the C Turing bracket in `c/fftz.c`, `c/speed/fftzsp.c`, `c/inproc/fftzr.c`, `c/inproc/fftzr8.c`, `c/closepair/fftzcp.c`, `c/closepair/tz/fftztz.c` (`arb_mul_d(..., 0.059)` then `arb_add_d(..., 2.067)`).

## Primary source

T. Trudgian, *Improvements to Turing's method*, Math. Comp. **80** (2011), 2259–2279 (doi:10.1090/S0025-5718-2011-02470-1; arXiv:0903.1885v3).

Theorem 2.2 (arXiv PDF, §2.1): If \(t_2 > t_1 > 168\pi\), then \(\bigl|\int_{t_1}^{t_2} S(t)\,dt\bigr| \le 2.067 + 0.059\log t_2\).

The paper notes these constants minimise the right-hand side of the Lehman–Brent Gram-block criterion near \(t_1 > 2\pi\cdot 10^{12}\), which is the height range of this repo's windows. A footnote records that the historical \(168\pi\) floor may itself be replaceable by a smaller number; our windows sit far above either floor.

## Cross-checks

1. arXiv:0903.1885v3 PDF text extracted 2026-10-18: Theorem 2.2 states the same pair \((a,b)=(2.067,0.059)\) and the same \(t_0=168\pi\).
2. Journal reference matches the arXiv abstract (Math. Comp. 80, 2259–2279, 2011).
3. Earlier triples in the same paper (Turing \(2.07,0.128\); Lehman \(1.7,0.114\)) are strictly weaker for large \(t\); the code does **not** use those.

## What this does / does not settle

- Settles: the numeric constants and the hypothesis \(t>168\pi\) as coded.
- Still open: independent reimplementation of the Turing bracket; referee of hand-derived double-precision FFT / phase-model bounds in `c/speed/fftzsp.c` (see `docs/SPEED.md`).
