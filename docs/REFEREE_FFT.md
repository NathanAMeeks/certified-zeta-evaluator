# Adversarial referee of fftb.py float error bound
Stress test: the float Z at ALL 131,072 grid points was compared against the fully-Arb enclosure (fftarb.py) at 1e6 (two densities) and 2.68e11. The test measured observed |Z_float − mid(Z_arb)| against the float bound minus Gabcke.
- Observed max error: 7.1e-11 / 1.6e-11 (1e6) and 1.7e-11 (2.68e11).
- Minimum slack, bound/observed: 294×, 1194×, 31,197×. Violations: 0 / 0 / 0.
Flaws found by re-deriving each term. None is violated in practice, but each one makes the float bound unproven as written:
1. Dead term: `rnd` ends with `...*0`. The x_n-rounding contribution is therefore covered only by `xerr`, which undercounts the needed 2πP·|δx| ≈ 3Pδ·lnN·u by ~1.5×. np.log's 1-ulp error is not budgeted separately.
2. Point mismatch: F is computed at the exact pδ, but θ is evaluated at the float u_p = u0+δ·p. The difference, ≈u·|u_p|·θ' (≈7e-13 at u=524, t=2.68e11, times 2S0 ≈ 6e-10), is not in the bound. It is comparable to Gabcke (1.3e-9) at that height.
3. FFT rounding: 5u·log2M·√M·‖G‖₁ assumes Cooley–Tukey with twiddles accurate to O(u). That accuracy is not documented or guaranteed for numpy/pocketfft, so the constant is unproven.
4. C0 slack 1e3·a^{-1/2}·64u(1+4πa) is ad hoc. It is not derived from the conditioning of cos(2π(p²−p−1/16))/cos(2πp).
5. S0 and the truncation tail are computed in float, without directed rounding (negligible, but unproven).
Conclusion: the float bound is empirically very conservative but NOT a proof. Superseded by the Arb path below.
