# Gabcke remainder constant: sources and form check [NUM]

**Status (2026-10-18):** second published sources found for the numeric constant. The code previously used a looser, mis-typed form; it now matches Gabcke.

## Bound used in this repo
For t ≥ 200, after the C0 term of the Riemann–Siegel expansion of Z(t),

\[
|R_0(t)| \le 0.127\, t^{-3/4}.
\]

Primary source: W. Gabcke, *Neue Herleitung und explizite Restabschätzung der Riemann-Siegel-Formel*, Dissertation, Göttingen, 1979 (re-typeset PDF: doi:10.53846/goediss-5113). Gabcke's own preface states |R₀(t)| < 0.127 t^{-3/4} for t ≥ 200.

## Independent secondary sources (same numeric constant and t^{-3/4} form)
1. X. Gourdon and P. Sebah, *Numerical evaluation of the Riemann Zeta-function* (2003), §1.3: for t ≥ 200, |R₀(t)| ≤ 0.127 t^{-3/4}, |R₁(t)| ≤ 0.053 t^{-5/4}, |R₂(t)| ≤ 0.011 t^{-7/4}, citing Gabcke via Borwein–Bradley–Crandall.
   http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf
2. D. J. Platt and T. S. Trudgian, *An explicit van der Corput estimate for ζ(1/2+it)*, arXiv:1507.01261, Lemma 2.1: remainder piece 0.127 t^{-3/4} for t ≥ 200 (alongside a C0-size term 1.48 t^{-1/4}).
3. J. Arias de Reyna, *High precision computation of Riemann's zeta function by the Riemann-Siegel formula, I*, Math. Comp. 80 (2011), 995–1009 (doi:10.1090/S0025-5718-2010-02426-3): cites Gabcke (1979) for exact bounds on the first ten remainder terms on the critical line, and derives independent sharp bounds for general s.
4. D. J. Platt and T. S. Trudgian, *An improved approximate formula for the height of the Riemann zeta-function*, J. Number Theory 147 (2015), 842–851: for t ≥ 200 they use "explicit error bounds due to Gabcke" for interval Riemann–Siegel evaluations.

5. A. Kuznetsov, *Series expansions for the Riemann zeta function*, arXiv:2312.03261v2 (2024), discussion after eq. (1): with the main sum plus the first correction term (order t^{-1/4}) included, Gabcke proved the error is less than 0.127 t^{-3/4} for all t ≥ 200. This matches how the code applies the bound: `c0term` already carries the (t/2π)^{-1/4} factor, so 0.127 t^{-3/4} is added directly to Z.

## Re-check after the form correction (2026-10-18)
Window T = 3,000,000,010,320 (closepair/tz build, same arguments as `results/closepair/tz/w_c1.json`): still PROVEN, found 44,169 = Turing upper count, same sign changes (44,873) and suspect gaps (1,809). The Gabcke term Rsup dropped from 2.21e-10 to 5.57e-11 (factor ≈ 3.97 = (2π)^{3/4}), as expected. Since the corrected bound is strictly smaller, every earlier PROVEN window stays valid.

## Form correction (2026-10-18)
Earlier code and some docs wrote |R| ≤ 0.127 (t/2π)^{-3/4}. That is larger than Gabcke's published bound by a factor (2π)^{3/4} ≈ 3.97, so it was still *sound* (strictly more conservative) but not the constant as stated. All `gabcke` / `gab` helpers now use 0.127 t^{-3/4}.

Do **not** cite the RH-to-3×10¹² paper (Platt–Trudgian, Bull. Lond. Math. Soc. 53 (2021), arXiv:2004.09765) as the source of the 0.127 figure: it uses Gabcke-style RS machinery but does not restate that constant.
