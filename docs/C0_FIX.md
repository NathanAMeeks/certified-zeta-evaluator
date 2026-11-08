# C0 removable-singularity fix in v9p — [LEM] candidate, unrefereed
Grok for Nathan A. Meeks, Oct 28 2026. Patched copy: c/interp/fftzi_v9.c (function `psi_c0`, used by `c0term`). /workspace/cze_repo untouched.
Not RH; this only removes an evaluator blind spot for verifying zeros in windows.

**Problem.** C0(t) = (−1)^{N−1} (t/2π)^{−1/4} Ψ(p), Ψ(p)=cos(2π(p²−p−1/16))/cos(2πp), p = √(t/2π) − N. v9 returned "not certified"
when |cos 2πp|<1e-3 (|p−p0| ≲ 1.6e-4 for p0 ∈ {1/4,3/4}; ≈ ±1.4e3 in t at 3e12), so windows containing such t could never be PROVEN.

**Lemma C0FIX.** Let p0 ∈ {1/4, 3/4}, s = −1 for p0=1/4, s=+1 for p0=3/4, h = p − p0, 0<|h|<1/2. Then
  Ψ(p) = sin(πh(1+2sh)) / sin(2πh) = ((1+2sh)/2) · sinc(πh(1+2sh)) / sinc(2πh),  sinc x = sin x / x,
and the right side is analytic for |h|<1/2 with value 1/2 at h=0 (the continuous extension of Ψ).
*Proof.* p=3/4+h: p²−p−1/16 = −1/4 + h/2 + h², so cos(2π(·)) = cos(−π/2 + πh + 2πh²) = sin(πh(1+2h)); cos(2π(3/4+h)) = cos(3π/2+2πh) = sin(2πh).
p=1/4+h: p²−p−1/16 = −1/4 − h/2 + h², numerator = sin(−πh + 2πh²) = −sin(πh(1−2h)); denominator cos(π/2+2πh) = −sin(2πh).
Divide numerator and denominator by 2πh to get the sinc form. sinc(2πh) has no zero for |2πh|<π, so the quotient is analytic on |h|<1/2. ∎
No cancellation: both sincs are ≈1, and the h factors cancel symbolically, not numerically.

**Enclosure.** For |mid(p) − p0| < δ = 0.01 the code computes h = p − p0 (p0 is an exact dyadic), the two sincs with Arb's `arb_sinc` (rigorous for balls
that contain 0), and the quotient as Arb balls. It marks the result uncertified only if the sinc(2πh) ball contains 0 (impossible: sinc(2πh) ≥ sinc(π/50) > 0.99 for |h| ≤ 0.01
plus any radius < 0.4). Outside δ, it uses the original quotient unchanged. There |cos 2πp| ≥ sin(0.02π) > 0.06, so the old 1e-3 cutoff never applies, and the δ=0.01
region covers the old cutoff region (|h| ≲ 1.6e-4) about 60 times over. The branch choice (by midpoint) only selects between two exact identities, so it is rigorous either way.
No Taylor remainder is needed. Env `C0OLD=1` restores the old formula, for A/B tests only.

**Dedekind-cut enclosure (Nathan A. Meeks' idea).** Following Nathan A. Meeks's suggestion to treat this as a Dedekind cut, Ψ over an input ball is returned
as a cut L ≤ Ψ ≤ U where both bounds are proven, built only from the sinc form above (the raw 0/0 quotient is never formed within δ). `psi_cut` bisects the h-ball
on both sides (it straddles p0 when the ball contains it), evaluates each piece in Arb, and takes L = min of the lower bounds and U = max of the upper bounds. It keeps splitting while that
tightens the hull by more than 10% and the width is above C0TGT (1e-15), up to depth C0MAXD=10.
Cut width vs CP margin: point inputs give a width of 1.0e-28 (exactly 0 at p=p0, value 1/2). For the evaluator's block balls (H=8 in t means p-radius ≈ 9.2e-7 at 3e12)
the cut width is 1.0009·2r, i.e. the true variation of Ψ over the ball (|Ψ'(p0)|=1). Bisection can't go below that, and it doesn't need to: C0r ≈ 1.8e-6/2 · (t/2π)^(-1/4) ≈ 1.1e-9,
against the dbl-path acceptance C0r < 1e-3 and against EPH ≈ 9.8e-8 in the point budget Ent. So C0 adds about 1% of the per-point error, which is well inside what Lemma CP needs.

**Validation [NUM].** c/interp/v9/c0check/val_psi.py (build t_psi.c there): 182,630 checks (dense p-grids ±0.012 and ±2e-4 around 1/4 and 3/4, exactly p0, h=±1e-16…1e-3, δ·(1±1e-12), plus 999 points on (0,1);
prec 96 and 256; point balls and radius 1e-12 and 1e-9 with both ball ends checked) against mpmath at 110 digits (original quotient; limit at p0): **0 violations**,
0 uncertified. Max cut width within δ at prec 96 (point inputs): 1.0e-28 near both p0 (re-run after the cut change: still 0 violations).
Windows: 308 (contains p=3/4, t≈3000003916629.1) and T=3000008253240 (contains p=1/4, t≈3000008258241.0) both PROVEN on the first attempt with all grid points certified.
Benchmarks (10 bench + 2 chain + 2e13+20800, ×2 reps): v9p 26/26 PROVEN, identical counts to v9, wall 13.28 s vs 13.26 s (sinc version), 12.90 s vs 12.94 s (cut version): no speed change.
