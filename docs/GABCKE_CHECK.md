# Referee note: Gabcke's constant

Primary source: W. Gabcke, "Neue Herleitung und explizite Restabschätzung der Riemann-Siegel-Formel", Dissertation, Göttingen 1979, doi:10.53846/goediss-5113. Formula (6) writes Z(t) = 2 Σ_{n≤N} n^{-1/2} cos(θ(t) − t log n) + (−1)^{N−1} (2π/t)^{1/4} Σ_{k≤K} C_k(z) (2π/t)^{k/2} + R_K(t), with R_K(t) = O(t^{-(2K+3)/4}). Formula (8): |R_0(t)| < 0.127 t^{-3/4} for t ≥ 200 (optimal for K ≤ 4).

Second source restating it: A. Kuznetsov, "Series expansions for the Riemann zeta function", arXiv:2312.03261 (2023), p.1-2, with the same normalisation (C0 term with (2π/t)^{1/4} factor, error < 0.127 t^{-3/4}, t ≥ 200).

Finding 1 (soundness; RESOLVED 2026-10-18, applied to all evaluators incl. c/interp v7/v8 on 2026-10-28: code now uses 0.127 t^{-3/4}, see docs/GABCKE_SOURCES.md): the earlier code (gabcke() in c/*.c, gab() in python/cert2.py) uses 0.127 (t/2π)^{-3/4} = 0.127 (2π)^{3/4} t^{-3/4} ≈ 0.504 t^{-3/4}. That is about 3.97× LARGER than Gabcke's bound, so every enclosure stays valid (conservative). Using the exact form 0.127 t^{-3/4} would tighten the margin ~4× (not required).

Finding 2 (citation): docs/CERT_TAYLOR.md says the constant is "restated in Platt–Trudgian, Math. Comp. 90 (2021)". The Platt–Trudgian paper is Bull. London Math. Soc. 53(3) (2021) 792–797, doi:10.1112/blms.12460, and its text (arXiv:2004.09765) does not mention Gabcke or 0.127. README.md and docs now quote the bound as 0.127 t^{-3/4}, which is what the code uses.

Caveat: this is a citation check (two texts agree on the statement), not an independent re-proof of Gabcke's bound.
