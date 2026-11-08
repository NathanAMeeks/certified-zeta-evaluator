# Preparation tests (2026-11-07 MST)

Nathan A. Meeks & Grok. Candidate/unrefereed; no RH claim.

- **C0_regression**: 182630 mpmath comparison checks, 0 violations, 0 not_ok (not rigorous independent certification)
- **N_step_regression**: PASS; straddling ball refused; both sides accepted
- **evaluator_smoke**: PROVEN, 132508 signs/count match, Na/Nb match first committed window
- **rebuilt_binary**: byte-identical to live evaluator SHA256 7c48db1452fae772b56b9dcdb666d6fc93accb5524d178d85a6bb7448ceb182d
- **independent_Arb**: first bracket signs +1/-1 only; not whole-window gate
- **raw_sample**: 65 stratified files: counts, sorted/non-overlap, 27 available historical raw hashes PASS
- **not_run**: ['full database endpoint sign/Turing re-certification', 'full referee-window0 Arb command', 'fast-tier portability rebuild/independent audit']
- **compiler_warnings**: GCC allocation-size-range warning for source P allocation; build succeeded; no source change made
- **production_safety**: No live files/locks changed, no run started/restarted; original cze_repo remains clean/read-only.
- **Full bulk artifact integrity**: all 19 chunk checksums passed; all 165,476 tar entries read and hashed; all 165,468 compressed bracket files matched the source manifest. No symlinks, absolute paths, or failed-window bracket included. This is not a decompressed sign/Turing proof rerun.
