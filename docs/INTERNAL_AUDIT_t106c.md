# Audit re-audit t106c — 2026-11-02 15:00 MST (read-only; scratch runs in /workspace/audit_scratch/t106c)
VERDICT: PASS for live (XZLEVEL=1 V9DIRECT=1 GATE_PROCS=7 GRS_TO_MIN=30 HYB_TMESH=1 HYB_VC=0, PC_FAST=0). 3 non-blocking nits below.

## 1 Byte identity (own check vs t106b/bench/t103_v2)
Code's runs final400_to30, S1_pcfast_to30 AND my reruns off/on: verdicts 400/400 (19 IND-gated), cols1-12 sha 43cc5aed, raw brackets (sorted concat, xz -dc) 701ccaa5, zeros_turing sum 17,667,837. All identical.

## 2 Soundness
- xz -1: only the .xz container changes; brackets_sha hashes raw bytes; CRC64 check kept; old -6 files decompress. Sound. Default in driver is now '1'.
- GATE_PROCS=7: scheduling only (Pool size), identical output. Sound.
- HYB_TMESH=1 search-only; HYB_VC forced '0' at driver line 11. Sound.
- Consecutive-timeout rule (indgates._rsbatch1): count++ on timeout, reset to 0 on any completed call, drop at GRS_TO_MAX=2, probe every 20 groups, killpg, points stay 0 -> Arb; per-batch reset kept. Matches spec. tests/timeout/test_consec.py exists.
- GRS_TO_MIN=30: timeout only affects speed (timed-out points -> Arb), never verdict. Still caps a hung child. OK. NIT a: code default is still '5'; live relies on env. Set default 30 in indgates.py so a missed env var can't bring back the Arb-fallback slowdown.

## 3 S1 (PC_FAST)
Off by default (gates.py:160). When on, close-pair signs use rs_fast512/rs_fast and that binary is removed from the G-RS tier list (indgates ~216), so the referee never uses the same evaluator. Sign 0 -> Arb, t<=200 -> Arb. Independence kept at the binary level. NIT b: both sides are still Riemann-Siegel family (rs_fast vs rs_tb_fix/rsbatch_t100); before S1 the close-pair used Arb. Fine as a certified evaluator, but say so in the referee notes. _TIERS is fixed at import, so PC_FAST has to be in the env before the process starts (it is, via env VAR=).

## 4 WIN_LEN / deferrals
10400 is supported as "no change" (lowest CPU/zero, lowest redo, lowest RSS). NIT c: the sweep ran with GRS_TO_MIN=5, so 2-9 timeouts per length added Arb CPU noise (for example, 31000 had 2 timeouts and the second-best CPU). The ranking isn't clean. If longer windows still matter, rerun 10400 and 31000 at TO_MIN=30. That's not blocking because the default doesn't change.
Deferrals hold: S2 has 0 retries at every length, so the payoff is about 0 and a certified splice needs a proof. S3 needs a new C build, a Z'' bound and a mismatch proof. S5 needs an audit of the statics first. V9DIRECT already gets the wrapper savings.

## 5 Rerun (load 0.9-6.3)
off 46.4 s wall / 272.7 cpu-s; on 44.6 s / 278.5 cpu-s; 0 GRS timeouts. Confirms Code's 46.1 and 42.5 within noise. S1 gives about 4-8% wall and 0 CPU, so it's marginal and fine to leave off.
Speed: about 0.115 s/window wall and about 27 z/ms/GHz/core CPU. xz -1 halves CPU. GRS_TO_MIN=30 is what avoids the 3.5x Arb-fallback wall cost under load.
