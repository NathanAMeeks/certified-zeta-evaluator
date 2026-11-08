# rh_extend — contiguous certified zeta-zero verification above 3e12 (v9)
[NUM] candidate, unrefereed. Credit: Nathan A. Meeks & Grok (v9 = v7 speed + v8 gated self-correcting loop, directed by Nathan).
This extends a *verified height* (all zeros in [3e12, covered_hi] are simple and on the critical line, assuming Platt–Trudgian below 3e12).
It does **not** prove RH and must never be described as doing so.

## Run
- `./run_batch.sh` — no args, one stdout line, safe for a daily routine. Exit codes: 0 OK, 2 window failed, 3 halted (open FAILED.md), 4 skipped (lock held), 5 driver crash.
- `./status.sh` — one line: height covered, zeros since 3e12, windows, failures, last run.
- Batch size: `BATCH_WINDOWS` env (default 1300 ≈ 15 min wall at ~0.68 s/window incl. ladder escalations on one pinned core, CPU 7).

## Windows
- Evaluator (since 2026-10-29, window 8758 resolution): `/workspace/rh_extend/v9p/v9/v9run.sh T P b` (previously cze_repo copy) (fixes merged as local commit 8ebc65b; binary byte-identical to v9p build; v9p kept as fallback) — a copy of cze_repo v9 (3c24e7f) whose only change is the C0 removable-singularity fix (C0_FIX.md, [LEM] candidate); cze_repo is untouched. Rebuild: `sh v9p/build.sh`. Original driver behaviour ( runs `fftzi_v9 T delta P 80 b 60 /dev/null dbl` with the v9 retry ladder:
  delta 0.24 at caller P, then P doubled, then delta 0.2/0.16/0.12). If the ladder fails, the driver reruns the whole ladder once with `MAXR=200` (refine-round cap, an existing fftzi_v9 env knob; default 20) — search-only, same Lemma CP acceptance. Brackets come from the evaluator's `ZDUMP` env hook (no source change).
- P=65536, b=10400 (SPEED2 v9 benchmark config). Window = [T+80, T+10400], width 10320; step 10280 → overlap 40. Windows idx 0–35 used P=131072, b=31040 (step 30920); idx 36 failed there (UNPROVEN_budget after the full ladder) but the same span proved with three 65536/10400 windows, so the driver was switched (see logs/FAILED_resolved_2026-10-28.md).
- First window T=3e12−80, so coverage starts exactly at t=3e12. No earlier contiguous ledger above 3e12 existed (PUSH_HEIGHT.md / SPEED2 only had isolated windows), so it starts fresh.
- AHCACHE = private copy `ahcache.bin` of `meeks_test/refine/speed2/ahcache_5e13.bin` (covers ≤5.2e13).

## Lemma CP (PROVEN rule)
A window is PROVEN iff v9 reports `status=PROVEN`, the certified sign-change count `found` equals the Turing count `Nb_upper − Na_lower`,
and the number of bracket lines equals `found`. Then Na, Nb are exact values of N(lo), N(hi).
Chain check per window: `overlap_N = N(prev hi) − N(this lo) ≥ 0` and equals the brackets of this window inside the overlap (±1 for a bracket straddling the edge).
Total zeros since 3e12 = N(covered_hi) − N(3e12) (no double counting of overlaps).

## Files
- `checkpoint.json` — next window (T arg, index), N(3e12), N(covered_hi), covered_hi, totals; replaced atomically only after a PROVEN window.
- `ledger.csv` — one row per attempted window (idx, T_arg, lo, hi, zeros_turing, cert_sign_changes, Na, Nb, status, delta, P_used, attempts, sec_v9, wall_s, overlap_N, overlap_brackets, timestamp).
- `brackets/wNNNNNNN_T<T>.txt.xz` — certified sign-change brackets, lines `u_lo u_hi` with t = T + u (≈80 KB/window, ≈100 MB per 1300-window batch).
- `logs/windows.jsonl` (raw v9 JSON per window), `logs/wNNNNNNN.err` (only when v9 escalated/failed), `logs/runs.log` (summary lines), `logs/driver.err`.
- `FAILED.md` — open failures; `run.lock` — flock so runs cannot overlap.

## After a failure
The failed window is logged in ledger.csv + FAILED.md, checkpoint stays on it, and every later run halts (exit 3) until resolved — no window is ever skipped.
To resume: inspect `logs/wNNNNNNN.err` and the jsonl row, fix the cause (e.g. try larger P / finer delta by hand for that same T and confirm PROVEN),
then `mv FAILED.md logs/FAILED_resolved_<date>.md` and rerun `./run_batch.sh`; it retries the same window first.
A crash between a ledger append and the checkpoint save can leave a duplicate ledger row for one window; checkpoint.json is authoritative.

## C0 fix (Oct 28)
Window 308 (t≈3000003916629, p=3/4) failed because v9's c0term refused |cos 2πp|<1e-3. v9p evaluates Ψ near p=1/4,3/4 with an exact sinc rewrite (C0_FIX.md). Validation: v9p/val_psi.py; bench in scratch/bench_v9p.json.

## N-step fix (Oct 28)
Window 519 crossed the Riemann–Siegel N step t_N=2π(N+1)² and v9 (fixed N per window) refused points past it. v9p's c0term now switches to N'=N+k past t_N and adds the missing main-sum terms (NSTEP_FIX.md). The driver pre-checks each window for a t_N crossing, logs it in logs/nstep.log, and skips the slow MAXR rung for such windows.

## Gates (Oct 30)
driver.py now runs two HEUR gates (gates.py, GATES.md): C0 0/0 Dedekind-cut logging, and the close-pair centre-out chunked-reuse resolver (Nathan A. Meeks). Ledger gained a `notes` column (rows ≥ 10202).
