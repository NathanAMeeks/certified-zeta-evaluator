# Window 10202 resolution — [NUM] candidate, unrefereed (not RH)
Method: Nathan A. Meeks (centre-out close-pair search). Execution: Nathan A. Meeks & Grok, 2026-10-30 (MST). Evaluator: v9p only; cze_repo untouched; nothing pushed.

- Failure: idx 10202, T=3000105619520, window [3000105619600, 3000105629920]: 44167 certified vs Turing 44169 (UNPROVEN_budget after ladder + MAXR=200), 2026-10-29 22:40 MST.
- Per Nathan's instruction: NO full-window fine rerun. Chunked reuse (GATES.md): the original 44167 certified brackets kept; deficient chunk located by bisection on
  exact N at chunk boundaries (each probe a PROVEN small v9p run, Lemma CP on the chunk): probes u=5244 (deficit right), 7983 (left), 6777 (right), chunk u∈[7542,7866] deficit 2.
- Only that 324-unit chunk rerun (delta 0.02, P=32768, PROVEN, N endpoints match the probes) and spliced in place of the 7542–7866 original brackets.
- Missing pair (stepped over by the coarse mesh, in the gap between original brackets (7552.56,7552.68) and (7552.932,7553.04)); centre-out refinement at delta 5e-4:
  zero 1 in [3000105627072.8280, 3000105627072.8285], zero 2 in [3000105627072.8290, 3000105627072.8295]; gap < 0.0015 (≥0, ~0.001; ≈0.4% of mean spacing 0.234).
- Lemma CP whole window: 44167 − 14 (old, in chunk) + 16 (chunk) = 44169 = Turing (Nb 12362855111362 − Na 12362855067193); brackets sorted, non-overlapping; overlap with 10201: ovN=ovB=170 → PROVEN [NUM].
- Recorded by the gated driver: ledger row (PROVEN, notes HEUR-gate:close_pair), brackets/w0010202_T3000105619520.txt.xz, checkpoint → idx 10203. FAILED.md → logs/FAILED_10202_resolved.md. Scratch: scratch/w10202/.
