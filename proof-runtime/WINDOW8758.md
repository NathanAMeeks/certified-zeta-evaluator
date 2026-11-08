# Window 8758 resolution — [NUM] candidate, unrefereed (not RH)
Method: Nathan A. Meeks (start at the centre of the suspect span, test adjacent points on an adaptive fine mesh expanding outward, shrink step near small |Z|, close-pair finder). Execution: Nathan A. Meeks & Grok, 2026-10-29 (MST). Evaluator: v9p only (scratch/w8758/); cze_repo untouched.

- Failure: idx 8758, T=3000090775200, window [3000090775280, 3000090785600]: found 44169 vs Turing 44171 (UNPROVEN_budget after ladder + MAXR=200).
- Centre-out fine mesh: subwindow T=3000090783500 [..783580, ..784100], delta 0.05 and 0.02: PROVEN (2225 = Turing). Closest pairs nearest centre (784.87, 848.38) refined at delta 5e-4: both already resolved by the delta 0.24 run — not the missing pair.
- Close-pair finder (diff of delta-0.24 brackets vs fine brackets): only discrepancy is coarse bracket u∈(8711.544, 8711.76) holding 3 zeros. Missing pair, refined at delta 5e-4:
  zero 1 in [3000090783911.6323, 3000090783911.6328], zero 2 in [3000090783911.6328, 3000090783911.6333]; gap < 0.001 (~0.0005; ≈0.2% of mean spacing ≈0.234). Third zero in [..911.7148, ..911.7153].
  NOTE: the pair lies ~31.6 above the reported suspect span end (783880); it is inside window 8758.
- Full window rerun: fftzi_v9 3000090775200 0.05 262144 80 10400 60 /dev/null dbl, MAXR=400 → status PROVEN, found 44171 = Turing (Nb 12362791578249 − Na 12362791534078) = brackets 44171; overlap with 8757: ovN=ovB=173. Lemma CP holds → PROVEN [NUM].
- Recorded: ledger.csv row (new PROVEN row appended after the UNPROVEN row), brackets/w0008758_T3000090775200.txt.xz, checkpoint advanced. FAILED.md → logs/FAILED_8758_resolved.md. Driver now points at v9p/v9/v9run.sh.
