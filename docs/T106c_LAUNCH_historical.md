# t106c go-live steps (Nathan A. Meeks & Grok; candidate, unrefereed). NOT for use until Audit PASS and Nathan's OK; Math runs the switch.
Same procedure as ../t106b/LAUNCH.md (stop at batch boundary with flock on run.lock -> loop rc=4; back up; copy; release; setsid nohup loop), with these differences.
Files to copy from scratch/t106c into live: driver_t106.py, gates.py, spotcheck/indgates.py (run.py is the bench timer wrapper, do NOT copy).
sha256: driver_t106.py e79fd0f9…1822, gates.py ffa8a81a…b99a, spotcheck/indgates.py 49c439a8…f86a (GRS_TO_MIN default 30) (write full list: `sha256sum driver_t106.py gates.py spotcheck/indgates.py > LAUNCH_SHA256.txt` and diff after copy).
Env for the live loop (add to the batch env via `env VAR=`, never export): XZLEVEL=1 V9DIRECT=1 GATE_PROCS=7 GRS_TO_MIN=30 HYB_TMESH=1 (already) HYB_VC=0 (must stay 0). PC_FAST=0 by default (optional S1, identical output, ~8% wall). WIN_LEN unset (=10400).
Pre-checks: `python3 -m py_compile` on the three files; tests/ (v9direct, timeout/test_consec, test_timeout, gt_batch_eq) all PASS; brackets_sha unchanged because it hashes uncompressed bytes (xz level only changes the .xz files; old -6 files still decompress).
First-batch checks: startup verify bad=0; no rc=3; ledger cols 1-12 format unchanged; log shows V9DIRECT and XZLEVEL=1; GRS_TIMEOUT count low.
Rollback: copy *.preT106c backups back at a batch boundary, drop the env vars.

## Referee note (Audit nit 2)
With PC_FAST=1 (S1), the close-pair signs and G-RS both use Riemann-Siegel evaluators (rs_fast512 in pair_cert; G-RS keeps a different RS evaluator and Arb). Independence is only at the implementation level, not the method level. Keep PC_FAST=0 for referee hand-offs unless the referee accepts this.
