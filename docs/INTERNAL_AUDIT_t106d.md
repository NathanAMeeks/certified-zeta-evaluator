# t106d re-check 2 (Audit worker, 2026-11-03 ~08:25 MST) — VERDICT: PASS (F2a, F2b closed; F1, F3 still PASS)
Integrity: gates.py bc6e79aa…4497, normalize_ledger.py 1d0cd20d…ca2c, backfill e4c89e16…da38 all OK; LAUNCH_SHA256 5/5 OK; driver_t106.py + spotcheck/indgates.py byte-identical to t106c; gates.py diff vs t106c = guard/F1 pin/F2 header logic/docstrings only (no gate math).
F2a PASS (copy of live 10-col fixes_ledger, 11,362 rows; gates copy with only the live-path pin redirected to scratch):
 - env unset -> 10-col row appended, file stays all-10, no err.  NCOLS=12 on 10-col -> refused + logged, file unchanged.  NCOLS=10 -> appended.
 - unknown header (3 cols) -> refused + logged.  New file: NCOLS=10 -> 10-col header; default -> 12-col. Header read + append are under the same flock; file reopened each write (sees an installed swap).
F2b PASS: normalize_ledger: same-dir tmp (O_EXCL), fsync, os.replace, dir fsync; input header must be C10 or FLCOLS; lossless (11,362/11,362 rows, all original fields preserved, lemma cols ''), idempotent (12->12 identical sha e72039d8…). OUT under rh_extend (incl. ledger.csv) refused; --install refuses any basename != fixes_ledger.csv, runs under DST.lock (same lock as gates), post-sha==SRC sha verified; symlink DST replaces the link, target untouched. LAUNCH.md: no cp onto live (cp only makes backups).
Tests (scratch copy): test_f2 PASS, test_guard_f PASS, test_guard PASS.
Speed: header read ~10 µs vs ~197 µs per flog write (fsync-dominated); writes fire only on gate events -> negligible; per-process inode+size cache NOT needed (and re-read is safer across installs).
Nits (non-blocking): normalize guards use assert (disabled under python -O) -> prefer explicit raise; install pre-sha is printed but not compared to an expected value; lemma cols '' vs writer 'n/a' (Nathan's spec).
Certified ledger.csv untouched (mtime 2026-11-02 22:38 MST). Header choice (10 vs 12) still Nathan's.
