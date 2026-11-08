# Referee gates: reproducible, isolated verification
Nathan A. Meeks & Grok. [NUM]/[LEM] candidate, unrefereed; no RH claim.

## Trust and acceptance

The standalone evaluator is not an entire database/referee pipeline. Internal Audit is an engineering/numerical check, not external peer review. Review the mathematical assumptions in `proof-runtime/WHITEPAPER.md`, `GATES.md`, `C0_FIX.md`, `NSTEP_FIX.md`, and included audit reports before treating a candidate numerical record as a theorem.

For ordinary windows, whole-window Lemma CP needs certified signs and a matching Turing count; repaired blocks additionally require independently pinned G-Turing counts, disjoint opposite-sign G-RS brackets, no conflicts and no unresolved cuts (or the documented G-Platt containment route). Failure/ambiguity means refusal; never skip a gap. Preserve the 1e-8 GT_BATCH arithmetic slack (never below 1e-9).

## Portable conservative reference path

Install python-flint and extract the proof database (DATABASE.md). The supplied `spotcheck/indgates.py` works without custom fast binaries by selecting its existing Arb fallback. Disable all optional pinned/path-specific accelerators; do **not** edit hash allowlists or libraries just to make a build pass. Run from a new scratch directory; only the supplied bracket file is read. For an entire committed example window:

```sh
env PLATT=0 GTURING_P=0 GRS_T100=0 GRS_TIER1=0 GT_BATCH=0 GRS_ESC=0   HYB_VC=0 HYB_TMESH=1 PC_FAST=0   python3 /ABS/REPO/proof-runtime/spotcheck/indgates.py   2999999999920 3000000000000 3000000030960   12362403019661 12362403152169   /ABS/DATA/proof-data/brackets/w0000000_T2999999999920.txt.xz   external-referee 1 > referee-window0.json
```

This checks every bracket in the window, not a sampled subset. Expect `result=RESOLVED`, need=132508, matching endpoint counts, 0 conflicts/unresolved cuts; **the whole command was not rerun during preparation** (can take hours on the Arb fallback). Inspect JSON explicitly; the CLI's exit status is not the gate verdict. `GRS_ESC=0` also prevents the historical optional residue/Platt escalation: an endpoint still unresolved after the existing G-cut checks is refused, never accepted with fewer checks. This may refuse more windows; it does not weaken acceptance. Independent Arb signs on the first bracket were checked as a bounded smoke test. The complete mathematical database was not freshly recertified by this preparation.

To verify any other committed window, take `T_arg`, `lo`, `hi`, `Na`, `Nb` and filename from `committed_ledger.csv`; pass exactly those values. Brackets store offsets u, with exact t=T_arg+u; preserve decimal strings, never round large t to floating-point before sign evaluation. G-cut is invoked for uncertain endpoints. A refusal requires investigation, not converting the result to PROVEN.

## Historical accelerated paths (sources provided, not portability promises)

`grs-src/` contains `rs_fast.c`, `rs_tb_fix.c`, `rsc.h`, `rscv.h` and the historical compiler flags in `grs-build-historical.sh`. `spotcheck/rsbatch_t100.c` supplies the older batch tier. The fast tiers require suitable AVX CPU features and exact pinned binary/library hashes. Rebuilds on a different compiler/library can be rejected and fall back to Arb; this is intentional. No need to enable them for the conservative verification above. G-Turing-P sources are in `v9p/gturing_p/`; the Platt wrapper source is `platt_wrap.c` (requires the corresponding newer FLINT API, historical Python FLINT 3.6). Hardcoded old paths are historical and are not rewritten in this snapshot. Do not use production drivers to replay data, change live locks, or restart a held run.

Keep **HYB_VC=0, PC_FAST=0** for referee handoff. Historic hybrid speed reports sometimes used shared caching; those timings must not be represented as independent G-Turing/G-RS verification. Original production files and hashes are preserved unchanged; a future portability patch requires separate code ownership and independent audit.
