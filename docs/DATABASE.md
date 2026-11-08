# Complete committed database: assets and verification

Candidate/unrefereed numerical proof record, Nathan A. Meeks & Grok. Not RH. The committed extent is in `data/BOUNDARY.json`; window165468 and failed attempt rows are excluded.

## Asset layout / free-storage plan

Source/docs/chart and the compressed-file manifest belong in ordinary git. The ~20 GB bulk database is prepared as **one-GiB `proof-data.tar.partNNN` release assets** plus `SHA256SUMS`. No Git LFS and no paid storage. GitHub's release documentation (checked 2026-11-07 MST) allows 1000 assets, each <2 GiB, without a total release-size/download-bandwidth quota: https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases . Actual asset publication is separate from preparation; no download URL or tag is asserted here before it exists.

After a release is actually approved and published, obtain every `proof-data.tar.part*` and `SHA256SUMS` from that same release (not GitHub's automatic source-only ZIP). Verify before extracting:

```sh
sha256sum -c SHA256SUMS
cat proof-data.tar.part* | tar -xf -
```

Allow >21 GB free for the extracted compressed database in addition to downloaded assets; much more space is needed if decompressing all bracket text. Fixed-width three-digit part numbering preserves concatenation order. All tar entries are regular files under `proof-data/`, no symlinks or absolute paths.

Contents: all 165,468 committed `.txt.xz` bracket files; normalized `committed_ledger.csv`; exact `checkpoint.json`; 12-column `fixes_ledger.csv` filtered to committed indices; committed-only `windows.jsonl` and `gates.jsonl`; projected `commit_indices.jsonl` containing commit index/time (not giant scheduler finished lists); `BOUNDARY.json`; `BRACKETS_SHA256.tsv`.

The legacy ledger header has 17 columns while later rows add a raw bracket SHA and unquoted notes. Normalization preserves core 17 fields, records the SHA when present, joins remaining note fragments, and writes a valid 19-column CSV. The live ledger is never changed. Selected rows have status PROVEN and idx<checkpoint.next_idx and are exactly 0..165467, without duplicate committed rows. Coverage and N-count chaining, overlap counts and net total were checked against the checkpoint.

`BRACKETS_SHA256.tsv` gives relative path, compressed byte length, compressed SHA256, historical raw/decompressed SHA256 if present, and ledger bracket count. Verify every compressed file against the manifest. After decompression, compare the raw SHA where available; 98,218 legacy rows have no historical raw SHA, so the newly computed compressed hash is **only an artifact-integrity checksum, not retroactive proof provenance**. Preparation's actual check scope is recorded in the report; never confuse hashing with mathematical sign/count validation.

A complete fresh referee check must independently validate all endpoints and counts plus overlap/coverage chaining; no full rerun is claimed. A sample or audit PASS cannot substitute for external mathematical review.

## Check every compressed artifact against the per-file manifest

Run inside the extracted `proof-data/` directory (Python >=3.11):

```sh
python3 - <<'PYVERIFY'
import csv, hashlib
from pathlib import Path
rows = list(csv.DictReader(open('BRACKETS_SHA256.tsv'), delimiter='\t'))
for r in rows:
    p = Path(r['path'])
    assert not p.is_absolute() and '..' not in p.parts
    assert p.stat().st_size == int(r['bytes']), p
    with p.open('rb') as f:
        assert hashlib.file_digest(f, 'sha256').hexdigest() == r['compressed_sha256'], p
print('PASS', len(rows), 'compressed bracket files; integrity only, not mathematical re-certification')
PYVERIFY
```

Historical failed attempts for indices which were eventually committed can remain in the diagnostic windows/gates journals; they never replace the canonical PROVEN row. Scheduler `finished_so_far` arrays are omitted from the projected commit journal because they are scheduling diagnostics, not proof evidence.
