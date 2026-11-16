# Data retention and restoration

Post-event handoff, 16 November 2026. The final presentation, ranked candidates, code, compact result tables, provenance and compressed experiment archives are retained. The cleanup removed local reproducible downloads, scratch environments and verified expanded duplicates. It did not rewrite Git history or access event infrastructure.

## What remains available

| Material | Location | Recovery boundary |
| --- | --- | --- |
| Final manual two-slide pitch and exact previews | [Presentation](../output/selection-slide/README.md) | Tracked in the handoff revision; do not regenerate it with the older four-slide builder |
| Final 30-row ranked CSV | [epitope-hunters.csv](../output/submission-draft/epitope-hunters.csv) | Identical to the frozen ranking export and installed-validator receipt |
| Selection inputs, scores and validation | [Experiment index](EXPERIMENT_INDEX.md), `results/`, `output/submission-draft/` | Enough for deterministic selection replay without original datasets or API access |
| Raw structures, trained neural weights and installed abguard package/reference | [Archive guide](../results/archives/README.md) | All 15 archive hashes and listed member hashes were checked before cleanup |
| Original five-file benchmark inputs | Local ignored `data/raw/retained-inputs/` | One 8,418,305-byte compressed archive plus SHA-256 manifest. Local only, not distributed by this repository |
| Original development environment | Main checkout `.venv`, shared by this worktree | Retained; rebuilding it is separate from restoring datasets |

Unique scratch outputs without verified archival equivalents were retained. Organizer credentials, other projects' caches and shared runtime installations were left untouched. Historical guides sometimes refer to deleted scratch paths; restore those inputs or pass explicit replacement paths before rerunning a full experiment.

## What was removed locally

- The downloaded ESM-2 650M checkpoint and contact-regression cache. Both upstream download endpoints returned HTTP 200 during the audit. Calling the loader in `src/metrics.py` downloads them again; it requires network access and disk space.
- Downloaded AbLang2 assets and its isolated temporary environment. Restore using the pinned source/weight URLs and hashes in the [AbLang protocol](research/ablang-feature-protocol.md) and `requirements-ablang.txt`.
- Downloaded ThermoMPNN-D weights/archive, its isolated temporary environment and verified duplicate inputs/results. Restore from the pinned sources and checksum manifests in the [thermal pilot](../results/thermal-pilot/README.md). The modified source tree and unique run records were retained.
- Both local copies of the two RFD4 metadata tables. Exact upstream revision, URLs and SHA-256 hashes remain in the [audit](../results/rfd4-training-metadata/audit.json); follow the [replay instructions](research/rfd4-training-data.md#replay).
- The uncompressed five-file training inputs and their temporary TAR/ZIP copies, after validating the compact local archive above.
- 725 expanded untracked files whose bytes matched the retained experiment archives, plus project Python/pytest caches. The retained archive manifests identify each restorable member.

Deleted file contents totalled **4.51 GB (4.20 GiB)** before allowing for the new 8.4 MB input archive. The filesystem reported approximately **4.08 GiB more available space** across the cleanup interval. These measurements differ because of APFS allocation, shared blocks and concurrent system activity; logical file size is not a guarantee of reclaimed physical space.

## Reproduce only what you need

The [README selection replay](../README.md#reproduce-the-final-selection) needs only committed evidence and RapidFuzz. It was run in a fresh environment and reproduced the final CSV byte for byte. It does not rerun abguard, generation or training.

For raw structural evidence, use the hash-verifying extraction commands in the [archive guide](../results/archives/README.md). Restore one relevant archive at a time. The archived raw evidence is enough to inspect prior predictions; reproducing provider inference may require models, hardware and access that are no longer available.

For the local retained training inputs, use Python 3.12 or newer from the repository root. This restores into a new temporary directory and verifies every input before use:

```bash
source .venv/bin/activate
python - <<'PYCODE'
import hashlib, json, tarfile, tempfile
from pathlib import Path
folder = Path('data/raw/retained-inputs')
manifest = json.loads((folder / 'manifest.json').read_text())
archive = folder / manifest['archive']
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['sha256']
destination = Path(tempfile.mkdtemp(prefix='her2-training-'))
with tarfile.open(archive) as bundle:
    for member in manifest['members']:
        data = bundle.extractfile(member['path']).read()
        assert len(data) == member['bytes']
        assert hashlib.sha256(data).hexdigest() == member['sha256']
    bundle.extractall(destination, filter='data')
print(destination)
PYCODE
```

A fresh clone does not contain that local-only archive. Obtain the exact organizer preparation from an authorized retained copy, then verify the input hashes recorded in [SPR metrics](../results/spr-benchmark/metrics.json) and the [original benchmark report](../results/experiment-20261114/benchmark_report.md). If those inputs cannot be recovered, the preserved selection can still be replayed, but full training reproduction is unavailable. The published dataset alone does not guarantee byte-identical organizer preprocessing.

Before cleanup, the local suite passed 311 tests with four optional skips. After removing expanded archives it passed 309 with six skips. The two additional skips are archive-dependent RF3/OpenDDE checks in `tests/test_structural_replicates.py`; their raw inputs can be restored from the retained archives. The same 22 subtests passed in both runs.
