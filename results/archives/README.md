# Structural run archive

## Control and dataset checkpoint

The later control experiment and sequence audit have separate archives so they can be restored independently.

| Manifest | Preserved content | Restore root |
| --- | --- | --- |
| [template-controls-manifest.json](template-controls-manifest.json) | All 48 no-MSA Protenix poses, confidence arrays, logs and receipts, plus the pinned original 1N8Z and 6OGE deposits | `results/` |
| [template-msa-controls-manifest.json](template-msa-controls-manifest.json) | All 24 target-MSA control poses and their runtime receipts | `results/` |
| [denovo-scale-manifest.json](denovo-scale-manifest.json) | All 64 generated structures from the two pre-existing team VM RFD4 runs and their provenance | `results/` |
| [neural-seed3-weights-manifest.json](neural-seed3-weights-manifest.json) | Selected fourth-seed large Transformer weights, metadata under `results/neural-runs/large-seed3/` | `results/` |
| [abguard-installed-package-manifest.json](abguard-installed-package-manifest.json) | Exact installed abguard Python package and bundled reference, with dependency versions recorded separately in `results/abguard-validation/raw/environment.json` | `results/tmp/abguard-installed/` |
| [external-spr-sources-manifest.json](external-spr-sources-manifest.json) | Five pinned Absci source files, including source README and Clear BSD license | `data/` |

For these manifests, select a filename below. Archive and member hashes are checked before use. The original reference deposits restore to `results/campaign/template-references/{1n8z,6oge}.cif`, allowing the template and candidate auditors to run without the original machine's `/tmp` files.

```bash
source .venv/bin/activate
python - <<'PY'
import hashlib, json, tarfile
from pathlib import Path
folder = Path('results/archives')
manifest = json.loads((folder / 'template-controls-manifest.json').read_text())
archive = folder / manifest['file']
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['sha256']
destination = Path(manifest['extract_under'])
destination.mkdir(parents=True, exist_ok=True)
with tarfile.open(archive) as bundle:
    for item in manifest['members']:
        data = bundle.extractfile(item['path']).read()
        assert len(data) == item['bytes']
        assert hashlib.sha256(data).hexdigest() == item['sha256']
    bundle.extractall(destination, filter='data')
print('Verified and restored', len(manifest['members']), 'files')
PY
```

## Earlier structural runs

The archive preserves 307 files exported from the completed VM pilots and JapanFold RF3 run, including structures, scores, detailed confidence arrays, input files and execution logs. The source organizer dataset, model weights and credentials are excluded. A duplicate JapanFold download zip is omitted because its extracted contents are included.

`structure-runs-manifest.json` records the archive SHA-256 and every member's relative path, byte count and SHA-256. The archive is about 83 MB because it retains compressed detailed confidence arrays needed for later audits. Extract it from the repository root into the ignored raw results directory:

```bash
source .venv/bin/activate
python - <<'PY'
import hashlib, json, tarfile
from pathlib import Path
folder = Path('results/archives')
manifest = json.loads((folder / 'structure-runs-manifest.json').read_text())
archive = folder / manifest['archive']
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['sha256']
destination = Path('results/campaign')
destination.mkdir(parents=True, exist_ok=True)
with tarfile.open(archive) as bundle:
    bundle.extractall(destination, filter='data')
for member in manifest['files']:
    path = destination / member['path']
    assert path.stat().st_size == member['bytes']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == member['sha256']
print('Verified', len(manifest['files']), 'files')
PY
```

Use Python 3.12 or newer. The archive includes all 21 Protenix context-panel predictions and all seven RF3 predictions. Completion and preservation do not establish biological validity; the matched-context structural interpretation remains a follow-up.
