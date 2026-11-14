# Structural run archive

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
