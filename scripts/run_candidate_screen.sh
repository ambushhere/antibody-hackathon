#!/usr/bin/env bash
# Root-owned first screen, only after the fixed control reproduction gate passes.
set -euo pipefail
input=${1:?absolute input directory}
out=${2:?absolute new output directory}
admission=${3:?path to the recorded control admission JSON}
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
source /home/hackathon/almaz/work/.venv/bin/activate
python - "$admission" "$input" <<'PY'
import hashlib, json, sys
from pathlib import Path
with open(sys.argv[1]) as stream:
    receipt = json.load(stream)
control_path = Path(receipt['control_audit_path'])
control = json.loads(control_path.read_bytes())
msa = receipt['selected_arm'] == 'complete_target_msa'
checks = [
    receipt['admitted'] is True,
    receipt['selected_arm'] in ('complete_no_msa', 'complete_target_msa'),
    receipt['expected_poses'] == 57,
    hashlib.sha256(control_path.read_bytes()).hexdigest() == receipt['control_audit_sha256'],
    hashlib.sha256((Path(sys.argv[2]) / 'manifest.json').read_bytes()).hexdigest() == receipt['input_manifest_sha256'],
    control['status'] == 'complete',
    control['msa_policy'] == ('target_only' if msa else 'none'),
    control['operational_rule']['arms']['complete']['qualifies'] is True,
    len(control['operational_rule']['arms']['complete']['passing_wt_top1_seeds']) >= 2,
]
if not all(checks):
    raise ValueError('Candidate screen does not have a matching qualified control audit')
PY
exec 9>/home/hackathon/almaz/work/gpu.lock
flock -n 9 || { echo 'GPU admission lock is held' >&2; exit 2; }
active=$(nvidia-smi --query-compute-apps=pid,used_gpu_memory --format=csv,noheader,nounits)
if [[ -n "$active" ]]; then
  [[ "${HER2_ALLOW_SHARED_GPU:-0}" == 1 ]] || { echo 'GPU has an active process' >&2; exit 2; }
  free_mib=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits)
  [[ "$free_mib" =~ ^[0-9]+$ && "$free_mib" -ge 40000 ]] || { echo 'Insufficient shared GPU headroom' >&2; exit 2; }
fi
mkdir -p "$out"
printf '%s\n' "$active" > "$out/preexisting-gpu-processes.csv"
printf '%s\n' "${HER2_ALLOW_SHARED_GPU:-0}" > "$out/shared-gpu-admission.txt"
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}' > "$out/image-id.txt"
sha256sum "$input"/{complete-target.cif,manifest.json,queries.json} > "$out/input.sha256"
if [[ -f "$input/target.a3m" ]]; then sha256sum "$input/target.a3m" >> "$out/input.sha256"; fi
cp "$admission" "$out/admission.json"
date -u +%FT%TZ > "$out/started.txt"
docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
  "$image" module-run predict --input "$input/queries.json" --output "$out" \
  --seeds 42 --num-diffusion-samples 1 --write-detailed-scores true \
  --template-paths "$input/complete-target.cif" > "$out/run.log" 2>&1
date -u +%FT%TZ > "$out/completed.txt"
