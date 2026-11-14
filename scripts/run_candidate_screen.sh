#!/usr/bin/env bash
# Root-owned first screen, only after the fixed control reproduction gate passes.
set -euo pipefail
input=${1:?absolute input directory}
out=${2:?absolute new output directory}
admission=${3:?path to the recorded control admission JSON}
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
python3 - "$admission" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    receipt = json.load(stream)
assert receipt['admitted'] is True
assert receipt['selected_arm'] == 'complete_no_msa'
assert receipt['passing_wt_top1_seeds'] >= 2
assert receipt['expected_poses'] == 57
PY
exec 9>/home/hackathon/almaz/work/gpu.lock
flock -n 9 || { echo 'GPU admission lock is held' >&2; exit 2; }
[[ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]] || { echo 'GPU has an active process' >&2; exit 2; }
mkdir -p "$out"
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}' > "$out/image-id.txt"
sha256sum "$input"/{complete-target.cif,manifest.json,queries.json} > "$out/input.sha256"
cp "$admission" "$out/admission.json"
date -u +%FT%TZ > "$out/started.txt"
docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
  "$image" module-run predict --input "$input/queries.json" --output "$out" \
  --seeds 42 --num-diffusion-samples 1 --write-detailed-scores true \
  --template-paths "$input/complete-target.cif" > "$out/predict.log" 2>&1
date -u +%FT%TZ > "$out/completed.txt"
