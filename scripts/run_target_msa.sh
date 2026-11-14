#!/usr/bin/env bash
# Root-owned target-MSA comparison, after the no-MSA control batch.
set -euo pipefail
input=${1:?absolute input directory}
out=${2:?absolute new output directory}
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
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
sha256sum "$input"/{complete-target.cif,manifest.json,target.a3m,queries.json} > "$out/input.sha256"
date -u +%FT%TZ > "$out/started.txt"
run="$out/complete"
docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
  "$image" module-run predict --input "$input/queries.json" --output "$run" \
  --seeds 42 --seeds 43 --seeds 44 --num-diffusion-samples 2 --write-detailed-scores true \
  --template-paths "$input/complete-target.cif" > "$run.log" 2>&1
date -u +%FT%TZ > "$out/completed.txt"
