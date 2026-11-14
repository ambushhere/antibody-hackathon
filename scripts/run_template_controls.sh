#!/usr/bin/env bash
# Root-owned sequential GPU admission. Existing outputs are never overwritten.
set -euo pipefail
input=${1:?absolute input directory}
out=${2:?absolute new output directory}
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
exec 9>/home/hackathon/almaz/work/gpu.lock
flock -n 9 || { echo 'GPU admission lock is held' >&2; exit 2; }
[[ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]] || { echo 'GPU has an active process' >&2; exit 2; }
mkdir -p "$out"
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}' > "$out/image-id.txt"
sha256sum "$input"/{complete-target.cif,manifest.json,masked-target.cif,queries.json} > "$out/input.sha256"
date -u +%FT%TZ > "$out/started.txt"
for seed in 42 43 44; do
  for arm in complete masked; do
    run="$out/$arm/seed-$seed"
    test ! -e "$run"
    mkdir -p "$out/$arm"
    docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
      -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
      "$image" module-run predict --input "$input/queries.json" --output "$run" \
      --seeds "$seed" --num-diffusion-samples 2 --write-detailed-scores true \
      --template-paths "$input/$arm-target.cif" > "$run.log" 2>&1
  done
done
date -u +%FT%TZ > "$out/completed.txt"
