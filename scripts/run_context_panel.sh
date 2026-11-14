#!/usr/bin/env bash
# Matched initial screen. Root owns the GPU and launches this script once.
set -euo pipefail
input=${1:?absolute directory with configs/her2-context files}
out=${2:?absolute new output directory}
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
mkdir -p "$out"
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}' > "$out/image-id.txt"
sha256sum "$input"/* > "$out/input.sha256"
for arm in split domain-iv ecd; do
  docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
    "$image" module-run predict --input "$input/$arm-queries.json" --output "$out/$arm" \
    --seeds 42 --num-diffusion-samples 1 --write-detailed-scores true \
    --template-paths "$input/$arm-target.cif" > "$out/$arm.log" 2>&1
done
