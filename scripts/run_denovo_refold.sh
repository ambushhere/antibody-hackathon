#!/usr/bin/env bash
# Matched target-only template diagnostic; no antibody pose is supplied.
set -euo pipefail
campaign=${1:?absolute campaign directory on team VM}
out="$campaign/denovo-target-only-protenix-v3"
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}' > "$campaign/denovo-refold-image-id.txt"
sha256sum "$campaign/denovo-refold-queries.json" "$campaign/her2_target_refold.cif" > "$campaign/denovo-refold-input-v3.sha256"
docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
  "$image" module-run predict --input "$campaign/denovo-refold-queries.json" --output "$out" \
  --seeds 42 --num-diffusion-samples 1 --write-detailed-scores true \
  --template-paths "$campaign/her2_target_refold.cif" \
  > "$campaign/denovo-refold-v3.log" 2>&1
