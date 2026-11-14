#!/usr/bin/env bash
# Run on team VM. Never overwrite an existing run.
set -euo pipefail
track=${1:?optimization or denovo}
input=${2:?absolute condition spec path}
out=${3:?absolute new output directory}
n=${4:-4}
seed=${5:-20261114}
case "$track" in
  optimization) tag=her2-novel-2026-11-09 ;;
  denovo) tag=antibody-2026-10-30 ;;
  *) exit 2 ;;
esac
[[ ! -e "$out" ]] || { echo 'Output already exists' >&2; exit 2; }
mkdir -p "$out"
image="quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-$tag"
printf '%s\n' "$image" > "$out/image.txt"
docker image inspect "$image" --format '{{.Id}}' > "$out/image-id.txt"
sha256sum "$input" > "$out/input.sha256"
printf 'track=%s\nn=%s\nseed=%s\n' "$track" "$n" "$seed" > "$out/parameters.txt"
docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
  "$image" module-run design --input "$input" --output "$out/design" \
  --num-replicates "$n" --seed "$seed" --recycle-strategy high/front_loaded \
  > "$out/design.log" 2>&1
