#!/usr/bin/env bash
set -euo pipefail
cd /home/hackathon/almaz/work/spr-adaptation-20261114
source /home/hackathon/almaz/work/.venv/bin/activate
export PYTHONPATH="$PWD"
export PYTHONUNBUFFERED=1
exec 9>/home/hackathon/almaz/work/gpu.lock
flock -n 9
if [ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
  echo 'GPU occupied; run not admitted.' >&2
  exit 3
fi
python -m scripts.benchmark_spr_finetune --data-dir data --checkpoint /home/hackathon/almaz/work/neural-20261114/runs/large-seed0/transformer/best.pt --output results/spr-finetune-seed0 --device cuda
