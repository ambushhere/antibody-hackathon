#!/usr/bin/env bash
set -euo pipefail
cd /home/hackathon/almaz/work/neural-20261114
source /home/hackathon/almaz/work/.venv/bin/activate
export PYTHONPATH="$PWD"
export PYTHONUNBUFFERED=1
exec 9>/home/hackathon/almaz/work/gpu.lock
flock -n 9
if [ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
  echo 'GPU occupied; run not admitted.' >&2
  exit 3
fi
python scripts/train_her2_neural.py --data-dir /home/hackathon/almaz/data/hackathon-data --output runs/large-seed3 --device cuda --size large --epochs 15 --seed 3 --models transformer
