#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
command -v uv >/dev/null || { echo 'Install uv first: https://docs.astral.sh/uv/getting-started/installation/'; exit 1; }
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements-benchmark.txt
[[ -f .env ]] || cp .env.example .env
chmod 600 .env
printf '%s\n' 'Environment ready. Fill .env from team credentials. GitHub Actions secrets cannot be downloaded.'
