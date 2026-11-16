#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
project_dir="$repo_root/videos/her2-pipeline"
asset_dir="$repo_root/docs/assets"
mkdir -p "$asset_dir"
(
  cd "$project_dir"
  npx --yes hyperframes@0.9.56 check .
  npx --yes hyperframes@0.9.56 render . --fps 24 --quality delivery --output "$asset_dir/her2-pipeline.mp4"
)
ffmpeg -hide_banner -loglevel error -y -i "$asset_dir/her2-pipeline.mp4" \
  -filter_complex '[0:v]fps=12,scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];[b][p]paletteuse=dither=sierra2_4a:diff_mode=rectangle' \
  -loop 0 "$asset_dir/her2-pipeline.gif"
ffmpeg -hide_banner -loglevel error -y -ss 6.5 -i "$asset_dir/her2-pipeline.mp4" \
  -frames:v 1 -update 1 "$asset_dir/her2-pipeline-poster.png"
ffprobe -v error -show_entries format=duration,size -of json "$asset_dir/her2-pipeline.mp4"
