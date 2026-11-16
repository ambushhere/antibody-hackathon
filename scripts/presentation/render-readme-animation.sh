#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
asset_dir="$repo_root/docs/assets/intuition"
mkdir -p "$asset_dir"
for fragment in loop mutations; do
  npx --yes hyperframes@0.9.56 check "$repo_root/videos/her2-intuition/$fragment"
done
snapshot_dir="$(mktemp -d)"
trap 'rm -rf "$snapshot_dir"' EXIT
npx --yes hyperframes@0.9.56 snapshot "$repo_root/videos/her2-intuition/loop" --at 5.7 --output "$snapshot_dir"
cp "$snapshot_dir/frame-00-at-5.7s.png" "$asset_dir/loop.png"
npx --yes hyperframes@0.9.56 render "$repo_root/videos/her2-intuition/mutations" --fps 24 --quality delivery --output "$asset_dir/mutations.mp4"
ffmpeg -hide_banner -loglevel error -y -i "$asset_dir/mutations.mp4" \
  -filter_complex '[0:v]fps=12,split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];[b][p]paletteuse=dither=sierra2_4a:diff_mode=rectangle' \
  -loop 0 "$asset_dir/mutations.gif"
ffmpeg -hide_banner -loglevel error -y -ss 5.7 -i "$asset_dir/mutations.mp4" -frames:v 1 -update 1 "$asset_dir/mutations.png"
ffprobe -v error -show_entries format=duration,size -of json "$asset_dir/mutations.mp4"
