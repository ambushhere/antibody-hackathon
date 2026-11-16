# Animated pipeline overview

An eight-second silent illustration for the repository README, authored at 1200 × 640 and 24 fps in HyperFrames 0.9.56. The README uses a smaller 960-pixel GIF at 12 fps. A static poster is available for readers who prefer no motion.

- [Animated GIF](../../docs/assets/her2-pipeline.gif)
- [MP4 video](../../docs/assets/her2-pipeline.mp4)
- [Static poster](../../docs/assets/her2-pipeline-poster.png)

The teal route depicts optimization of the ten-residue trastuzumab library window, informed by FACS screening classes and SPR binding rank. Its sequence letters remain the wild-type reference; highlighting marks design variables, not a proposed new sequence. The amber route depicts RFD4-Proteina generation of paired variable domains, followed by sequence checks and diversity. Both routes lead to the recorded 30-candidate list. Laboratory performance is not established by the animation. Domain outlines are schematic, not predicted structures.

## Render locally

Node.js 22 or newer, FFmpeg and a Chromium browser supported by HyperFrames are required. The composition bundles GSAP and local font files, so it does not fetch brand artwork, fonts or model data while rendering. Their upstream notices remain in `assets/`.

From the repository root:

```bash
bash scripts/presentation/render-readme-animation.sh
```

The script checks the composition, renders the MP4, creates the smaller GIF and extracts the poster at 6.5 seconds. It overwrites those three derived artifacts. No API credentials, model calls or GPU training are involved. HyperFrames CLI/browser installation may need network access on a new machine.

For review without changing the exports:

```bash
(cd videos/her2-pipeline && npm run check)
(cd videos/her2-pipeline && npx --yes hyperframes@0.9.56 snapshot --at 0,1.5,3,5,6.5,7.8)
```

The final two-slide PowerPoint is separate and unchanged. Company logos are sourced and attributed in the [event artwork directory](../../docs/assets/brands/README.md).
