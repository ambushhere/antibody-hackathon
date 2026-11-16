# CDRH3 intuition

Two small illustrations accompany the README: a static view of the design window and a quiet animation of four substitutions. Both use locally bundled Geist from the [official Vercel font repository](https://github.com/vercel/geist-font), with its SIL Open Font License. Original company marks remain unchanged.

- [Where the window sits](../../docs/assets/intuition/loop.png): whole IgG, one VH/VL binding site and the window within CDRH3. The six CDR loops are shown schematically, not as an atomic structure.
- [Substitutions](../../docs/assets/intuition/mutations.gif): four letters change while the scaffold stays fixed. [MP4](../../docs/assets/intuition/mutations.mp4) · [Static version](../../docs/assets/intuition/mutations.png).

The illustrative sequence starts at `WGGDGFYAMD`. G2A, D4N, Y7F and M9L produce `WAGNGFFALD`. These positions refer to the ten-residue library window, not IMGT numbering. The example is not a selected candidate and makes no claim of novelty eligibility or better binding.

## Reproduce

Requires Node.js 22+, FFmpeg and a Chromium browser supported by HyperFrames. The CLI is pinned in each project. Source, fonts and GSAP are local; initial CLI/browser installation may require a download. Font sources and hashes are recorded in [font-sources.json](font-sources.json).

```bash
bash scripts/presentation/render-readme-animation.sh
```

The script checks both illustrations, captures the static loop diagram, renders the mutation clip at 960 × 440 and 24 fps, creates a 12 fps GIF and extracts its static version. It overwrites the four derived files in `docs/assets/intuition/`. No model jobs, provider credentials, final slide edits or sequence edits are involved.

Inspect the substitution transitions and the settled frame before accepting edits. Short closing-line warnings at the loop reset concern already explained text, not a new caption. The drafts include briefs and timing plans; only the mutation illustration is animated in the README.
