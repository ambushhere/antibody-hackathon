# Pipeline explanations

The two four-page guides explain the generation and evaluation slides, the measured and synthetic datasets, current limitations and terminology. Both include diagrams, speaking notes and links to the underlying evidence. They describe a research snapshot on 14 November 2026, not laboratory validation.

- [Russian PDF](her2-pipeline-explained-ru.pdf), [Markdown source](../../docs/PIPELINE_EXPLAINED_RU.md)
- [English PDF](her2-pipeline-explained-en.pdf), [Markdown source](../../docs/PIPELINE_EXPLAINED_EN.md)
- [Two-slide presentation](../pipeline-draft/README.md)

## Rebuild

Activate the project environment and use Python with `reportlab` installed:

```bash
source .venv/bin/activate
python scripts/presentation/export_pipeline_guides.py
```

The exporter uses the macOS Arial fonts, including Cyrillic, from `/System/Library/Fonts/Supplemental`. On another host, pass `--font-dir /path/to/fonts` containing `Arial.ttf` and `Arial Bold.ttf`. Markdown content is read from the two source documents. The two diagram layouts and their labels are maintained explicitly in the exporter; update them alongside the Markdown diagrams when results change.

The PDFs embed fonts, contain searchable text and use vector diagrams. All eight rendered pages were visually reviewed. The source links point to the repository's `main` branch and require repository access.
