# AI x Antibodies Hackathon (Nucleate Norway, 14.11.2026)

Контекст и чек-лист: [CLAUDE.md](CLAUDE.md)

## Окружение
Intel Mac → PyTorch максимум 2.2.2, Python 3.12, numpy<2.

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
source .venv/bin/activate
```

## Важно: SSL для скачивания весов
У Homebrew-Python нет корневых сертификатов, `torch.hub` падает с `CERTIFICATE_VERIFY_FAILED`.
Перед запуском чего-либо, что качает веса:

```bash
export SSL_CERT_FILE=$(python -c "import certifi;print(certifi.where())")
export REQUESTS_CA_BUNDLE=$SSL_CERT_FILE
```

## Нумерация CDR
`abnumber` требует классический ANARCI с HMMER — не ставится без возни.
Вместо него `anarcii` (чистый Python/PyTorch). Обёртка: [src/cdr_tools.py](src/cdr_tools.py).

```bash
.venv/bin/python src/cdr_tools.py
```

## Структура
- `data/` — входные PDB/последовательности
- `notebooks/` — исследования
- `src/` — пайплайн
- `results/` — кандидаты, метрики, графики
