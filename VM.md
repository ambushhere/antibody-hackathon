# VM команды team-3

## Доступ

```bash
ssh -i ~/.ssh/team-3 hackathon@195.242.29.218
```

Ключ лежит в `~/.ssh/team-3` (права 600), в репозиторий не коммитится (`.gitignore`).
Пользователь и домашняя папка общие на всю команду, один GPU на всех. Своя папка: `~/almaz`.

## Железо

- NVIDIA H100 80GB HBM3, драйвер 580.173.02, CUDA 13.0
- 16 vCPU, 196 ГБ RAM, 364 ГБ свободно на диске
- Ubuntu, ядро 6.11 nvidia, Python 3.12.3 системный

## Что уже стоит

- `foundry-toolkit` (`/usr/local/bin`), nextflow 26.10.5, docker с рабочим `--gpus all`, uv 0.13.7, pixi 0.84.2
- llama.cpp в `/opt/llama.cpp`
- `~/foundry-assistant`: локальный помощник по докам Foundry на двух Qwen2.5-0.5B GGUF (оба файла на месте). Запуск по README, порт 8771, `NGL=0` в окружении значит CPU по умолчанию
- Креды quay/foundry прописаны в `~/.bashrc` как `FOUNDRY_*` и `QUAY_*`

Важно: `~/.bashrc` в неинтерактивном shell обрывается до блока с переменными, поэтому через
`ssh host 'команда'` они не подхватываются. Обходной путь:

```bash
eval "$(grep -E '^export (FOUNDRY|QUAY)' ~/.bashrc)"
```

`foundry-toolkit module list` без флагов падает с `no stored credential`, нужен явный
`--registry $FOUNDRY_REGISTRY --username $FOUNDRY_USERNAME --password $FOUNDRY_PASSWORD`.

## Доступные модули Foundry

Все образы уже скачаны в локальный docker, качать ничего не нужно.

| модуль | версии | что это |
|--------|--------|---------|
| apheris-rfd4-proteina | v0.3.0, `-antibody-2026-10-30`, `-her2-2026-11-07`, `-her2-novel-2026-11-09` | генерация дизайнов RFD4-Proteina, есть чекпойнты под антитела и под HER2 |
| apheris-boltz | v1.3.0, `-boltz1`, `-boltz2` | кофолдинг, Boltz-1 и Boltz-2 с confidence head |
| apheris-openfold3 | v1.3.0 и v1.1.4, суффиксы `-of3p2`, `-of3ob` | OpenFold3 preview2 и OpenBind |
| apheris-protenix | v1.3.0, `-v1-default`, `-v2` | Protenix |
| apheris-esmfold2 | v1.4.0, `-release` | ESMFold2 с бэкбоном ESMC-6B в bfloat16 |
| apheris-mol1 | v1.3.0, `-base` | Mol-1 |
| apheris-msa | v1.1.0, v1.0.0 | единый доступ к внешним MSA-провайдерам |
| apheris-data | v1.2.0, v1.1.3 | библиотека для биомолекулярных данных |

`workflow list` и `skill list` пустые, то есть nextflow-пайплайнов и agent skills нам не выдали,
модули запускаются напрямую через docker.

Наличие чекпойнтов `her2` и `her2-novel` от 7 и 9 ноября намекает, что задача дня
может быть связана с HER2. Не факт, но стоит держать в голове.

### Схема вызова модуля

У каждого модуля контракт `foundry.cli.v1`: задача, `--input` (файл или папка),
`--output` (папка), коды выхода (0 успех, 2 плохой вызов, 3 плохие входные данные,
4 нет нужного ресурса, 5 результат не сохранился). Пример для RFD4:

```bash
apheris-rfd4-proteina design --input ./in --output ./run
```

Схему конкретного модуля целиком:

```bash
foundry-toolkit module describe --registry $FOUNDRY_REGISTRY \
  --username $FOUNDRY_USERNAME --password $FOUNDRY_PASSWORD \
  --name apheris-rfd4-proteina --version v0.3.0-her2-novel-2026-11-09
```

## Наша папка на VM

```
~/almaz/work/.venv     uv venv на Python 3.12, torch 2.15.0+cu132 видит H100
~/almaz/work/repo      копия src/ из этого репозитория
~/almaz/data           входные данные
~/almaz/out            результаты
```

Запуск:

```bash
cd ~/almaz/work/repo && . ../.venv/bin/activate && python src/mini_task.py
```
