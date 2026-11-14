# Рунбук: готовые команды

Схемы вытащены из `foundry-toolkit module describe` для всех 14 модулей и из файлов внутри
образов. Форматы входа в схемах описаны только как «file or directory», реальные поля
восстановлены из документации внутри образов.

## Решающая деталь про чекпойнты HER2

Описания из метаданных модулей:

| Сборка | Описание из метаданных | `--model` | `sequence_tsr_k` |
|---|---|---|---|
| v0.3.0 (release) | общая модель, RFD4-Proteina v0.0.3 | `release` | 4 |
| antibody-2026-10-30 | antibody flow checkpoint, из NVIDIA | `antibody` | 4 |
| **her2-2026-11-07** | **antibody checkpoint, дообучен на trastuzumab CDRH3 FACS binders** | `her2-conservative` | 4 |
| | то же, другая резкость сэмплинга | `her2-balanced` | 2 |
| | то же, ещё свободнее | `her2-exploratory` | 1 |
| **her2-novel-2026-11-09** | **antibody checkpoint, preference-trained на парах HER2 high против low binder** | `her2-novel` | 4 |

Из этого следует почти наверняка: задача про **оптимизацию CDR-H3 трастузумаба под HER2**.
Чекпойнт `her2-2026-11-07` дообучен ровно на этих данных, а `her2-novel` обучен на
предпочтениях между сильными и слабыми связывателями.

Это независимо подтверждает стратегию из `STRATEGY.md`: оптимизация существующего
связывателя, а не дизайн с нуля. Организаторы положили нам чекпойнт именно под это.

Важно: пресеты `her2-conservative`, `her2-balanced` и `her2-exploratory` используют **один
и тот же файл весов** и отличаются только параметром `sequence_tsr_k` (4, 2, 1). Это
резкость сэмплинга последовательностей: меньше значит разнообразнее и рискованнее.
Разумный план — прогнать все три и сравнить, это дешёвый источник разнообразия.

В каждой сборке лежит ровно один flow-чекпойнт, поэтому `--model` должен соответствовать
образу, иначе ошибка при запуске.

## Порядок пайплайна

```
design  →  refold-prepare  →  cofolding predict  →  refold-score
```

`refold-prepare` берёт выход `design`, отбирает лучшие (`--top-n` или `--top-fraction`,
по умолчанию 0.2) и пишет готовые query-файлы для кофолдинга. `refold-score` принимает
предсказания и выдаёт CSV с метриками и pass/fail. Для антител есть
`--criterion antibody`.

## Общие правила Foundry

- Entrypoint у образов пустой, CLI внутри вызывается как `module-run <задача>`.
- `--output` обязателен везде, тип directory. **Непустая папка, созданную не инструментом,
  вызывает ошибку**, поэтому на каждый запуск своя свежая папка.
- `--dry-run true` проверяет только флаги и конфиг. Ошибки в query, MSA и чекпойнте он не
  ловит, модель не загружается.
- Коды выхода одинаковы: 0 успех, 2 неверный вызов, 3 плохие входные данные, 4 нет
  ресурса, 5 результат не сохранился.
- Веса уже внутри образов, монтировать не нужно: rfd4 в `/weights`, остальные в `/app/weights`.
- **rfd4 принимает только CIF и BCIF, PDB отклоняется.**
- У кофолдинга по умолчанию `missing_msas=fail`: запрос без `msa_path` и без совпадения в
  `--msa-dir` будет отвергнут.

## Команды

Генерация дизайнов:

```bash
mkdir -p in out
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/in:/input:ro" -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07 \
  module-run design --input /input/design.yaml --output /output/run \
  --model her2-conservative --num-replicates 10 --seed 42
```

Отбор лучших и подготовка запросов на переcворачивание:

```bash
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07 \
  module-run refold-prepare --input /output/run --output /output/refold --top-n 50
```

Кофолдинг (Protenix v2, по данным лучший на антителах):

```bash
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2 \
  module-run predict --input /output/refold/queries --output /output/pred_protenix \
  --seeds 42 43 44 --num-diffusion-samples 5 --write-detailed-scores true
```

OpenFold3 (обязателен для денежного приза, плюс второй голос в консенсусе):

```bash
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-openfold3-v1.3.0-of3ob \
  module-run predict --input /output/refold/queries --output /output/pred_of3 --seeds 42
```

Скоринг переcворачивания:

```bash
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07 \
  module-run refold-score --input /output/refold --predictions /output/pred_protenix \
  --output /output/score_protenix --model protenix-v2 --criterion antibody
```

MSA (нужен URL сервера, уточнить у организаторов):

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/in:/input:ro" -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-msa-v1.1.0 \
  module-run fetch-msa --input /input/seqs.fasta --output /output/msa \
  --input-format fasta --server.type colabfold --server.url <URL>
```

Метрики структур без GPU:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/in:/input:ro" -v "$PWD/out:/output" \
  quay.io/apheris/oslo-hackathon:apheris-data-v1.2.0 \
  module-run compute-metrics --input /input/manifest.yaml --output /output/metrics
```

## Параметры design, которые важны

| Флаг | Диапазон | По умолчанию | Смысл |
|---|---|---|---|
| `--num-replicates` | 1-200 | 1 | дизайнов на каждый вход. Всего = входы × replicates |
| `--batch-size` | ≥ 1 | 1 | входов в батче |
| `--seed` | ≥ 0 | случайный | **без значения сид случайный, для воспроизводимости задавать** |
| `--recycle-strategy` | `thinking_low/medium/high`, либо `low|medium|high` × `flat|front_loaded|midpeak` | `thinking_medium` | бюджет recycling. Для VHH upstream берёт `high/front_loaded` |
| `--resample-hotspots` | bool | false | выбирает 2-3 экспонированных атома рядом с hotspot на каждый дизайн |
| `--partial-trajectory-fraction` | (0, 1] | выкл | частичная диффузия из готового `*_generated.cif` |
| `--inverse-folding` | bool | false | новые последовательности на той же CA-трассе |
| `--sequence-tsr-k` | > 0 | 4 | резкость сэмплинга, переопределяет пресет `--model` |
| `--compile` | bool | false | окупается примерно от 8 дизайнов на длину |

Отдельных полей «длина CDR» и «формат антитела» нет. Длина задаётся expandable-сегментами
в `contig`, формат задаётся входным каркасом и условием `C_CATH` (для иммуноглобулинового
фолда это `"2.60.40"`).

## Формат condition spec

Ключ верхнего уровня это имя дизайна. Поля: `input` (путь к CIF, относительный от файла
spec), `input_as_annotated_cif` (true для частичной диффузии), `contig`, `conditions`.
Зарезервированный ключ `global_args` задаёт значения по умолчанию.

Токены `contig`: `/0` разрыв цепи, `A42-48` сохраняемый сегмент с исходной нумерацией,
`E1:2-20` именованный expandable-сегмент, `20-30` безымянный expandable.

Условия, которые нам нужны:

| Аннотация | Значения | Уровень | Смысл |
|---|---|---|---|
| `C_HOT` | true, false | атомы | **хотспоты, то есть эпитоп** |
| `C_CRD` | true, false | атомы | фиксировать координаты |
| `C_SEQ` | true, false | остатки | фиксировать последовательность |
| `C_CRDGRP` | целое или null | атомы | жёсткая группа без абсолютных координат |
| `C_CATH` | `"2.60.40"` | остатки или цепь | CATH-метка фолда |
| `C_RAS` | BURIED, INTERMEDIATE, EXPOSED | атомы | доступность |

Пример из образа (`condition_spec_ppi.yaml`, дизайн связывателя PD-L1):

```yaml
pdl1_binder:
  input: ../../../examples/inputs/5o45_cropped.cif.zst
  contig: "50-120, /0, A17-131"
  conditions:
    C_CRDGRP:
      1:
        - select: "sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()"
    C_SEQ:
      true:
        - select: "sel('A')"
    C_HOT:
      true:
        - select: "sel('A/*/56/[CG,OH]')"
        - select: "sel('A/*/115/[CG,SD]')"
        - select: "sel('A/*/123/[CD2,OH]')"
```

Синтаксис выбора: `sel('цепь/*/остаток/[атомы]')`. Хотспоты задаются **на уровне атомов**,
не остатков, то есть из карты эпитопов надо брать конкретные атомы боковых цепей.

Готовые примеры антител внутри образа: `/opt/rfd4/examples/antibody_design/` с
`vhh.yaml`, `partial.yaml` и 8 структурами в `structures/vhh/*.cif.zst` (мишени PMVK и
MZB1, четыре каркаса). Их можно подавать в `--input` прямо из контейнера.
Туториал: `/opt/rfd4/docs/examples/antibody_design_tutorial.ipynb`.

## Developability: пороги для VHH

**TNP** (OPIG) это единственный инструмент, откалиброванный на нанотелах, по 36 клиническим.
TAP рассчитан на пары VH и VL и к одиночному VHH не применим.

| Метрика TNP | amber | red |
|---|---|---|
| суммарная длина CDR | 20-24 и 37.25-39 | < 20 или > 39 |
| длина CDR3 | 5-8 и 21.25-23 | < 5 или > 23 |
| компактность CDR3 | 0.56-0.81 и 1.57-1.61 | < 0.56 или > 1.61 |
| PSH (гидрофобные патчи CDR) | 73.40-79.59 и 126.83-155.47 | < 73.40 или > 155.47 |
| PPC (положительный заряд патчей) | 0.39-1.18 | > 1.18 |
| PNC (отрицательный заряд патчей) | 1.47-1.88 | > 1.88 |

Amber и red это 5% хвосты распределения клинических нанотел и выход за наблюдаемый
клинический диапазон. Это не пороги годности, а ориентиры.

Установка: `conda install -c conda-forge -c bioconda tnp`, запуск
`TNP --output <dir> --file sequences.fasta`. Веса NanoBodyBuilder2 скачиваются при первом
запуске, поэтому прогнать заранее.

Прочее: instability index < 40 (ExPASy), AbNatiV2 nativeness ≥ 0.8 (`abnativ score -nat
VHH2 -i file.fa -align`, лицензия CC BY-NC-SA, только некоммерческое использование).

Порогов для pI, GRAVY и суммарного заряда VHH в литературе нет, направление по pI
противоречиво. Записывать, но не назначать порог.

## Мотивы, которых избегать в CDR

По таблице LAP (PLOS Comput Biol 2024), уровни из самой таблицы:

| Мотив | Серьёзность | Что происходит |
|---|---|---|
| N[GS] (NG, NS) | высокая | дезамидирование, NG около суток при pH 7.4 и 37 °C |
| D[DGHST] (DG, DS, DT, DH, DD) | высокая | изомеризация через сукцинимид |
| DP | высокая | кислотно-лабильный разрыв цепи |
| N[^P][ST] | высокая | N-гликозилирование, гликан в CDR мешает связыванию |
| непарный Cys | высокая | свободный тиол, агрегация |
| M, W | средняя | окисление |
| N[AHNT], NP, TS | средняя | дезамидирование и фрагментация |
| GPR, RGD, RYD, LDV, DGE, KGD, NGR | низкая | мотивы связывания интегринов |

В среднем на одиночную цепь приходится 3-3.4 таких мотива, поэтому цель «ноль»
нереалистична. Сравнивать с референсом, а не с нулём.

Из изомеризации D следует предпочитать E вместо D в CDR.

## Полиреактивность

Из работы по нанотелам (Nat Commun 2022), одиночные признаки по AUC:

| Признак | AUC |
|---|---|
| pI | 0.79 |
| число Arg | 0.73 |
| длина CDR3 | 0.58 |
| индекс гидрофобности | 0.57 |

То есть **заряд предсказывает полиреактивность сильнее, чем гидрофобность**. Arg в любом
CDR повышает, Lys, Trp или Tyr в CDR3 тоже. Кислые остатки в CDR2 и CDR3 снижают.
Числовых порогов по числу Arg в литературе нет.

## Нумерация VHH по IMGT

| Регион | Позиции | Примечание |
|---|---|---|
| FR1 | 1-26 | первый Cys = 23 |
| CDR1 | 27-38 | |
| FR2 | 39-55 | консервативный Trp = 41 |
| CDR2 | 56-65 | |
| FR3 | 66-104 | второй Cys = 104 |
| CDR3 | 105-117 | |
| FR4 | 118-129 | WGQG |

Hallmark FR2: IMGT 42, 49, 50, 52 (Kabat 37, 44, 45, 47). В VH это V, G, L, W, в
классическом VHH становится F/Y, E, R, L/G/F. Эти позиции закрывают поверхность, которая
у VH была интерфейсом с лёгкой цепью, поэтому их возврат к гидрофобным делает домен
нерастворимым без лёгкой цепи.

Канонический дисульфид: Cys 23 и Cys 104 по IMGT. Дополнительный мост CDR1-CDR3
встречается примерно у четверти VHH, но эффект в литературе противоречив, и **намеренно
вводить его не стоит**, тем более что у Adaptyv бесклеточная бактериальная экспрессия.
Если мост есть в шаблоне, сохранить.

Длина CDR3 у классических VHH альпаки около 15.95 ± 0.01 остатков, у неклассических
13.93 ± 0.02.

## Что осталось выяснить

1. URL сервера MSA. Без него кофолдинг с `missing_msas=fail` не пройдёт. Либо получить
   URL, либо `--missing-msas ignore`, либо считать MSA заранее и подавать `--msa-dir`.
2. Противоречие в документации rfd4: `yaml-json.md` допускает только `.cif`, `.cif.gz`,
   `.bcif`, а примеры используют `.cif.zst`. Проверить на первом запуске.
3. Лицензия на веса Protenix v2 (в образе файла лицензии нет).
4. GPU один на команду. Договориться об очереди до первого запуска.
