"""Мини-задача (репетиция пайплайна) на комплексе антитело-антиген из PDB.

Прогон:
    export SSL_CERT_FILE=$(.venv/bin/python -c "import certifi;print(certifi.where())")
    PYTHONPATH=src .venv/bin/python src/mini_task.py 1N8Z B

Шаги: структура -> цепи -> CDR -> остатки на интерфейсе ->
скан точечных мутаций этих остатков через ESM-2 -> метрики -> results/.
По умолчанию 1N8Z: Fab трастузумаба в комплексе с HER2.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from cdr_tools import mutate
from metrics import liability_flags, score_mutations, sequence_metrics
from structure import cdr_res_ids, disulfide_res_ids, epitope_summary, load_complex

RESULTS = Path(__file__).resolve().parent.parent / "results"
TOP_N = 10


def run(pdb_id: str = "1N8Z", ab_chain: str | None = None) -> pd.DataFrame:
    cx = load_complex(pdb_id)

    print(f"=== {cx.pdb_id}")
    for cid, ch in cx.chains.items():
        print(f"  цепь {cid}: {ch.kind:8} {len(ch.seq)} остатков")

    if ab_chain is None:  # по умолчанию тяжёлая цепь: CDR-H3 определяет специфичность
        heavy = [c for c, ch in cx.antibody_chains.items() if ch.kind == "heavy"]
        ab_chain = heavy[0] if heavy else next(iter(cx.antibody_chains))
    chain = cx.chains[ab_chain]
    if not chain.is_antibody:
        raise SystemExit(f"цепь {ab_chain} не распознана как антитело")

    print(f"\n=== цепь {ab_chain} ({chain.kind})")
    for name, cdr in chain.numbered.cdrs.items():
        print(f"  {name}: {cdr}")

    iface = epitope_summary(cx, ab_chain)
    all_cdr = cdr_res_ids(chain)
    print(f"\n=== остатки CDR в контакте с антигеном (<4.5 Å)")
    for name, rids in iface.items():
        print(f"  {name}: {rids or '—'} из {len(all_cdr[name])}")

    # Номера остатков структуры -> индексы в последовательности цепи
    res_to_idx = {int(r): i for i, r in enumerate(chain.res_ids)}
    hot = sorted({res_to_idx[r] for rids in iface.values() for r in rids if r in res_to_idx})
    if not hot:
        raise SystemExit("на интерфейсе нет остатков CDR — проверьте цепь антигена")
    print(f"\n{len(hot)} позиций для скана, {len(hot) * 19} вариантов, один forward pass ESM-2")

    scores = score_mutations(chain.seq, hot)

    # Цистеины в дисульфидах не считаем «свободными»
    ss = disulfide_res_ids(cx, ab_chain)
    bonded_idx = {i for i, r in enumerate(chain.res_ids) if int(r) in ss}
    if ss:
        print(f"дисульфиды в цепи {ab_chain}: остатки {sorted(ss)}")

    rows = []
    for ms in scores:
        variant = mutate(chain.seq, ms.position, ms.mut)
        res_id = int(chain.res_ids[ms.position])
        rows.append(
            {
                # метка по номеру остатка в структуре, а не по индексу в строке
                "mutation": f"{ms.wt}{res_id}{ms.mut}",
                "res_id": res_id,
                "esm_score": round(ms.esm_score, 3),
                **sequence_metrics(variant),
                "liabilities": "; ".join(
                    liability_flags(variant, (ms.position - 2, ms.position + 3), bonded_idx)
                ),
            }
        )
    df = pd.DataFrame(rows)

    wt = sequence_metrics(chain.seq)
    print(f"\n=== WT ({ab_chain}): " + "  ".join(f"{k}={v}" for k, v in wt.items()))
    print(f"\n=== топ-{TOP_N} по ESM wt-marginal (>0 = правдоподобнее WT)")
    cols = ["mutation", "res_id", "esm_score", "charge_pH7", "gravy", "instability", "liabilities"]
    print(df.head(TOP_N)[cols].to_string(index=False))

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"{cx.pdb_id}_{ab_chain}_scan.csv"
    df.to_csv(out, index=False)
    print(f"\nполный скан: {out.relative_to(RESULTS.parent)}  ({len(df)} вариантов)")
    print("ВАЖНО: ESM-скор — это правдоподобие последовательности, НЕ аффинность.")
    print("Для аффинности нужна структура комплекса варианта или обученный на связывании предиктор.")
    return df


if __name__ == "__main__":
    pdb = sys.argv[1] if len(sys.argv) > 1 else "1N8Z"
    ch = sys.argv[2] if len(sys.argv) > 2 else None
    run(pdb, ch)
