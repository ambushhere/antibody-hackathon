"""Post-benchmark distance diagnostic of frozen SPR predictions; no model fitting."""
import argparse
import json
from pathlib import Path

from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from scripts import rank_spr_candidates as rank
from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr

MODELS = (facs.SCORE, 'position_ridge')


def audit(data_dir, oof_path, metrics_path):
    benchmark, cohort = rank.load_benchmark(metrics_path)
    primary, inputs, splits = rank.load_primary(data_dir, benchmark)
    rows, oof_input = rank.load_oof(oof_path, MODELS)
    _, aligned = rank.original_folds(primary, rows, cohort)
    rank.max_gap(primary[spr.TARGET], aligned[spr.TARGET], f'OOF {spr.TARGET}', rank.OOF_TOLERANCE)
    if not aligned[facs.SCORE].between(0, 1).all():
        raise ValueError('OOF FACS scores must be between 0 and 1')

    facs_training_sources = list(splits['train'] | splits['val'])
    rows['nearest_facs_trainval'] = [
        process.extractOne(seq, facs_training_sources, scorer=Levenshtein.distance)[1]
        for seq in rows.cdrh3
    ]
    rows['nearest_spr_outer_train'] = [
        process.extractOne(row.cdrh3, rows.loc[rows.fold != row.fold, 'cdrh3'].tolist(), scorer=Levenshtein.distance)[1]
        for row in rows.itertuples()]
    if not rows.nearest_spr_outer_train.gt(spr.MAX_DISTANCE).all():
        raise ValueError('Primary OOF rows are too close to their outer SPR training fold')

    return rows, summarize(rows, [oof_input, *inputs.values(), rank.source(metrics_path)], benchmark['seed'])


def summarize(rows, inputs, seed):
    output = {'scope': 'post_benchmark_diagnostic_no_tuning', 'subsets': {},
              'hashes': {entry['source']: entry['sha256'] for entry in inputs}}
    for name, frame in [('primary', rows), ('distance_at_least_3_to_both_training_sources', rows[rows.nearest_facs_trainval >= 3])]:
        y = frame[spr.TARGET].to_numpy()
        scores = {model: frame[model].to_numpy() for model in MODELS}
        intervals, deltas = spr.group_bootstrap(y, scores, {}, frame.group.to_numpy(), seed)
        output['subsets'][name] = {'n': len(frame), 'groups': int(frame.group.nunique()),
                                   'spearman': {model: spr.spearman(y, values) for model, values in scores.items()},
                                   'intervals': intervals, 'paired_deltas': deltas}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--oof', type=Path, default=Path('results/spr-benchmark/oof_predictions.csv'))
    parser.add_argument('--metrics', type=Path, default=Path('results/spr-benchmark/metrics.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows, output = audit(args.data_dir, args.oof, args.metrics)
    args.output.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.output / 'distance_rows.csv', index=False)
    (args.output / 'report.json').write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: {x: v for x, v in value.items() if x in ['n', 'groups', 'spearman', 'paired_deltas']} for k, value in output['subsets'].items()}, indent=2))


if __name__ == '__main__':
    main()
