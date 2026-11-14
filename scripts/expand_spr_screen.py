"""Add SPR-prioritized proposals missed by the earlier FACS-only prescreen."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
from scripts.rank_spr_candidates import predict_from_parameters
from src.her2_candidates import WT_CDRH3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pool', type=Path, required=True)
    parser.add_argument('--previous-preflight', type=Path, required=True)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--facs', type=Path, required=True)
    parser.add_argument('--spr', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--top-n', type=int, default=250)
    args = parser.parse_args()
    if args.top_n <= 0:
        raise ValueError('top-n must be positive')
    model = json.loads(args.model.read_text())
    pool = pd.read_csv(args.pool)
    previous = json.loads(args.previous_preflight.read_text())
    if not pool.cdrh3.is_unique:
        raise ValueError('Pool has duplicate CDRH3s')
    pool = pool[pool.wt_hamming_distance >= 4].copy()
    pool['spr_prescreen_score'] = predict_from_parameters(model, pool.cdrh3)
    chosen = pool.sort_values(['spr_prescreen_score', 'cdrh3'], ascending=[False, True]).head(args.top_n)
    chosen = chosen[~chosen.cdrh3.isin([row['cdrh3'] for row in previous])]
    args.output.mkdir(parents=True, exist_ok=False)
    input_csv = args.output / 'additional_candidates.csv'
    input_json = args.output / 'additional_preflight.json'
    chosen.to_csv(input_csv, index=False)
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(root / 'src/her2_preflight.py'),
                    '--candidates', str(input_csv.resolve()), '--facs', str(args.facs.resolve()),
                    '--spr', str(args.spr.resolve()), '--output', str(input_json.resolve())],
                   check=True, env={**os.environ, 'PYTHONPATH': str(root / 'src')})
    extra = json.loads(input_json.read_text())
    combined = previous + extra
    assert len({row['cdrh3'] for row in combined}) == len(combined)
    (args.output / 'combined_preflight.json').write_text(json.dumps(combined, indent=2) + '\n')
    report = {
        'purpose': 'post_benchmark_exploratory_prescreen_not_new_performance_estimate',
        'pool_with_at_least_four_wt_changes': len(pool), 'spr_top_n': args.top_n,
        'additional_screened': len(extra),
        'additional_strict_passes': sum(row['strict_local_gate'] for row in extra),
        'combined_screened': len(combined),
        'combined_strict_passes': sum(row['strict_local_gate'] for row in combined),
        'wt_score_diagnostic': float(predict_from_parameters(model, pd.Series([WT_CDRH3]))[0]),
        'score_caveat': 'WT was excluded from training; absolute KD calibration is not established.',
        'inputs': {name: {'source': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                   for name, path in [('pool', args.pool), ('previous_preflight', args.previous_preflight),
                                      ('model', args.model), ('facs', args.facs), ('spr', args.spr)]}}
    (args.output / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main()
