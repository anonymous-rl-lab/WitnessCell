#!/usr/bin/env python3
"""Audit frozen Dixit RNA arithmetic from archived per-pair/fold error summaries.

This does not recompute predictions, endpoint projections or the original fit.
The separate REPRODUCE.py refits the later G2M development experiment.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    source = ROOT / 'primary_rna_archive'
    for rel, expected in json.loads((source / 'SOURCE_HASHES.json').read_text()).items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == expected, rel
    evaluation = json.loads((source / 'EVALUATION.json').read_text())
    decision = json.loads((source / 'DECISION.json').read_text())
    methods = evaluation['methods']
    pairs = list(csv.DictReader((source / 'PER_PAIR.csv').open()))
    assert len({row['pair'] for row in pairs}) == 45
    fields = ['full_sse', 'projected_sse', 'endpoint_sse',
              'standardized_full_sse', 'standardized_projected_sse',
              'program_full_sse', 'program_projected_sse']
    max_relative_error = 0.0
    for name, scores in methods.items():
        for field in fields:
            total = sum(fold['methods'][name][field] for fold in evaluation['fold_results'])
            error = abs(total - scores[field]) / max(1.0, abs(scores[field]))
            max_relative_error = max(max_relative_error, error)
            assert error < 2e-12, (name, field, total, scores[field])
        mean_mse = sum(float(row['non_target_mse']) for row in pairs if row['method'] == name) / 45
        assert abs(mean_mse - scores['full_mse']) < 2e-14, name
        assert abs(scores['full_sse'] - scores['projected_sse'] - scores['endpoint_sse']) < 2e-10
    def gain(method, reference, field):
        return 100 * (1 - methods[method][field] / methods[reference][field])
    values = {
        'complete_core_vs_additive_full_gain_percent': gain('zero', 'additive', 'full_sse'),
        'selected_vs_core_full_gain_percent': gain('selected', 'zero', 'full_sse'),
        'selected_vs_core_partner_gain_percent': gain('selected', 'zero', 'projected_sse')}
    assert abs(values['selected_vs_core_full_gain_percent'] - decision['secondary_full_gain_percent']) < 1e-11
    assert abs(values['selected_vs_core_partner_gain_percent'] - decision['primary_projected_gain_percent']) < 1e-11
    assert values['selected_vs_core_partner_gain_percent'] < 0
    report = {'status': 'PASS', 'scientific_decision': 'PRIMARY_RNA_PARTNER_NOT_CONFIRMED',
              'source': evaluation['source'], 'pairs': 45, 'non_target_readouts': 23519,
              'methods_checked': len(methods), 'folds_checked': 3,
              'max_fold_sum_relative_difference': max_relative_error, **values,
              'reproduction_level': 'Per-pair full-MSE and fold-SSE arithmetic audit only.',
              'limitations': 'Original full-RNA prediction arrays, projection refit and full source runner are not included here. G2M REPRODUCE.py is a distinct, later development fit and cannot replace this frozen primary outcome.'}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
