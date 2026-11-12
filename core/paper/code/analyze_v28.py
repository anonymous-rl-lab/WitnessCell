#!/usr/bin/env python3
"""Sensitivity of the frozen leaderboard and recovery of a historical protocol.

No estimator fitting or score regeneration: rerank the archived six metric scores.
"""
from pathlib import Path
import csv, gzip, hashlib, io, json, zipfile, os
import numpy as np

R = Path(__file__).resolve().parents[1]
D = R / 'revision'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def min_rank(values, higher_is_better=False):
    v = -np.asarray(values) if higher_is_better else np.asarray(values)
    return 1 + (v[:, None] > v[None, :]).sum(axis=1)

source = D / 'source/benchmark_per_operation.csv.gz'
with gzip.open(source, 'rt') as f:
    rows = [r for r in csv.DictReader(f) if r['DEG'] == '100' and '+' in r['op']]
groups = {}
for r in rows:
    groups.setdefault((r['DataSet'], r['op']), {})[r['method']] = r
assert len(groups) == 540 and all(len(g) == 16 for g in groups.values())
metrics = [('cor','Rank_pcc'), ('mse_score','Rank_mse'),
           ('edistance_score','Rank_edistance'), ('sym_score','Rank_sym'),
           ('was_score','Rank_was'), ('DEG_score','Rank_deg')]
expected = {'none':3.5148148148148146, 'linearModel':3.340740740740741,
            'trainMean':3.3666666666666667, 'linearModel+trainMean':3.1796296296296296}
output = []; operations = []; original_metric_checks = 0
for removed in expected:
    exclude = set() if removed == 'none' else set(removed.split('+'))
    methods = sorted(set(next(iter(groups.values()))) - exclude)
    ranks = []
    for (dataset, op), table in sorted(groups.items()):
        six = np.column_stack([min_rank([float(table[m][c]) for m in methods], True)
                               for c, _ in metrics])
        rank = min_rank(six.mean(axis=1))
        if removed == 'none':
            saved = np.array([[float(table[m][c]) for _, c in metrics] for m in methods])
            assert np.array_equal(six, saved), (dataset, op, 'metric ranks')
            assert np.array_equal(rank, [float(table[m]['Rank']) for m in methods])
            original_metric_checks += int(six.size)
        ranks.append(rank)
        for i, method in enumerate(methods):
            operations.append(dict(removed=removed,dataset=dataset,operation=op,
                                   method=method,mean_metric_rank=float(six[i].mean()),
                                   operation_rank=int(rank[i])))
    means = np.mean(ranks,axis=0); positions = min_rank(means)
    part = [dict(removed=removed,method=m,mean_operation_rank=float(means[i]),
                 position=int(positions[i]),methods=len(methods),operation_records=540)
            for i,m in enumerate(methods)]
    output.extend(sorted(part,key=lambda x:(x['mean_operation_rank'],x['method'])))
    wc = next(x for x in part if x['method']=='WitnessCell')
    assert abs(wc['mean_operation_rank']-expected[removed])<1e-12 and wc['position']==1
write_csv(D/'BASELINE_EXCLUSION_SENSITIVITY.csv',output)
write_csv(D/'BASELINE_EXCLUSION_OPERATION_RANKS.csv',operations)
summary = dict(status='PASS',source_sha256=sha(source.read_bytes()),
    metric_score_columns=[m[0] for m in metrics],higher_scores_are_better=True,
    tie_rule='competition / minimum rank at both metric and operation levels',
    original_metric_ranks_exactly_reconstructed=original_metric_checks,
    original_operation_ranks_exactly_reconstructed=540,
    all_four_WitnessCell_ranks_match_expert_review_within=1e-12,
    conditions=[r for r in output if r['method']=='WitnessCell'],
    interpretation='Sensitivity to comparator inclusion on fixed scores. Mean ranks across different-sized method sets are not performance improvements. Historical score-generation differences remain unresolved.',
    no_training=True)
(D/'BASELINE_EXCLUSION_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')

middle = Path(os.environ.get('WITNESSCELL_LEGACY_ROOT', str(R.parent/'legacy')))/'experiments/26_horlbeck_two_stage/WitnessCell_Horlbeck_TwoStage_Repro_20261018.zip'
inner = 'WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/Jurkat/protocol/v1_before_availability_amendment.zip'
with zipfile.ZipFile(middle) as z: protocol=z.read(inner)
assert len(protocol)==24655 and sha(protocol)=='b93be08507a73a521f8a63fbe88b4db2363ff57c7a38ed6fb4aeca92af6d0e67'
target=R/'reproduction/protocol/v1_before_availability_amendment.zip';target.parent.mkdir(exist_ok=True);target.write_bytes(protocol)
print(json.dumps({'baseline_exclusion':'PASS','original_metric_ranks_checked':original_metric_checks,'protocol_recovered':len(protocol),'conditions':summary['conditions']}))
