#!/usr/bin/env python3
"""Frozen-feature seen ablation: separate scoring and evaluation commands."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import hashlib
import json
import math
import shutil
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PLAN_PATH = HERE.parent / 'analysis_plan.json'
PLAN = json.loads(PLAN_PATH.read_text())
CFG = PLAN['seen_comparison']
METHODS = list(CFG['methods'])
SOURCE = Path('/workspace/scratch/7517238e1d04/v19_inputs/step1_unpacked/witnesscell_step1_run')
EXPECTED = {'Norman': [73, 66, 69], 'Wessels': [63, 54, 64]}


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError('No records for ' + str(path))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def tie(dataset, pair):
    return hashlib.sha256((CFG['tie_protocol_id'] + '|' + dataset + '|' + pair).encode()).hexdigest()


def quota(n, coverage):
    frac = Fraction(str(coverage))
    return int(n) * frac.numerator // frac.denominator


def feature(row, name):
    if name == 'seen0':
        return float(int(row['seen']) == 0)
    if name == 'seen2':
        return float(int(row['seen']) == 2)
    return float(row[name])


def fit_ridge(rows, names):
    weights = np.array([float(r['weight']) for r in rows])
    weights /= weights.sum()
    losses = np.array([float(r['loss_mse']) for r in rows])
    order = np.argsort(losses, kind='stable')
    mid = min(np.searchsorted(np.cumsum(weights[order]), .5 * weights.sum(), side='left'), len(order)-1)
    median = max(float(losses[order][mid]), 1e-12)
    y = np.log1p(losses / median)
    x = np.array([[feature(r, n) for n in names] for r in rows])
    mean = weights @ x
    sd = np.sqrt(weights @ ((x - mean) ** 2))
    sd[sd == 0] = 1.
    x = (x - mean) / sd
    intercept = float(weights @ y)
    beta = np.linalg.solve(x.T @ (weights[:, None] * x) + np.eye(len(names)),
                           x.T @ (weights * (y - intercept)))
    return {'features': names, 'mean': mean.tolist(), 'sd': sd.tolist(),
            'coef': beta.tolist(), 'intercept': intercept, 'm': median}


def predict(model, row):
    x = np.array([feature(row, n) for n in model['features']])
    return float(((x - model['mean']) / model['sd']) @ model['coef'] + model['intercept'])


def order_rows(rows, method, dataset):
    if method == 'R':
        key = lambda i: (-int(rows[i]['seen']), float(rows[i]['risk_B']), tie(dataset, rows[i]['pair']), rows[i]['pair'])
    else:
        key = lambda i: (float(rows[i]['risk_' + method]), tie(dataset, rows[i]['pair']), rows[i]['pair'])
    return np.array(sorted(range(len(rows)), key=key), dtype=int)


def prepare(source):
    if (HERE / 'all_new_decisions_frozen.json').exists():
        raise ValueError('A completed score freeze exists; do not overwrite it.')
    start = time.monotonic()
    source_inventory, replay_checks, freezes = [], [], []
    dump(HERE / 'analysis_execution_lock.json', {
        'created_utc': utc(), 'plan_sha256': sha(PLAN_PATH), 'script_sha256': sha(__file__),
        'historical_access': PLAN['historical_data_access'],
        'this_stage_reads_evaluation_losses': False, 'source_root': str(source)})
    for dataset in CFG['datasets']:
        for seed in CFG['seeds']:
            src = source / 'results' / f'{dataset}_seed{seed}'
            out = HERE / 'results' / f'{dataset}_seed{seed}'
            out.mkdir(parents=True, exist_ok=True)
            inputs = {}
            for name in ('training_oof_records.csv', 'risk_models.json', 'frozen_query_outputs.csv',
                         'evaluation_metadata.json', 'split_manifest.json', 'frozen_selections.csv'):
                target = out / ('original_' + name)
                shutil.copyfile(src / name, target)
                inputs[name] = {'sha256': sha(target), 'relative_copy': str(target.relative_to(HERE))}
            source_inventory.append({'dataset': dataset, 'seed': seed, 'files': inputs})
            oof = read_csv(out / 'original_training_oof_records.csv')
            counts = Counter(r['pair'] for r in oof)
            for row in oof:
                if not math.isclose(float(row['weight']), 1 / counts[row['pair']], rel_tol=1e-12):
                    raise ValueError('OOF pair weights changed')
                if int(row['missing_endpoints']) != 2 - int(row['seen']):
                    raise ValueError('Inconsistent seen encoding')
            original_models = json.loads((out / 'original_risk_models.json').read_text())
            original_rows = read_csv(out / 'original_frozen_query_outputs.csv')
            if len(original_rows) != EXPECTED[dataset][seed - 1]:
                raise ValueError('Unexpected candidate count')
            if any('loss_mse' in r or 'true_effect' in r for r in original_rows):
                raise ValueError('Evaluation labels present in scoring input')
            for name in ('B', 'S'):
                check = fit_ridge(oof, original_models[name]['features'])
                max_score_difference = max(abs(predict(check, r) - float(r['risk_' + name])) for r in original_rows)
                if max_score_difference > 1e-12:
                    raise ValueError('Original risk-head replay mismatch')
                replay_checks.append({'dataset': dataset, 'seed': seed, 'method': name,
                                      'max_abs_score_difference': max_score_difference})
            models = {
                'C': fit_ridge(oof, ['effect_score', 'variance_score', 'seen0', 'seen2']),
                'C0': fit_ridge(oof, ['effect_score', 'variance_score', 'missing_endpoints'])}
            dump(out / 'new_risk_models.json', models)
            scores = []
            for r in original_rows:
                s = dict(r)
                s.update({('risk_' + name): predict(model, r) for name, model in models.items()})
                scores.append(s)
            write_csv(out / 'new_frozen_scores.csv', scores)
            decisions = []
            old_decisions = {(float(r['coverage']), r['method'], r['pair']): int(r['accepted'])
                             for r in read_csv(out / 'original_frozen_selections.csv')}
            for scope in CFG['scopes']:
                pool = scores if scope == 'all_queries' else [r for r in scores if int(r['seen']) == 1]
                orders = {m: order_rows(pool, m, dataset) for m in METHODS}
                if scope == 'seen1' and not np.array_equal(orders['R'], orders['B']):
                    raise AssertionError('Within seen1 the R order must equal B')
                for coverage in CFG['coverages']:
                    k = quota(len(pool), coverage)
                    for method in METHODS:
                        chosen = set(orders[method][:k].tolist())
                        for i, row in enumerate(pool):
                            accepted = int(i in chosen)
                            if scope == 'all_queries' and method in ('S', 'B'):
                                if accepted != old_decisions[(coverage, method, row['pair'])]:
                                    raise AssertionError('Original S/B decisions changed')
                            decisions.append({'dataset': dataset, 'seed': seed, 'scope': scope,
                                              'pair': row['pair'], 'coverage': coverage, 'method': method,
                                              'accepted': accepted, 'quota': k, 'candidate_count': len(pool)})
            write_csv(out / 'new_frozen_decisions.csv', decisions)
            freezes.append({'dataset': dataset, 'seed': seed,
                            'scores_sha256': sha(out / 'new_frozen_scores.csv'),
                            'decisions_sha256': sha(out / 'new_frozen_decisions.csv'),
                            'models_sha256': sha(out / 'new_risk_models.json')})
    dump(HERE / 'source_inventory.json', source_inventory)
    dump(HERE / 'original_risk_replay_checks.json', replay_checks)
    dump(HERE / 'all_new_decisions_frozen.json', {
        'frozen_utc': utc(), 'script_sha256': sha(__file__), 'plan_sha256': sha(PLAN_PATH),
        'source_inventory_sha256': sha(HERE / 'source_inventory.json'),
        'read_evaluation_losses_in_this_stage': False, 'elapsed_seconds': time.monotonic() - start,
        'batches': freezes})
    print(json.dumps({'stage': 'prepare', 'status': 'FROZEN', 'batches': len(freezes),
                      'seconds': time.monotonic() - start}), flush=True)


def interval(values):
    values = np.asarray(values)
    good = values[np.isfinite(values)]
    return {'ci95': np.quantile(good, [.025, .975]).tolist() if len(good) else None,
            'valid': int(len(good)), 'invalid': int(len(values) - len(good))}


def geometric_ratio(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    valid = np.all(np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b > 0), axis=0)
    out = np.full(a.shape[1:], np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        vals = np.exp(np.mean(np.log(a / b), axis=0))
    return np.where(valid, vals, out)


def selection_weights(counts, order, coverage):
    multiplicity = np.asarray(counts)[:, order]
    frac = Fraction(str(coverage))
    k = multiplicity.sum(axis=1) * frac.numerator // frac.denominator
    before = np.cumsum(multiplicity, axis=1) - multiplicity
    selected = np.minimum(multiplicity, np.maximum(k[:, None] - before, 0))
    if not np.allclose(selected.sum(axis=1), k, atol=1e-10, rtol=0):
        raise AssertionError('Bootstrap quota mismatch')
    return selected, k


def risks(rows, counts, coverage, dataset, field):
    vals = np.array([float(r[field]) for r in rows])
    out = {}
    for method in METHODS:
        order = order_rows(rows, method, dataset)
        weights, k = selection_weights(counts, order, coverage)
        numerator = weights @ vals[order]
        out[method] = np.divide(numerator, k, out=np.full(len(k), np.nan), where=k > 0)
    return out


def evaluate(source):
    t = time.monotonic()
    freeze = json.loads((HERE / 'all_new_decisions_frozen.json').read_text())
    if sha(__file__) != freeze['script_sha256'] or sha(PLAN_PATH) != freeze['plan_sha256']:
        raise ValueError('Frozen analysis source or plan changed')
    evaluation_path = HERE / 'evaluation'
    if (evaluation_path / 'results.json').exists():
        raise ValueError('Completed evaluation exists; refuse overwrite')
    evaluation_path.mkdir(parents=True, exist_ok=True)
    dump(evaluation_path / 'evaluation_started.json', {
        'started_utc': utc(), 'freeze_sha256': sha(HERE / 'all_new_decisions_frozen.json'),
        'historical_access': PLAN['historical_data_access']})
    batches = {}
    for item in freeze['batches']:
        dataset, seed = item['dataset'], item['seed']
        out = HERE / 'results' / f'{dataset}_seed{seed}'
        for key, filename in [('scores_sha256', 'new_frozen_scores.csv'),
                              ('decisions_sha256', 'new_frozen_decisions.csv'),
                              ('models_sha256', 'new_risk_models.json')]:
            if sha(out / filename) != item[key]:
                raise ValueError('Frozen artifact changed')
        shutil.copyfile(source / 'results' / f'{dataset}_seed{seed}' / 'scored_queries.csv',
                        out / 'original_scored_queries.csv')
        old = {r['pair']: r for r in read_csv(out / 'original_scored_queries.csv')}
        rows = read_csv(out / 'new_frozen_scores.csv')
        for row in rows:
            previous = old[row['pair']]
            for f in ('seen', 'effect_score', 'variance_score', 'risk_B', 'risk_S'):
                if float(previous[f]) != float(row[f]):
                    raise ValueError('Score/label join mismatch')
            for field in CFG['losses'] + ['true_effect']:
                row[field] = float(previous[field])
        write_csv(out / 'new_scored_queries.csv', rows)
        batches[(dataset, seed)] = {'rows': rows,
            'metadata': json.loads((out / 'original_evaluation_metadata.json').read_text())}
    all_results, seed_rows, diagnostics, bootstrap_arrays, effect_results = [], [], [], {}, []
    for ds_index, dataset in enumerate(CFG['datasets']):
        ids = sorted({r['pair'] for s in CFG['seeds'] for r in batches[(dataset, s)]['rows']})
        index = {p: i for i, p in enumerate(ids)}
        rng = np.random.default_rng(CFG['bootstrap']['seed'] + ds_index)
        draws = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)),
                               size=CFG['bootstrap']['replicates'])
        dataset_counts = {}
        for seed in CFG['seeds']:
            rows = batches[(dataset, seed)]['rows']
            dataset_counts[seed] = np.vstack([np.ones(len(rows), dtype=int),
                                              draws[:, [index[r['pair']] for r in rows]]])
            for method in METHODS:
                order = order_rows(rows, method, dataset)
                chosen = set(order[:quota(len(rows), .5)].tolist())
                for seen in (0, 1, 2):
                    pool_ids = [i for i, r in enumerate(rows) if int(r['seen']) == seen]
                    accepted = [i for i in pool_ids if i in chosen]
                    diagnostics.append({'dataset': dataset, 'seed': seed, 'method': method,
                        'seen': seen, 'candidate_count': len(pool_ids), 'accepted_count': len(accepted),
                        'coverage': len(accepted)/len(pool_ids) if pool_ids else None,
                        'accepted_mse': float(np.mean([rows[i]['loss_mse'] for i in accepted])) if accepted else None,
                        'accepted_true_effect': float(np.mean([rows[i]['true_effect'] for i in accepted])) if accepted else None})
        for scope in CFG['scopes']:
            for coverage in CFG['coverages']:
                for field in CFG['losses']:
                    by_method = {m: [] for m in METHODS}
                    scope_info = []
                    for seed in CFG['seeds']:
                        full_rows = batches[(dataset, seed)]['rows']
                        mask = np.array([scope == 'all_queries' or int(r['seen']) == 1 for r in full_rows])
                        # Secondary losses share one method-common finite pool.
                        mask &= np.array([math.isfinite(float(r[field])) for r in full_rows])
                        rows = [r for i, r in enumerate(full_rows) if mask[i]]
                        counts = dataset_counts[seed][:, mask]
                        loss_by_method = risks(rows, counts, coverage, dataset, field)
                        scope_info.append({'seed': seed, 'n': len(rows), 'quota': quota(len(rows), coverage)})
                        for method in METHODS:
                            by_method[method].append(loss_by_method[method])
                            seed_rows.append({'dataset': dataset, 'seed': seed, 'scope': scope,
                                'coverage': coverage, 'loss': field, 'method': method,
                                'n': len(rows), 'quota': quota(len(rows), coverage),
                                'mean_selected_loss': float(loss_by_method[method][0])})
                    for num, den in CFG['comparisons']:
                        ratio = geometric_ratio(by_method[num], by_method[den])
                        info = interval(ratio[1:])
                        row = {'dataset': dataset, 'scope': scope, 'coverage': coverage,
                               'loss': field, 'comparison': num + '/' + den,
                               'ratio': float(ratio[0]), 'gain': float(1-ratio[0]),
                               'ratio_ci95': info['ci95'],
                               'gain_ci95': [1-info['ci95'][1], 1-info['ci95'][0]] if info['ci95'] else None,
                               'valid_bootstrap': info['valid'], 'invalid_bootstrap': info['invalid'],
                               'seed_counts': scope_info,
                               'seed_ratios': {str(s): float(by_method[num][i][0]/by_method[den][i][0])
                                               for i, s in enumerate(CFG['seeds'])}}
                        all_results.append(row)
                        bootstrap_arrays[f'{dataset}__{scope}__{coverage}__{field}__{num}_{den}'] = ratio[1:]
        # Effect-size control uses original training-only cutpoints and original population weighting.
        aggregate = {m: [] for m in ('S', 'C')}
        for seed in CFG['seeds']:
            rows = batches[(dataset, seed)]['rows']
            cuts = batches[(dataset, seed)]['metadata']['oof_effect_quartile_edges']
            quartiles = np.searchsorted(cuts, [float(r['effect_score']) for r in rows], side='right')
            counts = dataset_counts[seed]
            totals = counts.sum(axis=1)
            numerators = {m: np.zeros(len(totals)) for m in ('S', 'C')}
            valid = np.ones(len(totals), dtype=bool)
            for q in range(4):
                mask = quartiles == q
                subset = [r for i, r in enumerate(rows) if mask[i]]
                nc = counts[:, mask].sum(axis=1)
                rr = risks(subset, counts[:, mask], .5, dataset, 'loss_mse')
                for m in ('S', 'C'):
                    valid &= (nc == 0) | np.isfinite(rr[m])
                    numerators[m] += np.where(nc == 0, 0., nc * rr[m])
            for m in ('S', 'C'):
                aggregate[m].append(np.where(valid & (totals > 0), numerators[m]/totals, np.nan))
        eratio = geometric_ratio(aggregate['S'], aggregate['C'])
        ci = interval(eratio[1:])
        effect_results.append({'dataset': dataset, 'comparison': 'S/C',
            'scope': 'all_queries_within_original_TRAIN_OOF_effect_quartiles',
            'coverage_per_quartile': .5, 'ratio': float(eratio[0]), 'gain': float(1-eratio[0]),
            'gain_ci95': [1-ci['ci95'][1], 1-ci['ci95'][0]] if ci['ci95'] else None,
            'valid_bootstrap': ci['valid'], 'invalid_bootstrap': ci['invalid']})
        bootstrap_arrays[dataset + '__effect_quartile_control_S_C'] = eratio[1:]
        print(json.dumps({'stage': 'evaluation', 'dataset_complete': dataset,
                          'elapsed_seconds': time.monotonic()-t}), flush=True)
    summary = {'protocol_id': PLAN['protocol_id'], 'completed_utc': utc(),
               'design': PLAN['study_status'], 'results': all_results,
               'effect_quartile_control': effect_results, 'elapsed_seconds': time.monotonic()-t,
               'bootstrap': CFG['bootstrap'], 'freeze_sha256': sha(HERE / 'all_new_decisions_frozen.json')}
    dump(evaluation_path / 'results.json', summary)
    csv_results = []
    for row in all_results:
        csv_results.append({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()})
    write_csv(evaluation_path / 'comparisons.csv', csv_results)
    write_csv(evaluation_path / 'seed_results.csv', seed_rows)
    write_csv(evaluation_path / 'seen_allocation.csv', diagnostics)
    np.savez_compressed(evaluation_path / 'bootstrap_ratios.npz', **bootstrap_arrays)
    primary = [r for r in all_results if r['coverage'] == .5 and r['loss'] == 'loss_mse'
               and r['comparison'] in ('S/C', 'S/R', 'S/B', 'S/C0')]
    dump(evaluation_path / 'primary_and_key_comparisons.json', primary)
    print(json.dumps({'stage': 'evaluate', 'status': 'COMPLETE',
                      'seconds': time.monotonic()-t, 'key_results': primary}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'evaluate'])
    parser.add_argument('--source', type=Path, default=SOURCE)
    args = parser.parse_args()
    {'prepare': prepare, 'evaluate': evaluate}[args.stage](args.source)
