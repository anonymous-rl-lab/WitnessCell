"""Lossless input adapter for program-only follow-up experiments.

Only original condition means and identities are materialized. Cell halves and
cell variances are not needed by these learners and are intentionally absent.
"""
from pathlib import Path
from types import SimpleNamespace
import csv
import numpy as np


def load_task(compact_root, gm, bridge):
    root = Path(compact_root)
    path = root / 'Dixit_program_source_means.npz'
    with np.load(path, allow_pickle=False) as z:
        s = {key: z[key] for key in z.files}
    ci = {str(c): i for i, c in enumerate(s['conditions'])}
    nodes = sorted(c for c in ci if c != 'control' and '+' not in c)
    keep = ~np.isin(s['gene_symbols'].astype(str), nodes)
    genes = s['genes'].astype(str)[keep]
    symbols = s['gene_symbols'].astype(str)[keep]
    labels = root / 'Dixit_labels_only_threefold_candidate.csv'
    with labels.open() as handle:
        rows = list(csv.DictReader(handle))
    pairs = [bridge.canonical(r['gene_a'] + '+' + r['gene_b']) for r in rows]
    foldids = np.array([int(r['fold']) for r in rows])
    blocks = [b for f in range(3) for b in bridge.rectangles(pairs, np.flatnonzero(foldids == f))]
    primary = [[pairs[i] for i in b] for b in blocks]
    used, training_blocks = set(), []
    for block in sorted(blocks, key=lambda ids: gm.digest([pairs[i] for i in ids])):
        if not used.intersection(block):
            training_blocks.append([pairs[i] for i in block])
            used.update(block)
    folds = []
    for f in range(3):
        tr, q = np.flatnonzero(foldids != f), np.flatnonzero(foldids == f)
        inner = gm.preserved_folds(tr, pairs, training_blocks)
        folds.append({'fold': f,
                      'train_pairs': [pairs[i] for i in tr],
                      'test_pairs': [pairs[i] for i in q],
                      'inner_test_pairs': [[pairs[i] for i in ids] for ids in inner],
                      'test_rectangle_blocks': [[pairs[i] for i in b]
                                               for b in bridge.rectangles(pairs, q)]})
    values = {c: s['full_means'][ci[c], keep].copy() for c in ['control', *nodes, *pairs]}
    control = values['control']
    single = {g: values[g] - control for g in nodes}
    d = gm.make_dataset(pairs, single, np.stack([values[p] - control for p in pairs]), training_blocks)
    return SimpleNamespace(name='Dixit', d=d, values=values, genes=genes,
                           symbols=symbols, pairs=pairs, nodes=nodes, folds=folds,
                           primary_blocks=primary, source_paths=[path, labels])
