"""Exact specialization of the three features used by locked inc_amp_fp_r02.

No gate, fitted coefficient, single-perturbation input or output coordinate is
removed. The original broad probe bank calculates unused candidate features.
This adapter is intentionally valid only for the frozen v14 inc_amp_fp recipe.
"""
from types import SimpleNamespace
import numpy as np
from identity_head import self_statistic

def build_bank(target,effects,genes,gene2go,go_top_k,support_sizes,neighbor_sizes):
    gene_index={str(gene):index for index,gene in enumerate(genes)}
    matrix=np.stack(list(effects.values())).astype(np.float64)
    target_index=gene_index.get(target)
    self_anchor=np.zeros(len(genes),dtype=np.float64)
    fingerprint_tail_anchor=np.zeros(len(genes),dtype=np.float64)
    fingerprint_rms_anchor=np.zeros(len(genes),dtype=np.float64)
    if target_index is not None:
        self_anchor[target_index]=self_statistic(effects,gene_index,'mean')
        fingerprint=matrix[:,target_index]
        tail_quantile=.10 if self_anchor[target_index]<0. else .90
        fingerprint_tail_anchor[target_index]=float(np.quantile(fingerprint,tail_quantile))
        fingerprint_rms_anchor[target_index]=float(np.sign(self_anchor[target_index])*np.sqrt(np.mean(np.square(fingerprint))))
    return SimpleNamespace(self_anchor=self_anchor,fingerprint_tail_anchor=fingerprint_tail_anchor,fingerprint_rms_anchor=fingerprint_rms_anchor)

def install():
    import multihead_v14
    assert multihead_v14.RECIPE.family=='inc_amp_fp'
    assert multihead_v14.RECIPE.name=='inc_amp_fp_r02'
    multihead_v14.build_bank=build_bank
