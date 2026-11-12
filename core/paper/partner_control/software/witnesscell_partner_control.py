#!/usr/bin/env python3
"""Train-side amplitude fitting and outcome-free, fixed-panel partner control.

This module fits only a scalar using supplied out-of-fold predictions. It does
not fit the underlying response model. Losses are unweighted squared errors.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path
import numpy as np


def incidence(edges, nodes):
    edges = np.asarray(edges)
    if edges.ndim != 2 or edges.shape[1] != 2 or len(edges) == 0:
        raise ValueError('edges must be a nonempty N by 2 array')
    if not np.issubdtype(edges.dtype, np.integer):
        raise ValueError('edges must contain integer endpoint indices')
    if nodes < 2 or np.any(edges < 0) or np.any(edges >= nodes):
        raise ValueError('endpoint index out of range')
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError('two distinct endpoints are required')
    if len(np.unique(np.sort(edges, axis=1), axis=0)) != len(edges):
        raise ValueError('duplicate or reversed duplicate unordered pairs')
    x = np.zeros((len(edges), nodes), dtype=float)
    x[np.arange(len(edges))[:, None], edges] = 1.
    return x


def array2d(value, rows, name):
    a = np.asarray(value, dtype=float)
    if a.ndim == 1:
        a = a[:, None]
    if a.ndim != 2 or a.shape[0] != rows or a.shape[1] == 0 or not np.isfinite(a).all():
        raise ValueError(name + ' must be a finite N by R array')
    return a


def project(edges, nodes, values):
    x = incidence(edges, nodes)
    v = array2d(values, len(edges), 'values')
    endpoint, _, rank, _ = np.linalg.lstsq(x, v, rcond=None)
    endpoint = x @ endpoint
    partner = v - endpoint
    dim = len(edges) - int(rank)
    if dim == 0:
        endpoint = v.copy()
        partner = np.zeros_like(v)
    return endpoint, partner, dim


def fit_amplitude(y, baseline, prediction, edges, folds, nodes):
    """Fit alpha from archived held-combination predictions, fold by fold.

    The caller must generate baseline and prediction without the held fold's
    outcomes. This function cannot certify that provenance from arrays alone.
    """
    incidence(edges, nodes)
    y = array2d(y, len(edges), 'y')
    baseline = array2d(baseline, len(edges), 'baseline')
    prediction = array2d(prediction, len(edges), 'prediction')
    if y.shape != baseline.shape or y.shape != prediction.shape:
        raise ValueError('response arrays must have identical shapes')
    folds = np.asarray(folds)
    if folds.shape != (len(edges),) or len(np.unique(folds)) < 2:
        raise ValueError('at least two recorded inner folds are required')
    cross = energy = residual_energy = 0.
    details = []
    for f in np.unique(folds):
        mask = folds == f
        _, r, dim = project(np.asarray(edges)[mask], nodes, y[mask] - baseline[mask])
        _, d, _ = project(np.asarray(edges)[mask], nodes, prediction[mask] - baseline[mask])
        aa = float(np.sum(r*d)); bb = float(np.sum(d*d)); ee = float(np.sum(r*r))
        cross += aa; energy += bb; residual_energy += ee
        details.append({'fold':str(f),'pairs':int(mask.sum()),'contrast_dim':dim,'cross':aa,'correction_energy':bb,'residual_energy':ee})
    alpha = float(np.clip(cross/energy, 0, 1)) if energy > 1e-14 else 0.
    rho = cross / np.sqrt(energy*residual_energy) if energy > 1e-14 and residual_energy > 1e-14 else None
    return dict(alpha=alpha, cross=cross, correction_energy=energy,
                residual_energy=residual_energy, rho=rho,
                inner_gain_pp=100*(2*cross-energy)/residual_energy if residual_energy > 1e-14 else None,
                inner_adjusted_gain_pp=100*(2*alpha*cross-alpha*alpha*energy)/residual_energy if residual_energy > 1e-14 else None,
                folds=details, rule='bounded least squares on foldwise projected internal errors',
                status='calibration fit, not a forecast guarantee', query_outcomes_used=False)


def predict_panel(baseline, prediction, edges, nodes, alpha):
    """Use identities and predictions only; no query response argument exists."""
    baseline = array2d(baseline, len(edges), 'baseline')
    prediction = array2d(prediction, len(edges), 'prediction')
    if baseline.shape != prediction.shape:
        raise ValueError('prediction and baseline shapes differ')
    if not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError('alpha must lie in [0,1]')
    de, dp, dim = project(edges, nodes, prediction-baseline)
    return dict(baseline=baseline, original=prediction,
                whole_shrink=baseline+alpha*(prediction-baseline),
                endpoint_correction=de, partner_correction=dp,
                partner_only=baseline+de+alpha*dp, contrast_dim=np.array(dim))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fit'); f.add_argument('--input', required=True); f.add_argument('--output', required=True)
    q = sub.add_parser('predict'); q.add_argument('--input', required=True); q.add_argument('--calibration', required=True); q.add_argument('--output', required=True)
    args = p.parse_args()
    with np.load(args.input, allow_pickle=False) as a:
        if args.command == 'fit':
            required = {'y','baseline','prediction','edges','folds','nodes'}
            if set(a.files) != required: raise ValueError('fit input fields must be '+str(sorted(required)))
            result=fit_amplitude(a['y'],a['baseline'],a['prediction'],a['edges'],a['folds'],int(a['nodes']))
            result['input_sha256']=hashlib.sha256(Path(args.input).read_bytes()).hexdigest()
            Path(args.output).parent.mkdir(parents=True,exist_ok=True)
            Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
        else:
            required={'baseline','prediction','edges','nodes'}
            if set(a.files) != required: raise ValueError('predict accepts identities and predictions only; fields must be '+str(sorted(required)))
            alpha=float(json.loads(Path(args.calibration).read_text())['alpha'])
            result=predict_panel(a['baseline'],a['prediction'],a['edges'],int(a['nodes']),alpha)
            Path(args.output).parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(args.output,edges=a['edges'],**result)
    print(json.dumps({'command':args.command,'output':args.output,'status':'complete'}))

if __name__ == '__main__': main()
