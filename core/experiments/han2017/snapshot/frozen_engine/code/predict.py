"""Frozen inference and source-only coefficient arithmetic. No optimizer or target outcomes."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import argparse, copy, datetime, hashlib, json, sys
import numpy as np
import pandas as pd
import torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'vendor'))
from frozen_models import FittedModel
from partner_model_v4 import PartnerConfig, PartnerHeadV4
import sf12_legacy_shared_core as legacy
torch.set_num_threads(1)
torch.set_grad_enabled(False)
CELLS = ['HeLa', 'A549', '293T']
FAMILIES = ['rank1', 'spectral', 'sf12', 'v4']

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def arrays(p):
    with np.load(p, allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}

def project(edges, n, value):
    x = legacy.incidence(edges, n)
    return value-x@np.linalg.lstsq(x, value, rcond=None)[0]

def coefficient(blocks, n):
    A = B = Af = Bf = E = 0.
    for edges, error, delta in blocks:
        r = project(edges, n, error)
        d = project(edges, n, delta)
        A += float(np.sum(r*d)); B += float(np.sum(d*d))
        Af += float(np.sum(error*delta)); Bf += float(np.sum(delta*delta))
        E += float(np.sum(r*r))
    return dict(alpha=float(np.clip(A/B, 0, 1)) if B > 1e-14 else 0.,
                gamma=float(np.clip(Af/Bf, 0, 1)) if Bf > 1e-14 else 0.,
                partner_cross=A, partner_energy=B, partner_baseline_sse=E)

def sf_features(s, fixed_scale, embedding_coordinates, center, spread, edges):
    t = np.tanh(s/fixed_scale)
    a, b = t[edges[:, 0]], t[edges[:, 1]]
    ap, an = np.maximum(a, 0), np.maximum(-a, 0)
    bp, bn = np.maximum(b, 0), np.maximum(-b, 0)
    op = ap*bn+an*bp
    loc = np.column_stack([ap*bp, an*bn, op, op*np.tanh(a+b)])*fixed_scale
    mod = np.column_stack([np.ones(len(edges)),
        (embedding_coordinates[edges[:, 0]]*embedding_coordinates[edges[:, 1]]-center)/spread])
    return (loc[:, :, None]*mod[:, None, :]).reshape(len(edges), 12)

class Source:
    def __init__(self, cell):
        self.cell = cell
        self.path = p = ROOT/'pretrained'/cell
        inp = json.loads((p/'source_input.json').read_text())
        self.genes = inp['genes']; self.n = len(self.genes)
        self.lookup = {g: i for i, g in enumerate(self.genes)}
        self.singles = np.array(inp['singles'])
        self.smean = self.singles.mean(1)
        self.scales = np.std(self.singles, axis=0)
        self.h = 1./float(np.mean(1/self.scales))
        self.rms = float(np.sqrt(np.mean(self.smean**2)))
        self.train_edges = np.array([[self.lookup[g] for g in v['pair']] for v in inp['training']])
        self.query = np.array([[self.lookup[g] for g in pair] for pair in inp['queries']])
        self.models = [FittedModel.from_dict(v) for v in json.loads((p/'M0_models.json').read_text())]
        self.native_archive = arrays(p/'predictions.npz')
        self.phase_archive = arrays(p/'phase2_predictions.npz')
        self.sel = json.loads((p/'selection.json').read_text())
        self.phase_sel = json.loads((p/'phase2_selection.json').read_text())['selected']
        self.sf_state = arrays(p/'sf12_parameters.npz')
        knowledge = arrays(p/'knowledge.npz')
        assert np.array_equal(knowledge['genes'], np.array(self.genes))
        self.knowledge = knowledge['embedding']
        # Reconstruct the original source feature transform; do not refit using Han.
        a = self.knowledge.astype(float).copy(); a -= a.mean(0)
        gram = a@a.T; gram /= np.trace(gram)
        ev, u = np.linalg.eigh(gram); ix = np.argsort(ev)[::-1]
        x = u[:, ix[:2]]*np.sqrt(np.maximum(ev[ix[:2]], 0))
        signs = np.sign(x[np.argmax(abs(x), axis=0), np.arange(2)])
        signs[signs == 0] = 1; x *= signs; x /= np.sqrt(np.mean(x*x))
        self.x = x
        prod = x[self.train_edges[:, 0]]*x[self.train_edges[:, 1]]
        self.center = prod.mean(0)
        self.spread = np.maximum(np.sqrt(np.mean((prod-self.center)**2, axis=0)), 1e-8)
        kk = self.knowledge.astype(np.float32); kk -= kk.mean(0)
        kk /= np.maximum(np.sqrt(np.mean(kk*kk, axis=1, keepdims=True)), 1e-8)
        self.kt = torch.tensor(kk)
        self.networks = []
        ns = self.phase_sel['network_selected']
        if ns['key'] != 'M0':
            for seed in (17, 29):
                state = torch.load(p/f"{ns['key']}_seed{seed}.pt", map_location='cpu', weights_only=False)
                model = PartnerHeadV4(PartnerConfig(**state['config']))
                model.load_state_dict(state['state_dict']); model.eval()
                self.networks.append((model, state))
        self.meta = {}
        for family in FAMILIES:
            if family in ('rank1', 'spectral'):
                sel = self.sel[family]
                z = arrays(p/(sel['best_regularized']+'_inner_predictions.npz'))
                delta = z['pred']-z['baseline'] if sel['selected'] != 'G0' else np.zeros_like(z['pred'])
                blocks = [(z['edges'][z['fold'] == f], (z['y']-z['baseline'])[z['fold'] == f],
                           delta[z['fold'] == f]) for f in np.unique(z['fold'])]
                beta = 1. if sel['selected'] != 'G0' else 0.
                selected_key = sel['selected']
            else:
                sel = self.phase_sel['sf12_selected' if family == 'sf12' else 'network_selected']
                beta = sel['beta']; selected_key = sel['key']; blocks = []
                for f in range(3):
                    z = arrays(p/f'phase2_inner_{f}.npz')
                    delta = beta*z[sel['key']] if sel['key'] != 'M0' else np.zeros_like(z['M0'])
                    blocks.append((z['edges'], z['y']-z['M0'], delta))
            self.meta[family] = dict(selected_key=selected_key, native_beta=beta,
                                     **coefficient(blocks, self.n))

    def raw_delta(self, family, edges, singles_source_units):
        if family in ('rank1', 'spectral'):
            if self.sel[family]['selected'] == 'G0':
                return np.zeros(len(edges))
            return self.native_archive[family+'_G'][edges[:, 0], edges[:, 1]]
        sel = self.phase_sel['sf12_selected' if family == 'sf12' else 'network_selected']
        if sel['key'] == 'M0':
            return np.zeros(len(edges))
        if family == 'sf12':
            F = sf_features(singles_source_units, self.rms, self.x, self.center, self.spread, edges)
            shift = self.sf_state[sel['key']+'__shift']
            return F@self.sf_state[sel['key']+'__theta']-shift[edges[:, 0]]-shift[edges[:, 1]]
        out = np.zeros(len(edges))
        st = torch.tensor(singles_source_units/self.rms, dtype=torch.float32)
        a, b = torch.tensor(edges.T.copy(), dtype=torch.long)
        for model, state in self.networks:
            raw = model(self.kt[a], self.kt[b], torch.zeros((1, 64)), st[a], st[b], torch.zeros(())).numpy()*state['scale']
            shift = state['endpoint_shift']
            out += (raw-shift[edges[:, 0]]-shift[edges[:, 1]])/len(self.networks)
        return out

    def replay(self):
        checks = {}
        base = np.column_stack([m.predict(self.query)*s for m, s in zip(self.models, self.scales)])
        checks['M0'] = float(np.max(abs(base-self.native_archive['M0'])))
        checks['phase2_M0'] = float(np.max(abs(base-self.phase_archive['M0'])))
        for family in FAMILIES:
            raw = self.raw_delta(family, self.query, self.smean)
            if family in ('rank1', 'spectral'):
                ref = self.native_archive[family+'_selected']-base
            else:
                sel = self.phase_sel['sf12_selected' if family == 'sf12' else 'network_selected']
                ref = self.phase_archive[family+'_selected' if family == 'sf12' else 'network_selected']-base
            pred = self.meta[family]['native_beta']*raw[:, None]
            checks[family] = float(np.max(abs(pred-ref)))
        # Exact original 12-feature implementation comparison, without fitting.
        _, original = legacy.features(self.smean[:, None], self.knowledge, self.train_edges, self.query)
        ours = sf_features(self.smean, self.rms, self.x, self.center, self.spread, self.query)
        checks['sf12_feature_map'] = float(np.max(abs(original[:, 0]-ours)))
        assert max(checks.values()) < 1e-7, checks
        return checks

    def predict(self, genes, pairs, singles):
        local = {g: i for i, g in enumerate(genes)}
        edges = np.array([[self.lookup[g] for g in p] for p in pairs])
        le = np.array([[local[g] for g in p] for p in pairs])
        scale = float(np.std(singles)); assert scale > 0
        blocks = []
        for model in self.models:
            m = copy.deepcopy(model)
            for g in set(np.array(pairs).ravel()):
                m.a[self.lookup[g]] = singles[local[g]]/scale
            blocks.append(m.predict(edges))
        base = scale*np.mean(blocks, axis=0)
        sf_single = self.smean.copy()
        for g in set(np.array(pairs).ravel()):
            sf_single[self.lookup[g]] = singles[local[g]]*self.h/scale
        output = {'M0': base,
                  'single_additive': np.array([singles[local[a]]+singles[local[b]] for a, b in pairs])}
        for family in FAMILIES:
            raw = self.raw_delta(family, edges, sf_single)*scale/self.h
            meta = self.meta[family]; d = meta['native_beta']*raw
            pd = project(le, len(genes), d); ep = d-pd
            variants = {'raw_selected_direction':base+raw, 'native_selected':base+d,
                        'whole_shrink':base+meta['alpha']*d, 'partner_off':base+ep,
                        'partner_control':base+ep+meta['alpha']*pd,
                        'full_calibrated_whole':base+meta['gamma']*d}
            # Symmetric pair inputs are required for every model.
            reverse = self.raw_delta(family, edges[:, ::-1], sf_single)*scale/self.h
            assert np.max(abs(reverse-raw)) < 1e-5*max(1., np.max(abs(raw)))
            for arm, v in variants.items():
                assert np.isfinite(v).all(); output[family+'__'+arm] = v
        return output

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', default='results')
    args = ap.parse_args(); out = ROOT/args.output; out.mkdir(exist_ok=True, parents=True)
    if (out/'PREDICTIONS_LOCKED.json').exists():
        raise RuntimeError('Preserve frozen output; choose a new directory')
    files = sorted(p for p in (ROOT/'pretrained').rglob('*') if p.is_file())
    before = {str(p.relative_to(ROOT)): sha(p) for p in files}
    sources = {c: Source(c) for c in CELLS}
    checks = {c: s.replay() for c, s in sources.items()}
    old_dec = pd.read_csv(ROOT/'source_replay/v33_decisions.csv')
    for c, s in sources.items():
        ref = float(old_dec.loc[old_dec.feature_id == c+'_high', 'alpha'].iloc[0])
        assert abs(s.meta['rank1']['alpha']-ref) < 1e-9
    metadata = {c: s.meta for c, s in sources.items()}
    for readout in ('raw_repeats', 'published'):
        inp = json.loads((ROOT/'inputs'/f'{readout}.json').read_text())
        singles = np.array(inp['singles']); pairs = inp['pairs']; genes = inp['genes']
        payload = {'pairs': np.array(pairs), 'genes': np.array(genes), 'single_input':singles}
        rows = []
        for c, s in sources.items():
            byrep = [s.predict(genes, pairs, singles[:, r]) for r in range(singles.shape[1])]
            for arm in byrep[0]:
                v = np.column_stack([a[arm] for a in byrep]); payload[c+'__'+arm] = v
                for i, pair in enumerate(pairs):
                    for r in range(singles.shape[1]):
                        rows.append(dict(source_model=c, arm=arm, repeat=r+1, pair='|'.join(pair), prediction=v[i,r]))
        if readout == 'raw_repeats':
            old = arrays(ROOT/'source_replay/previous_han_predictions.npz')
            assert np.array_equal(old['pairs'], payload['pairs'])
            for c in CELLS:
                for now, prior in [('M0','baseline'), ('rank1__native_selected','original'),
                                    ('rank1__partner_control','partner_only'), ('single_additive','single_additive')]:
                    error = float(np.max(abs(payload[c+'__'+now]-old[c+'__'+prior])))
                    checks[c]['previous_han_'+now] = error
                    assert error < 1e-8, (c, now, error)
        np.savez_compressed(out/f'{readout}_predictions.npz', **payload)
        pd.DataFrame(rows).to_csv(out/f'{readout}_predictions.csv', index=False)
    assert before == {str(p.relative_to(ROOT)):sha(p) for p in files}
    (out/'SOURCE_REPLAY_CHECKS.json').write_text(json.dumps(checks, indent=2))
    (out/'SOURCE_DECISIONS.json').write_text(json.dumps(metadata, indent=2))
    lock = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                purpose='Exploratory comparison; earlier Han results already known.',
                model_fit_calls=0, target_calibration_calls=0, target_outcomes_read_by_prediction=False,
                source_weights_unchanged=True, model_sha256=before,
                code_sha256=sha(__file__), protocol_sha256=sha(ROOT/'PROTOCOL.json'),
                outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (out/'PREDICTIONS_LOCKED.json').write_text(json.dumps(lock, indent=2))
    print(json.dumps(dict(stage='frozen_inference_complete', checks=checks, source_decisions=metadata), indent=2))

if __name__ == '__main__':
    main()
