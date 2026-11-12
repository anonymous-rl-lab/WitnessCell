"""Fold-local H1 measurement calibration and original shrunk M0.

The only expression asset contains 10 singles and 21 allowed training pairs.
The four development-query expression vectors are absent. Every call selects
its allowed batch rows BEFORE numerical operations. Evaluation truth is not
returned by get_fold; fixed_reference is a separate training-pair target API.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
from pathlib import Path
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parent
HISTORY = ROOT.parent / 'h1_probe_20261026'
DEFAULT_ASSET = HISTORY / 'training_inputs/WitnessCell_Neural_Partner_Pilot_20261021/joung_assets/ready/TRAINING_BATCH_ONLY.npz'
DEFAULT_KNOWLEDGE = HISTORY / 'training_inputs/WitnessCell_Neural_Partner_Pilot_20261021/assets/GenePT_fixed64.npz'
PARAMS = dict(alpha=0.5, noise=0.3, gamma=0.14684699381999897)
OUTER_PAIRS = ('CDX1+HNF4A', 'CDX1+MSGN1', 'HNF4A+NHLH1', 'MSGN1+NHLH1')


def incidence(edges, n):
    edges = np.asarray(edges, dtype=int).reshape(-1, 2)
    out = np.zeros((len(edges), n), dtype=float)
    out[np.arange(len(edges))[:, None], edges] = 1.
    return out


def projection(edges, n):
    x = incidence(edges, n)
    p = np.eye(len(x)) - x @ np.linalg.pinv(x)
    return (p + p.T) / 2


class FoldData:
    """Minimal interface compatible with archived shared12 feature routines."""
    def __init__(self, d):
        self.d = d
        self.n, self.m, self.G = len(d['nodes']), len(d['train_edges']), len(d['genes'])
        self.X = incidence(d['train_edges'], self.n)
        self.Pinv = np.linalg.pinv(self.X)
        self.P = projection(d['train_edges'], self.n)
        self.R = d['train_residual']
        self.prank = int(round(np.trace(self.P)))
        self.baseline_max_abs = 0.
        self.feature_cache = {}

    def node_features(self, dims, source):
        key = (dims, source)
        if key in self.feature_cache:
            return self.feature_cache[key]
        s = self.d['single_centered_exact']
        a = self.d['knowledge_node'].astype(float)
        a -= a.mean(0)
        gs, gk = s @ s.T, a @ a.T
        gs /= np.trace(gs)
        gk /= np.trace(gk)
        gram = {'fused': .5 * (gs + gk), 'singles': gs, 'knowledge': gk}[source]
        ev, u = np.linalg.eigh(gram)
        ix = np.argsort(ev)[::-1]
        ev, u = ev[ix], u[:, ix]
        x = u[:, :dims] * np.sqrt(np.maximum(ev[:dims], 0))
        signs = np.sign(x[np.argmax(abs(x), axis=0), np.arange(dims)])
        signs[signs == 0] = 1
        x *= signs
        scale = np.sqrt(np.mean(x * x))
        x /= scale
        info = dict(dimension=dims, source=source,
                    feature_gram_energy_retained=float(ev[:dims].sum() / ev.sum()),
                    single_block_trace_normalized=True, knowledge_block_trace_normalized=True,
                    common_RMS_scale_before_normalization=float(scale), per_axis_whitening=False)
        self.feature_cache[key] = x, info
        return x, info


class FoldFactory:
    def __init__(self, asset_path=None, knowledge_path=None, asset=None, knowledge_node=None):
        self.asset_path = Path(asset_path or (ROOT / 'data/TRAINING_BATCH_ONLY.npz'
                               if (ROOT / 'data/TRAINING_BATCH_ONLY.npz').exists() else DEFAULT_ASSET))
        if asset is None:
            with np.load(self.asset_path, allow_pickle=False) as z:
                asset = {k: z[k].copy() for k in z.files}
        self.asset = asset
        self.nodes = [str(x) for x in asset['nodes']]
        self.genes = asset['genes'].copy()
        self.pairs = [str(x) for x in asset['train_pairs']]
        assert len(self.nodes) == 10 and len(self.pairs) == 21
        assert not set(OUTER_PAIRS).intersection(asset['condition_ids'].astype(str))
        self.edge_array = np.array([[self.nodes.index(g) for g in p.split('+')] for p in self.pairs])
        self.pair_index = {tuple(sorted(e)): i for i, e in enumerate(self.edge_array)}
        if knowledge_node is None:
            local = ROOT / 'data/KNOWLEDGE_NODE.npz'
            if local.exists() and knowledge_path is None:
                with np.load(local, allow_pickle=False) as z:
                    assert np.array_equal(z['nodes'], asset['nodes'])
                    knowledge_node = z['knowledge_node'].copy()
            else:
                path = Path(knowledge_path or DEFAULT_KNOWLEDGE)
                with np.load(path, allow_pickle=False) as z:
                    lut = {str(g): i for i, g in enumerate(z['genes'])}
                    knowledge_node = z['embedding'][[lut[g] for g in self.nodes]].copy()
        self.knowledge_node = knowledge_node
        self._cache = {}
        self._raw_halves = None

    def edges(self, values):
        """Resolve IDs, canonical strings or node-index pairs."""
        a = np.asarray(values)
        if a.size == 0:
            return np.empty((0, 2), int)
        if a.ndim == 2 and a.shape[1] == 2:
            out = a.astype(int)
        elif a.dtype.kind in 'iu':
            out = self.edge_array[a.astype(int)]
        else:
            out = np.array([[self.nodes.index(g) for g in str(p).split('+')] for p in a], int)
        out = np.sort(out, axis=1)
        assert np.all((out >= 0) & (out < len(self.nodes))) and np.all(out[:, 0] != out[:, 1])
        return out

    def pair_ids(self, values):
        return np.array([self.pair_index[tuple(e)] for e in self.edges(values)], int)

    def get_fold(self, train_edges, eval_edges=(), use_cache=True):
        """Fit calibration/M0 on support pairs and predict any distinct query.

        Held-out identities are allowed. Their expression and counts are never
        consumed by this computation. No evaluation response is returned.
        """
        train = self.pair_ids(train_edges)
        qe = self.edges(eval_edges)
        assert len(train) > 0 and len(set(train)) == len(train)
        te = self.edge_array[train]
        assert not set(map(tuple, te)).intersection(map(tuple, qe))
        key = (tuple(train), tuple(map(tuple, qe)))
        if use_cache and key in self._cache:
            return self._cache[key]
        a = self.asset
        n = len(self.nodes)
        allowed = np.r_[np.arange(n), n + train]
        obs_mask = np.isin(a['obs_condition'], allowed)
        oi = np.flatnonzero(obs_mask)
        qi, bi = a['obs_condition'][oi], a['obs_batch'][oi]
        values = a['batch_values'][oi]  # only allowed expression is accessed
        counts = a['counts'][allowed].astype(float)
        p = counts / counts.sum(1)[:, None]
        b = counts.shape[1]
        L = np.zeros((b, b))
        F = np.zeros((b, len(oi)))
        pool = np.zeros((len(allowed), len(oi)))
        for k, q in enumerate(allowed):
            pq = p[k]
            L += np.diag(pq) - np.outer(pq, pq)
            for o in np.flatnonzero(qi == q):
                F[:, o] = -pq[bi[o]] * pq
                F[bi[o], o] += pq[bi[o]]
                pool[k, o] = pq[bi[o]]
        assert np.linalg.matrix_rank(L, tol=1e-10) == b - 1
        A = np.linalg.pinv(L, rcond=1e-12) @ F
        A -= np.ones((b, 1)) @ (a['pi'] @ A)[None, :]
        C = pool - (p - a['pi']) @ A
        measured = C @ values
        singles, y = measured[:n], measured[n:]
        origin = singles.mean(0)
        centered = singles - origin
        xt, xq = incidence(te, n), incidence(qe, n)
        addt, addq = .5 * xt @ centered, .5 * xq @ centered
        st = origin + addt / (1 + PARAMS['alpha'] * abs(addt))
        sq = origin + addq / (1 + PARAMS['alpha'] * abs(addq))
        xc = xt - xt.mean(0)
        r = y - st
        u = r.mean(0)
        beta = np.linalg.solve(xc.T @ xc + 2 * PARAMS['noise'] * np.eye(n), xc.T @ (r - u))
        bt = st + u + PARAMS['gamma'] * xc @ beta
        bq = sq + u + PARAMS['gamma'] * (xq - xt.mean(0)) @ beta
        d = dict(nodes=np.array(self.nodes), genes=self.genes.copy(),
                 single_means=singles, single_centered_exact=centered,
                 knowledge_node=self.knowledge_node, train_ids=train.copy(),
                 train_edges=te.copy(), train_pairs=np.array(self.pairs)[train],
                 train_y=y, train_base=bt, train_residual=y-bt,
                 eval_edges=qe.copy(), query_edges=qe.copy(),
                 eval_pairs=np.array(['+'.join(self.nodes[i] for i in e) for e in qe]),
                 query_pairs=np.array(['+'.join(self.nodes[i] for i in e) for e in qe]),
                 eval_base=bq, query_base=bq, origin=origin,
                 support_counts=a['counts'][n+train].copy(),
                 single_counts=a['counts'][:n].copy(),
                 allowed_conditions=allowed, allowed_observation_rows=oi,
                 M0_parameters=dict(PARAMS), m0_alpha=np.array(PARAMS['alpha']),
                 m0_noise=np.array(PARAMS['noise']), m0_gamma=np.array(PARAMS['gamma']),
                 train_projection=projection(te,n), eval_projection=projection(qe,n))
        out = FoldData(d)
        assert np.isfinite(y).all() and np.isfinite(bt).all() and np.isfinite(bq).all()
        if use_cache:
            # Episode cross-fitting can visit many support graphs. Bound retained
            # full-RNA arrays rather than keeping an unbounded in-memory archive.
            while len(self._cache) >= 8:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = out
        return out

    def fixed_reference(self, target_edges):
        """Explicit scoring/episode-target API for allowed training pairs only.

        The fixed seven-library direct mean is not processed by the fitted
        calibration. No outer-query expression exists in this asset.
        """
        ids = self.pair_ids(target_edges)
        a = self.asset
        out = np.zeros((len(ids), len(self.genes)))
        for k, i in enumerate(ids):
            obs = np.flatnonzero(a['obs_condition'] == len(self.nodes) + i)
            w = a['pi'][a['obs_batch'][obs]]
            assert abs(w.sum() - 1) < 1e-12
            nz = obs[w > 0]
            out[k] = w[w > 0] @ a['batch_values'][nz]
        return out

    def manifest(self):
        return dict(training_asset=str(self.asset_path),
                    training_asset_sha256=hashlib.sha256(self.asset_path.read_bytes()).hexdigest(),
                    allowed_pairs=self.pairs, nodes=self.nodes, readouts=len(self.genes),
                    M0_parameters=dict(PARAMS),
                    query_expression_available=False,
                    technical_halves='Optional separate training-only raw-composition technical halves; calibrated/fixed7pi half arrays unavailable.')

    def reliability(self, train_edges):
        """Old-R-form shrinkage from actual support-side raw technical halves.

        These are raw cell-mixture halves, NOT calibrated/fixed7pi replicates.
        Their correspondence to the calibrated target is only an auxiliary
        control. Training pairs are selected before projection/weight fitting.
        """
        ids = self.pair_ids(train_edges)
        if self._raw_halves is None:
            path = ROOT / 'data/TRAINING_RAW_HALVES.npz'
            with np.load(path, allow_pickle=False) as z:
                assert np.array_equal(z['pairs'], self.pairs)
                assert np.array_equal(z['genes'], self.genes)
                self._raw_halves = z['halves'].copy()
        p = projection(self.edge_array[ids], len(self.nodes))
        a, b = [p @ h[ids] for h in self._raw_halves]
        denominator = np.sum(a*a + b*b, axis=0)
        numerator = 2*np.sum(a*b, axis=0)
        return np.divide(numerator, denominator, out=np.zeros_like(numerator),
                         where=denominator > 0).clip(0, 1)
