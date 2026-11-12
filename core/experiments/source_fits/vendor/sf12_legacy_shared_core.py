"""Real-RNA adapter: preserve supplied M0; refit protected shared coefficients.

This is algebraically equivalent to the archived SF12 at matching settings,
not a novel model or a dynamic-ATAC experiment.
"""
import numpy as np

RHO = .1
AMPLITUDE = .6717359152620326

def incidence(edges, n):
    x = np.zeros((len(edges), n))
    x[np.arange(len(edges))[:, None], np.asarray(edges, int)] = 1.
    return x

def features(singles, knowledge, train_edges, query_edges):
    scale = np.sqrt(np.mean(singles*singles, axis=0))
    floor = max(1e-10, float(np.median(scale[scale > 1e-12]))*1e-3)
    scale = np.maximum(scale, floor)
    t = np.tanh(singles/scale)
    a = knowledge.astype(float).copy()
    a -= a.mean(0)
    gram = a@a.T
    gram /= np.trace(gram)
    ev, u = np.linalg.eigh(gram)
    ix = np.argsort(ev)[::-1]
    x = u[:, ix[:2]]*np.sqrt(np.maximum(ev[ix[:2]], 0))
    signs = np.sign(x[np.argmax(abs(x), axis=0), np.arange(2)])
    signs[signs == 0] = 1
    x *= signs
    x /= np.sqrt(np.mean(x*x))
    prod = x[train_edges[:, 0]]*x[train_edges[:, 1]]
    center = prod.mean(0)
    ps = np.maximum(np.sqrt(np.mean((prod-center)**2, axis=0)), 1e-8)
    def make(edges):
        a, b = t[edges[:, 0]], t[edges[:, 1]]
        ap, an = np.maximum(a, 0), np.maximum(-a, 0)
        bp, bn = np.maximum(b, 0), np.maximum(-b, 0)
        op = ap*bn+an*bp
        loc = np.stack([ap*bp, an*bn, op, op*np.tanh(a+b)], axis=-1)*scale[None, :, None]
        mod = np.column_stack([np.ones(len(edges)), (x[edges[:, 0]]*x[edges[:, 1]]-center)/ps])
        return (loc[:, :, :, None]*mod[:, None, None, :]).reshape(len(edges), len(scale), 12)
    return make(train_edges), make(query_edges)

def fit_fixed_m0(train_residual, train_x, train_features, rho=RHO, amplitude=AMPLITUDE):
    """Residual is measured Y minus the unchanged real M0, including saturation."""
    pinv = np.linalg.pinv(train_x)
    p = np.eye(len(train_x))-train_x@pinv
    p = (p+p.T)/2
    d = np.einsum('ab,bgp->agp', p, train_features, optimize=True)
    m = np.einsum('egp,egq->pq', d, d, optimize=True)
    b = np.einsum('egp,eg->p', d, p@train_residual, optimize=True)
    m = (m+m.T)/2
    sc = np.sqrt(np.maximum(np.diag(m), 0))
    sc = np.maximum(sc, max(float(sc.max(initial=0))*1e-10, 1e-15))
    normal = m/sc[:, None]/sc[None, :]
    theta = np.linalg.solve(normal+rho*np.eye(len(b)), b/sc)/sc
    fitted = np.einsum('egp,p->eg', train_features, theta*amplitude)
    return dict(theta=theta, shift=pinv@fitted, amplitude=float(amplitude),
                train_delta=p@fitted, gram=m, rhs=b,
                independent_edge_contrasts=len(train_x)-int(np.linalg.matrix_rank(train_x)),
                projected_feature_rank=int(np.linalg.matrix_rank(normal)))

def predict(state, frozen_m0, query_x, query_features):
    delta = np.einsum('egp,p->eg', query_features, state['theta']*state['amplitude'])-query_x@state['shift']
    return frozen_m0+delta, delta

def score(truth, base, prediction, contrast):
    r, e = truth-base, truth-prediction
    pr, pe = contrast@r, contrast@e
    delta = contrast@(prediction-base)
    denom = np.linalg.norm(pr)*np.linalg.norm(delta)
    return dict(full_SSE=float(np.sum(e*e)), partner_SSE=float(pe@pe),
                full_gain_pct=float(100*(1-np.sum(e*e)/np.sum(r*r))),
                partner_gain_pct=float(100*(1-(pe@pe)/(pr@pr))),
                partner_direction_cosine=float(pr@delta/denom) if denom else None,
                per_query_full_gain_pct=(100*(1-np.sum(e*e, axis=1)/np.sum(r*r, axis=1))).tolist())
