"""Shared residual matrix models; unique undirected observations, free diagonal.

The fitted objective is mean((rbar-G[edges])**2) + tau/n*||G||_* +
nu/n**2*||G||_F**2.  If r contains two repetitions, their mean is used.
The omitted mean-squared difference between repetitions is constant in G.
All deadlines are absolute time.monotonic() values.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.optimize import minimize


class _Deadline(Exception):
    pass


def _inputs(n, edges, r):
    edges = np.asarray(edges, dtype=np.int64)
    r = np.asarray(r, dtype=float)
    if r.ndim == 2:
        r = r.mean(axis=1)
    if edges.ndim != 2 or edges.shape[1] != 2 or len(edges) == 0:
        raise ValueError("edges must be a nonempty m-by-2 array")
    if r.shape != (len(edges),) or not np.all(np.isfinite(r)):
        raise ValueError("responses must be finite and match edges")
    if np.any(edges < 0) or np.any(edges >= n):
        raise ValueError("edge endpoint outside matrix")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("diagonal measurements are not part of this protocol")
    if len(np.unique(np.sort(edges, axis=1), axis=0)) != len(edges):
        raise ValueError("each unordered pair must occur once")
    return edges, r


def _parts(G, edges, r, tau, nu):
    n = len(G)
    data = float(np.mean((G[edges[:, 0], edges[:, 1]]-r)**2))
    nuclear = float(tau/n*np.sum(np.abs(np.linalg.eigvalsh(G))))
    frobenius = float(nu/n**2*np.sum(G*G))
    return dict(objective=data+nuclear+frobenius, data_loss=data,
                penalty_nuclear=nuclear, penalty_frobenius=frobenius)


def data_loss_objective(G, edges, r, tau=0., nu=0.):
    G = np.asarray(G, dtype=float)
    edges, r = _inputs(len(G), edges, r)
    if G.shape != (len(G), len(G)) or not np.allclose(G, G.T, atol=1e-12, rtol=1e-12):
        raise ValueError("G must be square and symmetric")
    return _parts(G, edges, r, tau, nu)["objective"]


def _gradient(G, edges, r, nu):
    n = len(G)
    grad = (2*nu/n**2)*G
    err = (G[edges[:, 0], edges[:, 1]]-r)/len(edges)
    grad[edges[:, 0], edges[:, 1]] += err
    grad[edges[:, 1], edges[:, 0]] += err
    return grad


def smooth_gradient(G, edges, r, nu=0.):
    G = np.asarray(G, dtype=float)
    edges, r = _inputs(len(G), edges, r)
    return _gradient(G, edges, r, nu)


def tau_zero(n, edges, r):
    edges, r = _inputs(n, edges, r)
    grad = _gradient(np.zeros((n, n)), edges, r, 0.)
    return float(n*np.max(np.abs(np.linalg.eigvalsh(grad))))


def _prox_symmetric(Z, threshold):
    vals, vecs = np.linalg.eigh((Z+Z.T)/2)
    vals = np.sign(vals)*np.maximum(np.abs(vals)-threshold, 0.)
    out = (vecs*vals)@vecs.T
    return (out+out.T)/2


def fit_spectral(n, edges, r, tau, nu, max_iter=3000, tol=1e-6,
                 deadline=None, initial=None):
    """Monotone restarted FISTA with a full-matrix proximal stationarity check.

    stationarity_relative scales the norm of a valid objective subgradient
    by the loss gradient norm at zero.  With nu>0, objective_gap_bound is
    the strong-convexity upper bound ||subgradient||^2/(2*mu).
    """
    begun = time.monotonic()
    if tau < 0 or nu < 0:
        raise ValueError("penalties must be nonnegative")
    edges, r = _inputs(n, edges, r)
    X = np.zeros((n, n)) if initial is None else np.asarray(initial, float).copy()
    if X.shape != (n, n) or not np.all(np.isfinite(X)):
        raise ValueError("invalid initial matrix")
    X = (X+X.T)/2
    Y, momentum = X.copy(), 1.
    lipschitz = 1./len(edges)+2*nu/n**2
    step = 1./lipschitz
    reg = tau/n
    scale = max(float(np.linalg.norm(_gradient(np.zeros((n, n)), edges, r, 0.))), 1e-12)
    obj = _parts(X, edges, r, tau, nu)["objective"]
    converged, timed_out, restarts = False, False, 0
    stationarity, relative, gap_bound = np.inf, np.inf, np.inf
    iteration = 0
    for iteration in range(1, max_iter+1):
        if deadline is not None and time.monotonic() >= deadline:
            timed_out = True
            break
        grad_y = _gradient(Y, edges, r, nu)
        candidate = _prox_symmetric(Y-step*grad_y, step*reg)
        candidate_obj = _parts(candidate, edges, r, tau, nu)["objective"]
        # Acceleration occasionally overshoots; retain monotone descent.
        if candidate_obj > obj + 1e-14*max(1., abs(obj)):
            Y, momentum = X.copy(), 1.
            grad_y = _gradient(Y, edges, r, nu)
            candidate = _prox_symmetric(Y-step*grad_y, step*reg)
            candidate_obj = _parts(candidate, edges, r, tau, nu)["objective"]
            restarts += 1
        # A valid subgradient at candidate follows directly from prox optimality.
        subgradient = _gradient(candidate, edges, r, nu)-grad_y+(Y-candidate)/step
        stationarity = float(np.linalg.norm(subgradient))
        relative = stationarity/scale
        gap_bound = stationarity**2/(4*nu/n**2) if nu > 0 else np.inf
        previous = X
        X, obj = candidate, candidate_obj
        if relative <= tol:
            converged = True
            break
        next_momentum = (1+np.sqrt(1+4*momentum**2))/2
        Y = X+((momentum-1)/next_momentum)*(X-previous)
        momentum = next_momentum
    eigvals = np.linalg.eigvalsh(X)
    return dict(G=X, converged=bool(converged), timed_out=bool(timed_out),
                iterations=int(iteration), stationarity=stationarity,
                stationarity_relative=relative, objective_gap_bound=gap_bound,
                seconds=time.monotonic()-begun, restarts=restarts,
                numerical_rank=int(np.sum(np.abs(eigvals)>1e-8*max(1., np.max(np.abs(eigvals))))),
                eigenvalues=eigvals, **_parts(X, edges, r, tau, nu))


def fit_rank1(n, edges, r, tau, nu, max_iter=3000, tol=1e-6,
              deadline=None, initial=None):
    """Signed rank-one fit with exact soft-thresholded strength profiling.

    u has unit Euclidean norm without centering.  Six deterministic starts
    (four data eigenvectors and two random vectors) are used.  This is a
    nonconvex comparator; convergence denotes stationarity, not globality.
    """
    begun = time.monotonic()
    if tau < 0 or nu < 0:
        raise ValueError("penalties must be nonnegative")
    edges, r = _inputs(n, edges, r)
    i, j = edges.T
    m = len(edges)
    scale = max(float(np.mean(r*r)), 1e-12)
    M = -_gradient(np.zeros((n, n)), edges, r, 0.)
    vals, vecs = np.linalg.eigh(M)
    order = np.argsort(-np.abs(vals), kind="stable")[:min(4, n)]
    starts = [(f"spectral_{k}", vecs[:, ix]) for k, ix in enumerate(order)]
    rng = np.random.default_rng(20261018)
    starts += [(f"random_{k}", rng.normal(size=n)) for k in range(2)]
    if initial is not None:
        init = np.asarray(initial, float)
        if init.shape == (n, n):
            ev, evec = np.linalg.eigh((init+init.T)/2)
            init = evec[:, np.argmax(np.abs(ev))]
        if init.shape != (n,):
            raise ValueError("rank-one initial must be a vector or matrix")
        starts.insert(0, ("provided", init))

    def evaluate(w):
        norm = float(np.linalg.norm(w))
        if norm < 1e-15:
            raise ValueError("zero factor vector")
        u = w/norm
        v = u[i]*u[j]
        a = float(np.mean(v*v)+nu/n**2)
        b = float(np.mean(r*v))
        lam = float(np.sign(b)*max(abs(b)-tau/(2*n), 0.)/a) if a > 0 else 0.
        error = lam*v-r
        objective = float(np.mean(error*error)+tau/n*abs(lam)+nu/n**2*lam*lam)
        gu = np.zeros(n)
        coeff = 2*lam*error/m
        np.add.at(gu, i, coeff*u[j])
        np.add.at(gu, j, coeff*u[i])
        tangent = gu-u*np.dot(u, gu)
        return objective, tangent/norm, u, lam, float(np.linalg.norm(tangent))

    # Zero is an explicit admissible candidate, including when all starts
    # profile to zero.  It is globally optimal if tau >= tau_zero.
    zero = np.zeros((n, n))
    best = dict(G=zero, u=np.zeros(n), **{"lambda": 0.},
                objective=scale if np.any(r) else 0., stationarity=0.,
                stationarity_relative=0., converged=True, iterations=0,
                start="zero")
    best["objective"] = float(np.mean(r*r))
    records = []
    timed_out = False
    for name, start in starts:
        if deadline is not None and time.monotonic() >= deadline:
            timed_out = True
            break
        latest = [np.asarray(start, float).copy()]
        def fun(w):
            if deadline is not None and time.monotonic() >= deadline:
                raise _Deadline()
            latest[0] = w.copy()
            value, grad, _, _, _ = evaluate(w)
            # A fixed unit change makes L-BFGS-B stopping tolerances comparable
            # across response magnitudes.  Model selection and stationarity
            # are still evaluated with the unchanged original objective.
            return value/scale, grad/scale
        nit, success, message = 0, False, ""
        try:
            fitted = minimize(fun, start, jac=True, method="L-BFGS-B",
                              options=dict(maxiter=max_iter, ftol=1e-15,
                                           gtol=tol*.01, maxls=50))
            w, nit, success, message = fitted.x, int(fitted.nit), bool(fitted.success), str(fitted.message)
        except _Deadline:
            w, timed_out, message = latest[0], True, "deadline"
        objective, _, u, lam, tangent_norm = evaluate(w)
        relative = tangent_norm/scale
        stationary = bool(relative <= tol and np.isfinite(objective))
        record = dict(start=name, objective=objective, iterations=nit,
                      stationarity=tangent_norm, stationarity_relative=relative,
                      converged=stationary, scipy_success=success, message=message,
                      **{"lambda": lam})
        records.append(record)
        if objective < best["objective"]-1e-15:
            best = dict(G=lam*np.outer(u, u), u=u, **{"lambda": lam},
                        objective=objective, stationarity=tangent_norm,
                        stationarity_relative=relative, converged=stationary,
                        iterations=nit, start=name)
        if timed_out:
            break
    best.update(_parts(best["G"], edges, r, tau, nu))
    best.update(seconds=time.monotonic()-begun, timed_out=bool(timed_out),
                starts=records, starts_completed=len(records),
                numerical_rank=int(abs(best["lambda"])>1e-12))
    return best
