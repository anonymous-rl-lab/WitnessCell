"""Registered CPU-only Horlbeck prediction models (protocol Appendix B).

All inputs and predictions use the caller's scaled gamma units.  This module
never reads files, chooses data partitions, inspects test responses, or tunes
hyperparameters.  A deadline is an absolute ``time.monotonic()`` value.
"""
from dataclasses import dataclass, field
import time
import numpy as np
from scipy import linalg
from scipy.optimize import minimize

ALPHA_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)
ETA_GRID = (0.01, 0.1, 1.0)
RANKS = (0, 1, 2, 3)
START_SEEDS = (101, 202, 303)


class FitDeadlineExceeded(RuntimeError):
    """The registered global fitting allowance has expired."""


class FitConvergenceError(RuntimeError):
    """All registered starts failed. Caller applies the protocol fallback."""
    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


def _deadline_check(deadline):
    if deadline is not None and time.monotonic() >= deadline:
        raise FitDeadlineExceeded("Cumulative fitting deadline exceeded")


def _diagnostic_float(value):
    """Keep failed-start diagnostics valid JSON (null instead of NaN/Infinity)."""
    return float(value) if np.isfinite(value) else None


def _inputs(a, edges, y=None):
    a = np.asarray(a, dtype=np.float64)
    edges = np.asarray(edges)
    if a.ndim != 1 or not len(a) or not np.isfinite(a).all():
        raise ValueError("a must be a nonempty finite one-dimensional array")
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edges must have shape (m, 2)")
    if not np.issubdtype(edges.dtype, np.integer):
        raise ValueError("edge indices must be integers")
    edges = edges.astype(np.int64, copy=False)
    if edges.size and (edges.min() < 0 or edges.max() >= len(a)):
        raise ValueError("edge index outside a")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("Self pairs are excluded by this protocol")
    if y is not None:
        y = np.asarray(y, dtype=np.float64)
        if not len(edges) or y.shape != (len(edges),) or not np.isfinite(y).all():
            raise ValueError("Training y must be finite and match nonempty edges")
    return a, edges, y


def _mean_zero_basis(n):
    """Deterministic orthonormal Helmert columns spanning sum(v)=0."""
    result = np.zeros((n, max(0, n - 1)), dtype=np.float64)
    for j in range(n - 1):
        denom = np.sqrt((j + 1) * (j + 2))
        result[:j + 1, j] = 1.0 / denom
        result[j + 1, j] = -(j + 1) / denom
    return result


@dataclass
class FittedModel:
    a: np.ndarray
    supported: np.ndarray
    family: str
    rank: int
    alpha: float
    eta: float
    e0: float = 0.0
    beta: float = 1.0
    c: float = 0.0
    e: np.ndarray = field(default_factory=lambda: np.empty(0))
    d_beta: np.ndarray = field(default_factory=lambda: np.empty(0))
    d_c: np.ndarray = field(default_factory=lambda: np.empty(0))
    u: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    lambdas: np.ndarray = field(default_factory=lambda: np.empty(0))
    shared_theta: np.ndarray = field(default_factory=lambda: np.empty(0))
    diagnostics: dict = field(default_factory=dict)

    @property
    def converged(self):
        return bool(self.diagnostics.get("converged", False))

    def predict(self, edges):
        _, edges, _ = _inputs(self.a, edges)
        g, h = edges.T
        ag, ah = self.a[g], self.a[h]
        if self.family == "additive":
            return ag + ah
        if self.family == "endpoint":
            return ag + ah + self.e0 + self.e[g] + self.e[h]
        if self.family == "shared_quadratic":
            s = ag + ah
            t = self.shared_theta
            return t[0] + t[1] * s + t[2] * s * s + self.e[g] + self.e[h]
        if self.family != "horlbeck":
            raise ValueError("Unknown model family: " + self.family)
        result = (0.5 * (ag + ah + (self.beta + self.d_beta[g]) * ah
                         + (self.beta + self.d_beta[h]) * ag
                         + (self.c + self.d_c[g]) * ah * ah
                         + (self.c + self.d_c[h]) * ag * ag)
                  + self.e0 + self.e[g] + self.e[h])
        if self.rank:
            result = result + (self.u[g] * self.u[h]) @ self.lambdas
        return result

    def to_dict(self):
        result = {}
        for key, value in self.__dict__.items():
            result[key] = value.tolist() if isinstance(value, np.ndarray) else value
        return result

    @classmethod
    def from_dict(cls, values):
        result = dict(values)
        for key in ("a", "e", "d_beta", "d_c", "u", "lambdas", "shared_theta"):
            result[key] = np.asarray(result[key], dtype=np.float64)
        result["supported"] = np.asarray(result["supported"], dtype=np.int64)
        return cls(**result)


class _Design:
    """Linear constrained-baseline design plus latent-factor derivatives."""
    def __init__(self, a, edges, y, alpha, rank=0, eta=0.0):
        self.a, self.edges, self.y = _inputs(a, edges, y)
        self.alpha, self.eta, self.rank = float(alpha), float(eta), int(rank)
        self.supported = np.unique(edges)
        self.s = len(self.supported)
        self.basis = _mean_zero_basis(self.s)
        self.k = self.s - 1
        local = np.full(len(a), -1, dtype=int)
        local[self.supported] = np.arange(self.s)
        self.g, self.h = local[edges].T
        ag, ah = a[edges[:, 0]], a[edges[:, 1]]
        bg, bh = self.basis[self.g], self.basis[self.h]
        self.offset = 0.5 * (ag + ah)
        self.X = np.column_stack((np.ones(len(edges)), 0.5 * (ag + ah),
                                  0.5 * (ag * ag + ah * ah),
                                  bg + bh,
                                  0.5 * (bg * ah[:, None] + bh * ag[:, None]),
                                  0.5 * (bg * (ah * ah)[:, None]
                                         + bh * (ag * ag)[:, None])))
        self.p = self.X.shape[1]
        self.penalty = np.zeros(self.p)
        self.penalty[3:] = self.alpha / self.s

    def ridge(self, deadline=None):
        _deadline_check(deadline)
        m = len(self.y)
        penalized = np.flatnonzero(self.penalty)
        A = np.vstack((self.X / np.sqrt(m),
                       np.diag(np.sqrt(self.penalty))[penalized]))
        b = np.concatenate(((self.y - self.offset) / np.sqrt(m),
                            np.zeros(len(penalized))))
        coef = linalg.lstsq(A, b, lapack_driver="gelsy", check_finite=False)[0]
        _deadline_check(deadline)
        return coef

    def unpack_latent(self, raw):
        centered = raw - raw.mean(axis=0, keepdims=True)
        rms = np.sqrt(np.mean(centered * centered, axis=0))
        active = rms > np.sqrt(np.finfo(float).tiny)
        u = np.zeros_like(centered)
        u[:, active] = centered[:, active] / rms[active]
        return u, rms, active

    def objective_gradient(self, params):
        coef = params[:self.p]
        raw = params[self.p:self.p + self.s * self.rank].reshape(self.s, self.rank)
        lam = params[self.p + self.s * self.rank:]
        u, rms, active = self.unpack_latent(raw)
        pair_products = u[self.g] * u[self.h]
        residual = self.offset + self.X @ coef + pair_products @ lam - self.y
        loss = (np.mean(residual * residual) + np.dot(self.penalty, coef * coef)
                + self.eta * np.dot(lam, lam))
        dr = 2.0 * residual / len(self.y)
        gc = self.X.T @ dr + 2.0 * self.penalty * coef
        gl = pair_products.T @ dr + 2.0 * self.eta * lam
        gu = np.zeros_like(u)
        np.add.at(gu, self.g, dr[:, None] * u[self.h] * lam)
        np.add.at(gu, self.h, dr[:, None] * u[self.g] * lam)
        gv = gu - gu.mean(axis=0, keepdims=True) - u * np.mean(gu * u, axis=0)
        gv[:, active] /= rms[active]
        gv[:, ~active] = 0.0
        return float(loss), np.concatenate((gc, gv.ravel(), gl))

    def make_model(self, coef, raw=None, lam=None, diagnostics=None):
        n = len(self.a)
        e, db, dc = (np.zeros(n) for _ in range(3))
        e[self.supported] = self.basis @ coef[3:3 + self.k]
        db[self.supported] = self.basis @ coef[3 + self.k:3 + 2 * self.k]
        dc[self.supported] = self.basis @ coef[3 + 2 * self.k:]
        u = np.zeros((n, self.rank))
        if self.rank:
            u[self.supported] = self.unpack_latent(raw)[0]
        return FittedModel(a=self.a.copy(), supported=self.supported.copy(),
                           family="horlbeck", rank=self.rank,
                           alpha=self.alpha, eta=self.eta,
                           e0=float(coef[0]), beta=float(coef[1]), c=float(coef[2]),
                           e=e, d_beta=db, d_c=dc, u=u,
                           lambdas=np.asarray(lam if lam is not None else []),
                           diagnostics=diagnostics or {})


def fit_model(a, edges, y, rank, alpha, eta=0.0, deadline=None):
    """Fit one candidate; positive-rank failures are explicit, never zero gains.

    Non-grid penalties are permitted for numerical unit checks; the experiment
    runner is responsible for enumerating only the registered grids.
    """
    if rank not in RANKS or alpha < 0 or eta < 0:
        raise ValueError("Invalid rank or penalty")
    started = time.monotonic()
    _deadline_check(deadline)
    a, edges, y = _inputs(a, edges, y)
    design = _Design(a, edges, y, alpha, rank, eta)
    base_coef = design.ridge(deadline)
    if rank == 0:
        residual = design.offset + design.X @ base_coef - y
        objective = float(np.mean(residual * residual)
                          + np.dot(design.penalty, base_coef * base_coef))
        gradient = design.X.T @ (2 * residual / len(y)) + 2 * design.penalty * base_coef
        return design.make_model(base_coef, diagnostics={
            "converged": True, "method": "constrained_ridge_gelsy",
            "objective": objective, "gradient_inf": float(np.max(np.abs(gradient))),
            "iterations": 0, "seconds": time.monotonic() - started,
            "termination": "exact_constrained_least_squares", "starts": []})

    starts, candidates = [], []
    for seed in START_SEEDS:
        _deadline_check(deadline)
        start_time = time.monotonic()
        rng = np.random.default_rng(seed)
        raw = rng.standard_normal((design.s, rank))
        lam = 0.01 * rng.standard_normal(rank)
        params = np.concatenate((base_coef, raw.ravel(), lam))
        phases, total_iterations = [], 0

        def fg(x):
            _deadline_check(deadline)
            return design.objective_gradient(x)

        res = None
        for maxiter in (500, 1500):
            _deadline_check(deadline)
            res = minimize(fg, params, method="L-BFGS-B", jac=True,
                           options={"maxiter": maxiter, "ftol": 1e-10,
                                    "gtol": 1e-6, "maxfun": 100000})
            total_iterations += int(res.nit)
            finite = bool(np.isfinite(res.fun) and np.isfinite(res.x).all()
                          and np.isfinite(res.jac).all())
            phases.append({"iterations": int(res.nit), "objective": _diagnostic_float(res.fun),
                           "gradient_inf": _diagnostic_float(np.max(np.abs(res.jac))),
                           "success": bool(res.success), "status": int(res.status),
                           "termination": str(res.message)})
            if res.success or not finite:
                break
            params = res.x  # Continue this very same start, never reinitialize.
        info = {"seed": seed, "converged": bool(res.success and finite),
                "iterations": total_iterations, "seconds": time.monotonic() - start_time,
                "objective": _diagnostic_float(res.fun),
                "gradient_inf": _diagnostic_float(np.max(np.abs(res.jac))),
                "termination": str(res.message), "phases": phases}
        starts.append(info)
        if info["converged"]:
            candidates.append((float(res.fun), seed, res.x.copy()))
    _deadline_check(deadline)
    diagnostics = {"converged": bool(candidates), "method": "L-BFGS-B",
                   "seconds": time.monotonic() - started, "starts": starts}
    if not candidates:
        raise FitConvergenceError("All three positive-rank starts failed", diagnostics)
    objective, best_seed, params = min(candidates, key=lambda item: (item[0], item[1]))
    winning = next(x for x in starts if x["seed"] == best_seed)
    diagnostics.update({"objective": objective, "selected_start": best_seed,
                        "gradient_inf": winning["gradient_inf"],
                        "iterations": winning["iterations"],
                        "termination": winning["termination"]})
    return design.make_model(params[:design.p],
                             params[design.p:design.p + design.s * rank].reshape(design.s, rank),
                             params[design.p + design.s * rank:], diagnostics)


def fit_control(a, edges, y, kind, alpha=0.0, deadline=None):
    """Appendix B8 descriptive controls; caller performs training-only alpha CV."""
    started = time.monotonic()
    _deadline_check(deadline)
    a, edges, y = _inputs(a, edges, y)
    if kind not in ("additive", "endpoint", "shared_quadratic") or alpha < 0:
        raise ValueError("Invalid control family or penalty")
    supported = np.unique(edges)
    model = FittedModel(a=a.copy(), supported=supported, family=kind,
                        rank=0, alpha=float(alpha), eta=0.0,
                        e=np.zeros(len(a)), d_beta=np.zeros(len(a)), d_c=np.zeros(len(a)),
                        u=np.empty((len(a), 0)))
    if kind != "additive":
        basis = np.zeros((len(a), len(supported) - 1))
        basis[supported] = _mean_zero_basis(len(supported))
        g, h = edges.T
        s = a[g] + a[h]
        if kind == "endpoint":
            shared = np.ones((len(edges), 1))
            offset = s
        else:
            shared = np.column_stack((np.ones(len(edges)), s, s * s))
            offset = np.zeros(len(edges))
        X = np.column_stack((shared, basis[g] + basis[h]))
        p0 = shared.shape[1]
        penalty = np.zeros(X.shape[1])
        penalty[p0:] = alpha / len(supported)
        idx = np.flatnonzero(penalty)
        A = np.vstack((X / np.sqrt(len(y)), np.diag(np.sqrt(penalty))[idx]))
        b = np.concatenate(((y - offset) / np.sqrt(len(y)), np.zeros(len(idx))))
        coef = linalg.lstsq(A, b, lapack_driver="gelsy", check_finite=False)[0]
        model.e = basis @ coef[p0:]
        if kind == "endpoint":
            model.e0 = float(coef[0])
        else:
            model.shared_theta = coef[:p0]
    _deadline_check(deadline)
    residual = model.predict(edges) - y
    objective = float(np.mean(residual * residual) +
                      alpha * np.mean(model.e[supported] ** 2))
    model.diagnostics = {"converged": True, "method": "constrained_ridge_gelsy",
                         "objective": objective, "iterations": 0,
                         "seconds": time.monotonic() - started,
                         "termination": "exact_control", "starts": []}
    return model
