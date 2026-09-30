"""
Shared utilities for the C-space SDF CBF scripts (cspace_*.py):
cubic B-spline / Bezier helpers, exact polygon SDF, and the NNLS least-distance QP solver.
(The swarm scripts in swarm/ keep their own copies and do not depend on this file.)
"""
import numpy as np
import torch
from scipy.optimize import nnls

# Cubic B-spline span -> cubic Bezier control points (ICRA'27 eq. b2b)
Q_B2B = np.array([[1, 4, 1, 0],
                  [0, 4, 2, 0],
                  [0, 2, 4, 0],
                  [0, 1, 4, 1]], dtype=np.float64) / 6.0


def cubic_bspline_local(xi):
    """[N0..N3](xi) of paper eq. (basis), xi in [0,1] -> (n, 4)."""
    return np.stack([(1 - xi) ** 3,
                     3 * xi ** 3 - 6 * xi ** 2 + 4,
                     -3 * xi ** 3 + 3 * xi ** 2 + 3 * xi + 1,
                     xi ** 3], axis=-1) / 6.0


def cubic_bspline_local_deriv(xi):
    return np.stack([-3 * (1 - xi) ** 2,
                     9 * xi ** 2 - 12 * xi,
                     -9 * xi ** 2 + 6 * xi + 3,
                     3 * xi ** 2], axis=-1) / 6.0


def polygon_sdf(P, A, B):
    """
    Signed distance of points P (..., N, 2) (already in each robot's frame) to N polygons with
    edges A[n, e] -> B[n, e] (N, E, 2); padded edges are degenerate (A == B).
    Unsigned distance = min over edges; sign by the even-odd rule, so non-convex shapes are exact.
    Returns sdf (..., N) and its spatial gradient (..., N, 2) (unit vector).
    """
    P = P[..., None, :]                                                   # (..., N, 1, 2)
    d = B - A
    ap = P - A
    t = ((ap * d).sum(-1) / (d * d).sum(-1).clamp_min(1e-12)).clamp(0.0, 1.0)
    diff = ap - t[..., None] * d                                          # (..., N, E, 2)
    dist2, k = (diff * diff).sum(-1).min(dim=-1)                          # (..., N)
    diff = torch.gather(diff, -2, k[..., None, None].expand(*k.shape, 1, 2)).squeeze(-2)

    py = P[..., 1]
    ya, yb = A[..., 1], B[..., 1]
    cond = (ya > py) != (yb > py)
    x_int = A[..., 0] + (py - ya) * d[..., 0] / torch.where(cond, yb - ya, torch.ones_like(ya))
    inside = ((cond & (P[..., 0] < x_int)).sum(-1) % 2) == 1
    sign = 1.0 - 2.0 * inside.to(P.dtype)
    dist = dist2.sqrt()
    return sign * dist, sign[..., None] * diff / dist.clamp_min(1e-9)[..., None]


def _ldp(z0, G, h):
    """
    min ||z - z0||^2  s.t.  G z >= h   (least-distance programming via NNLS,
    Lawson & Hanson ch. 23). Returns None if infeasible.
    """
    f = h - G @ z0
    if np.all(f <= 0):                                                    # z0 already feasible
        return z0.copy()
    n = z0.size
    E = np.vstack([G.T, f[None, :]])
    e = np.zeros(n + 1)
    e[-1] = 1.0
    lam, _ = nnls(E, e, maxiter=50 * (n + 1))
    r = E @ lam - e
    if abs(r[-1]) < 1e-9:
        return None
    return z0 - r[:-1] / r[-1]


def solve_ldp_qp(u_nom, G, h, n_soft, z_prev=None, rho=1e4, ws_size=200, tol=1e-9, max_rounds=10):
    """
    Hard QP  min ||u - u_nom||^2  s.t.  G u >= h.
    Only if it is infeasible: soft version with one slack s on the first n_soft rows,
    min ||u - u_nom||^2 + rho s^2. Returns (u, slack).

    Warm start by constraint generation: solve the QP on a working set of rows (the ws_size
    rows closest to tight at the previous solution z_prev), then add every violated row and
    repeat. A relaxation's minimiser that satisfies *all* rows is the minimiser of the full QP,
    so the result is exact; if a relaxation is infeasible, so is the full QP.
    """
    feasible = True
    if z_prev is not None:
        r = G @ z_prev - h
        work = np.zeros(len(h), dtype=bool)
        work[np.argsort(r)[:ws_size]] = True
        for _ in range(max_rounds):
            u = _ldp(u_nom, G[work], h[work])
            if u is None:
                feasible = False
                break
            viol = G @ u - h < -tol * np.maximum(1.0, np.abs(h))
            if not viol.any():
                return u, 0.0
            work |= viol
    if feasible:
        u = _ldp(u_nom, G, h)
        if u is not None:
            return u, 0.0
    col = np.zeros((G.shape[0], 1))
    col[:n_soft, 0] = 1.0 / np.sqrt(rho)                                  # scaled slack s' = sqrt(rho) s
    z = _ldp(np.concatenate([u_nom, [0.0]]), np.hstack([G, col]), h)
    if z is None:
        raise RuntimeError("QP infeasible (speed polygon + slack should always be feasible)")
    return z[:-1], max(z[-1] / np.sqrt(rho), 0.0)
