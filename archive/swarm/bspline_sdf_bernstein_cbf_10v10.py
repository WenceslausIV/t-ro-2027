"""
B-spline SDF + Bernstein-coefficient CBF for swarm-vs-swarm shape avoidance (10 vs 10).

Robots are rigid polygons (convex or non-convex, translation only); a team's SDF target is the
exact union SDF  min_i sdf_i(x - x_i)  (--shapes disk uses 32-gons).

Changes from the BP-SDF version (timevaryingBP lse+enclosure+innercollisionavoid.py):

1) Shape representation: tensor-product *uniform cubic B-spline* SDF instead of a global
   tensor Bernstein polynomial. Local support (K cells per axis) gives much finer detail
   than a single degree-(n-1) Bernstein polynomial over the whole domain.

2) Safety condition: no closest-point sampling / level-set projection / LSE.
   Following the Bernstein-coefficient idea of the ICRA'27 semantic-CBF paper, every
   B-spline cell is a bicubic polynomial that is rewritten *exactly* in the Bernstein
   (Bezier) basis with the b2b matrix Q (paper eq. b2b, applied in x and y).

   Team shapes:   S_A = {x | phi_A(x) <= l_A},   S_B = {x | phi_B(x) <= l_B}
   Barrier:       g(x) = phi_A(x) + phi_B(x) - l_A - l_B - m
   If g(x) > 0 for all x then S_A and S_B are disjoint (a common point would give g <= -m).
   For true SDFs, min_x g = dist(S_A, S_B) - m, so g is exactly a "gap between the shapes".

   g is linear in the control weights W_A + W_B, so on every cell
       g(x) = sum_{i,j} beta_ij B_i^3(xi) B_j^3(eta),   beta = Q (W_A+W_B)_cell Q^T - (l_A+l_B+m)
   and the CBF polynomial Psi = g_dot + gamma*g has Bernstein coefficients
       beta_dot + gamma*beta,
   which are *affine in all robots' inputs* (W_A depends on robot positions through the
   fitted SDF, see BSplineSDF2D.fit_team).  Requiring all of them >= 0 gives
   Psi >= 0 on the whole cell (convex combination, paper Thm. 1) -> finite affine constraints.

3) One joint QP for all 20 robots:
       min ||u - u_nom||^2
       s.t. Bernstein coefficient constraints (inter-team shape CBF)
            pairwise bounding-circle CBF inside each team
            ||u_i|| <= u_max (inscribed polygon)
   solved exactly as a least-distance program with NNLS (Lawson-Hanson), no extra deps.
"""
import argparse
import time

import numpy as np
import torch
from numpy.lib.stride_tricks import sliding_window_view
from scipy.optimize import nnls

# Cubic B-spline span -> cubic Bezier control points (paper eq. b2b)
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


# ---------- robot shapes: simple polygons in the body frame (convex or non-convex, CCW) ----------
def regular_polygon(n, r, phase=0.0):
    a = phase + 2 * np.pi * np.arange(n) / n
    return r * np.stack([np.cos(a), np.sin(a)], axis=1)


def star_polygon(n, r_out, r_in, phase=np.pi / 2):
    a = phase + np.pi * np.arange(2 * n) / n
    r = np.where(np.arange(2 * n) % 2 == 0, r_out, r_in)
    return np.stack([r * np.cos(a), r * np.sin(a)], axis=1)


def centered(P):
    P = np.asarray(P, dtype=np.float64)
    return P - P.mean(axis=0)


SHAPES = {
    'disk':     regular_polygon(32, 0.3),
    'square':   regular_polygon(4, 0.32, np.pi / 4),
    'triangle': regular_polygon(3, 0.36, np.pi / 2),
    'hexagon':  regular_polygon(6, 0.3),
    'star':     star_polygon(5, 0.36, 0.16),
    'L':        centered([[0, 0], [0.5, 0], [0.5, 0.18], [0.18, 0.18], [0.18, 0.5], [0, 0.5]]),
    'cross':    centered([[0.17, 0], [0.33, 0], [0.33, 0.17], [0.5, 0.17], [0.5, 0.33], [0.33, 0.33],
                          [0.33, 0.5], [0.17, 0.5], [0.17, 0.33], [0, 0.33], [0, 0.17], [0.17, 0.17]]),
    'U':        centered([[0, 0], [0.5, 0], [0.5, 0.45], [0.34, 0.45], [0.34, 0.16], [0.16, 0.16],
                          [0.16, 0.45], [0, 0.45]]),
}


def rotate(P, th):
    c, s = np.cos(th), np.sin(th)
    return P @ np.array([[c, s], [-s, c]])


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


class RobotShapes:
    """N rigid polygonal robots (translation only). Holds padded edges and body-frame samples."""

    def __init__(self, polys, device, ds=0.02):
        self.polys = [np.asarray(p, dtype=np.float64) for p in polys]
        E = max(len(p) for p in self.polys)
        A = np.zeros((len(polys), E, 2))
        B = np.zeros((len(polys), E, 2))
        for n, p in enumerate(self.polys):
            q = np.roll(p, -1, axis=0)
            A[n, :len(p)], B[n, :len(p)] = p, q
            A[n, len(p):] = B[n, len(p):] = p[-1]                          # degenerate padding
        self.A = torch.as_tensor(A, device=device)
        self.B = torch.as_tensor(B, device=device)
        self.device = device
        self.bound_r = np.array([np.linalg.norm(p, axis=1).max() for p in self.polys])

        # body-frame samples: boundary every ds, plus an interior grid (for the level set l)
        bnd, own_b, inner, own_i = [], [], [], []
        for n, p in enumerate(self.polys):
            q = np.roll(p, -1, axis=0)
            for a, b in zip(p, q):
                m = max(int(np.ceil(np.linalg.norm(b - a) / ds)), 1)
                s = np.arange(m)[:, None] / m
                bnd.append(a + s * (b - a))
                own_b += [n] * m
            g = np.arange(-self.bound_r[n], self.bound_r[n], 2 * ds)
            gx, gy = np.meshgrid(g, g)
            cand = np.stack([gx.ravel(), gy.ravel()], axis=1)
            sd, _ = self.sdf_body(cand, n)
            inner.append(cand[sd < 0])
            own_i += [n] * int((sd < 0).sum())
        self.boundary = np.concatenate(bnd)
        self.boundary_owner = np.array(own_b)
        self.samples = np.concatenate([self.boundary] + inner)
        self.samples_owner = np.concatenate([self.boundary_owner, own_i]).astype(int)

    def sdf_body(self, pts, n):
        P = torch.as_tensor(pts, device=self.device)[:, None, :]
        sd, gr = polygon_sdf(P, self.A[n:n + 1], self.B[n:n + 1])
        return sd[:, 0].cpu().numpy(), gr[:, 0].cpu().numpy()

    def team_sdf(self, pts, pos):
        """Union SDF min_i sdf_i(p - x_i) at pts (G, 2) (torch); returns sdf, nearest idx, gradient."""
        p = torch.as_tensor(pos, device=self.device)
        sd, gr = polygon_sdf(pts[:, None, :] - p[None], self.A, self.B)   # (G, N), (G, N, 2)
        val, idx = sd.min(dim=1)
        grad = gr[torch.arange(len(idx), device=self.device), idx]
        return val, idx, grad


class BSplineSDF2D:
    """
    phi(x, y) = sum_{a,b} W[a, b] b_a(x) b_b(y)  on a square domain with K uniform cells per axis
    (K+3 control weights per axis).  The fit uses a fixed sample grid, so the ridge solution is
    a fixed linear map  W = M T M^T  of the grid SDF targets T (separable ridge regularisation).
    """

    def __init__(self, domain_min, domain_max, n_cells=48, n_samples=160, lam=1e-3):
        self.lo = float(domain_min)
        self.hi = float(domain_max)
        self.K = n_cells
        self.h = (self.hi - self.lo) / n_cells
        self.n_ctrl = n_cells + 3

        self.g1d = self.lo + (np.arange(n_samples) + 0.5) * (self.hi - self.lo) / n_samples
        gx, gy = np.meshgrid(self.g1d, self.g1d, indexing='ij')
        self.grid_pts = np.stack([gx.ravel(), gy.ravel()], axis=1)       # (ns*ns, 2)
        self.ns = n_samples

        Phi = self.design(self.g1d)                                       # (ns, K+3)
        A = Phi.T @ Phi + lam * np.eye(self.n_ctrl)
        self.M = np.linalg.solve(A, Phi.T)                                # (K+3, ns)

        # Offline rate constant for lossless pruning.  Cell (a,b), coefficient (i,j):
        #   beta_dot = sum_g K(g) T_dot(g),  K = (Q M_a)_i (x) (Q M_b)_j,  |T_dot(g)| <= u_max
        #   => |beta_dot| <= kappa * u_max,   kappa = ||(Q M_a)_i||_1 ||(Q M_b)_j||_1
        M_win = np.moveaxis(sliding_window_view(self.M, 4, axis=0), -1, 1)   # (K, 4, ns)
        r = np.abs(np.einsum('ip,apq->aiq', Q_B2B, M_win)).sum(axis=2)       # (K, 4)
        self.kappa = np.einsum('ai,bj->abij', r, r)                          # (K, K, 4, 4)

        self.dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.M_t = torch.as_tensor(self.M, device=self.dev)
        self.grid_t = torch.as_tensor(self.grid_pts, device=self.dev)

    # ---------- basis ----------
    def _cell(self, x):
        t = (np.asarray(x, dtype=np.float64) - self.lo) / self.h
        c = np.clip(np.floor(t).astype(int), 0, self.K - 1)
        xi = np.clip(t - c, 0.0, 1.0)
        return c, xi

    def design(self, x, deriv=False):
        """Dense 1-D design matrix (n, K+3) (or its x-derivative)."""
        c, xi = self._cell(x)
        B = cubic_bspline_local_deriv(xi) / self.h if deriv else cubic_bspline_local(xi)
        out = np.zeros((len(c), self.n_ctrl))
        rows = np.arange(len(c))
        for i in range(4):
            out[rows, c + i] = B[:, i]
        return out

    # ---------- evaluation ----------
    def predict(self, pts, W):
        Bx, By = self.design(pts[:, 0]), self.design(pts[:, 1])
        return ((Bx @ W) * By).sum(axis=1)

    def gradient(self, pts, W):
        Bx, By = self.design(pts[:, 0]), self.design(pts[:, 1])
        dBx, dBy = self.design(pts[:, 0], True), self.design(pts[:, 1], True)
        gx = ((dBx @ W) * By).sum(axis=1)
        gy = ((Bx @ W) * dBy).sum(axis=1)
        return np.stack([gx, gy], axis=1)

    def predict_grid(self, x1d, W):
        return self.design(x1d) @ W @ self.design(x1d).T                  # [ix, iy]

    def bernstein(self, W):
        """Exact Bezier (Bernstein) coefficients of every cell: (..., K, K, 4, 4)."""
        win = sliding_window_view(W, (4, 4), axis=(-2, -1))
        return np.einsum('ip,...abpq,jq->...abij', Q_B2B, win, Q_B2B, optimize=True)

    # ---------- team SDF fit + sensitivity to robot motion ----------
    def fit_team(self, shapes, pos):
        """
        Targets T(g) = min_i sdf_i(g - x_i) (union of the robots' polygons) on the fixed grid.
        Returns
            W       : (K+3, K+3) control weights (NumPy)
            dW      : (N, 2, K+3, K+3) torch tensor with  W_dot = sum_i dW[i] @ u_i
        dT(g)/dx_i = -grad sdf_i(g - x_i) on the grid points whose nearest robot is i.
        Float64 on the GPU when available.
        """
        g, M = self.grid_t, self.M_t
        T, idx, n = shapes.team_sdf(g, pos)
        W = M @ T.reshape(self.ns, self.ns) @ M.T

        N = pos.shape[0]
        S = torch.zeros(N, 2, g.shape[0], device=self.dev, dtype=torch.float64)
        S[idx, :, torch.arange(g.shape[0], device=self.dev)] = -n
        dW = (M @ S.reshape(N * 2, self.ns, self.ns) @ M.T).reshape(N, 2, self.n_ctrl, self.n_ctrl)
        return W.cpu().numpy(), dW


def adaptive_level(model, W, dW, shapes, pos, eps_active=0.01, max_active=32):
    """
    l = max of phi over the robots' polygons (boundary every 2 cm + interior grid): the smallest
    level whose sublevel set {phi <= l} contains every sampled robot point (no LSE).
    l is a max, so dl/dt = max over the active samples of d phi(q_k)/dt (Danskin), with
        d phi(q_k)/dt = grad phi(q_k) . u_owner(k)  +  psi(q_k)^T W_dot.
    Returns l and dv (n_active, N, 2) for every sample within eps_active of the max;
    the QP bounds all of them from above with an epigraph variable t >= dv_k . u.
    """
    N = pos.shape[0]
    owner = shapes.samples_owner
    pts = shapes.samples + pos[owner]

    Bx, By = model.design(pts[:, 0]), model.design(pts[:, 1])
    v = ((Bx @ W) * By).sum(axis=1)
    l = v.max()
    act = np.argsort(-v)[:max_active]
    act = act[v[act] >= l - eps_active]
    Bx, By = Bx[act], By[act]

    dv = np.zeros((len(act), N, 2))
    dv[np.arange(len(act)), owner[act]] = model.gradient(pts[act], W)     # moving sample points
    Bx_t = torch.as_tensor(Bx, device=dW.device)                          # field change
    By_t = torch.as_tensor(By, device=dW.device)
    dv += ((Bx_t[:, None, None, None, :] @ dW[None]).squeeze(-2) * By_t[:, None, None, :]
           ).sum(-1).cpu().numpy()
    return l, dv


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


def speed_polygon(n_robots, u_max, n_sides=12):
    """Rows for n_k^T u_i <= u_max cos(pi/P): polygon inscribed in the speed circle."""
    ang = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
    nk = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    G = np.zeros((n_robots * n_sides, 2 * n_robots))
    for i in range(n_robots):
        G[i * n_sides:(i + 1) * n_sides, 2 * i:2 * i + 2] = -nk
    h = -np.full(n_robots * n_sides, u_max * np.cos(np.pi / n_sides))
    return G, h


def intra_team_rows(pos, offset, n_total, radii, gamma_in, margin_in):
    """
    Pairwise CBF inside one team on bounding circles r_i of the polygons:
        2(x_i-x_j)^T(u_i-u_j) + gamma*(|x_i-x_j|^2-(r_i+r_j+m)^2) >= 0.
    """
    N = pos.shape[0]
    rows, rhs = [], []
    for i in range(N):
        for j in range(i + 1, N):
            d = pos[i] - pos[j]
            r = np.zeros(2 * n_total)
            r[2 * (offset + i):2 * (offset + i) + 2] = 2 * d
            r[2 * (offset + j):2 * (offset + j) + 2] = -2 * d
            rows.append(r)
            rhs.append(-gamma_in * (d @ d - (radii[i] + radii[j] + margin_in) ** 2))
    return np.array(rows), np.array(rhs)


def bernstein_shape_rows(model, W_A, dW_A, l_A, dv_A, W_B, dW_B, l_B, dv_B, margin, gamma, u_max):
    """
    Bernstein-coefficient CBF rows for g = phi_A + phi_B - l_A - l_B - margin over z = [u, t_A, t_B]:
        dbeta_field/du . u - t_A - t_B + gamma * beta >= 0    (every cell, all 16 coefficients)
        t_A >= dv_A[k] . u,  t_B >= dv_B[k] . u               (every active level-set sample)
    t_A, t_B upper-bound dl_A/dt, dl_B/dt, so the Bernstein rows hold for every argmax switch.
    Both t are boxed by |t| <= t_bar (t_bar = u_max max_k sum_i ||dv[k]_i||), which every
    admissible u already satisfies, so the box does not change the QP.

    Lossless two-stage pruning (a dropped row is satisfied by every admissible (u, t)):
      1) offline rate bound, needs beta only:
             gamma*beta > n_moving kappa u_max + t_bar_A + t_bar_B   (|beta_dot| <= this)
      2) exact per-row bound on the survivors:
             gamma*beta > u_max sum_i ||dbeta/du_i|| + t_bar_A + t_bar_B
    Returns (G_bern, h_bern, G_epi, h_epi, beta_all, n_stage1).
    """
    beta_all = model.bernstein(W_A + W_B) - (l_A + l_B + margin)                        # (K, K, 4, 4)
    N = dW_A.shape[0] + dW_B.shape[0]
    NA = dW_A.shape[0]
    tbar_A = u_max * np.linalg.norm(dv_A, axis=2).sum(axis=1).max()
    tbar_B = u_max * np.linalg.norm(dv_B, axis=2).sum(axis=1).max()

    # stage 1 (cell level: keep a cell if any of its 16 coefficients may be active);
    # each moving team contributes kappa*u_max (a static side, e.g. a maze, has no robots)
    n_moving = int(dW_A.shape[0] > 0) + int(dW_B.shape[0] > 0)
    cand = gamma * beta_all <= n_moving * u_max * model.kappa + tbar_A + tbar_B
    a, b = np.nonzero(cand.any(axis=(2, 3)))
    beta = beta_all[a, b].reshape(-1)                                                    # (16 C1,)
    dW = torch.cat([dW_A, dW_B], dim=0)                                                  # (N, 2, n, n)
    win = dW.unfold(2, 4, 1).unfold(3, 4, 1)[:, :, torch.as_tensor(a), torch.as_tensor(b)]
    Q = torch.as_tensor(Q_B2B, device=dW.device)
    Bd = (Q @ win @ Q.T).reshape(N, 2, -1)                                               # (N, 2, 16 C1)
    bound = (u_max * torch.linalg.norm(Bd, dim=1).sum(dim=0)).cpu().numpy()

    # stage 2
    keep = gamma * beta - (bound + tbar_A + tbar_B) <= 0
    Bd_keep = Bd[:, :, torch.as_tensor(keep, device=Bd.device)].cpu().numpy()
    G_bern = np.hstack([Bd_keep.reshape(2 * N, -1).T,
                        -np.ones((keep.sum(), 2))])                                      # (C_keep, 2N+2)
    h_bern = -gamma * beta[keep]

    rows = []
    for col, tbar in ((2 * N, tbar_A), (2 * N + 1, tbar_B)):                            # t <= t_bar
        r = np.zeros(2 * N + 2)
        r[col] = -1.0
        rows.append(r)
    h_box = [-tbar_A, -tbar_B]
    for k in range(len(dv_A)):
        r = np.zeros(2 * N + 2)
        r[:2 * NA] = -dv_A[k].ravel()
        r[2 * N] = 1.0
        rows.append(r)
    for k in range(len(dv_B)):
        r = np.zeros(2 * N + 2)
        r[2 * NA:2 * N] = -dv_B[k].ravel()
        r[2 * N + 1] = 1.0
        rows.append(r)
    G_epi = np.array(rows)
    h_epi = np.concatenate([h_box, np.zeros(len(rows) - 2)])
    return G_bern, h_bern, G_epi, h_epi, beta_all.reshape(-1), len(beta)


def robot_gap(shapes_A, pos_A, shapes_B, pos_B):
    """
    Min distance between any A polygon and any B polygon (negative if they overlap):
    boundary samples of one team (all vertices included) against the other team's exact SDF.
    For disjoint polygons the minimum is attained at a vertex, so this is exact.
    """
    dev = shapes_A.device
    pa = torch.as_tensor(shapes_A.boundary + pos_A[shapes_A.boundary_owner], device=dev)
    pb = torch.as_tensor(shapes_B.boundary + pos_B[shapes_B.boundary_owner], device=dev)
    return min(shapes_B.team_sdf(pa, pos_B)[0].min().item(),
               shapes_A.team_sdf(pb, pos_A)[0].min().item())


def contact_points(x1d, zA, zB, dist):
    """
    Where the team shapes {zA <= 0} and {zB <= 0} nearly touch: Z = max(zA, zB) is how far a
    point is from being inside both shapes. Each connected region with Z < dist gets one point,
    at its minimum of Z (the pinch point between the two boundaries).
    """
    from scipy import ndimage
    Z = np.maximum(zA, zB)
    lab, n = ndimage.label(Z < dist)
    pts = []
    for k in range(1, n + 1):
        idx = np.argmin(np.where(lab == k, Z, np.inf))
        ix, iy = np.unravel_index(idx, Z.shape)
        pts.append((x1d[ix], x1d[iy]))
    return pts


def formation(center, n_cols, n_rows, spacing):
    xs = (np.arange(n_cols) - (n_cols - 1) / 2) * spacing
    ys = (np.arange(n_rows) - (n_rows - 1) / 2) * spacing
    return np.array([[center[0] + x, center[1] + y] for y in ys for x in xs])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--frames', type=int, default=700)
    ap.add_argument('--headless', action='store_true', help='no window (metrics only)')
    ap.add_argument('--gif', type=int, default=0, help='save the first N frames as a GIF')
    ap.add_argument('--plot-every', type=int, default=1)
    ap.add_argument('--shapes', choices=['mixed', 'disk'], default='mixed')
    ap.add_argument('--cells', type=int, default=64, help='B-spline cells per axis')
    ap.add_argument('--samples', type=int, default=220, help='fit grid samples per axis')
    ap.add_argument('--out', default='bspline_bernstein_10v10.gif')
    ap.add_argument('--scale', type=float, default=1.0, help='robot shape size factor')
    ap.add_argument('--red-levels', action='store_true', help='draw red-team SDF level sets 0.2n, n=1..20')
    ap.add_argument('--contact-dist', type=float, default=0.1,
                    help='mark (star) every place where the two team shapes are closer than this [m]')
    args = ap.parse_args()

    import matplotlib
    if args.headless:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    np.random.seed(0)
    domain = 5.5
    model = BSplineSDF2D(-domain, domain, n_cells=args.cells, n_samples=args.samples, lam=1e-3)

    if args.shapes == 'mixed':
        names_A = ['star', 'L', 'square', 'cross', 'triangle', 'U', 'hexagon', 'star', 'L', 'cross']
        names_B = ['U', 'triangle', 'star', 'L', 'hexagon', 'cross', 'square', 'U', 'star', 'triangle']
    else:
        names_A = names_B = ['disk'] * 10
    rng = np.random.default_rng(1)
    shapes_A = RobotShapes([rotate(args.scale * SHAPES[s], rng.uniform(0, 2 * np.pi)) for s in names_A], model.dev)
    shapes_B = RobotShapes([rotate(args.scale * SHAPES[s], rng.uniform(0, 2 * np.pi)) for s in names_B], model.dev)

    spacing = 0.9 * args.scale
    pos_A = formation((-0.4, 3.6), 5, 2, spacing)
    goals_A = formation((-0.4, -3.6), 5, 2, spacing)
    pos_B = formation((0.4, -3.6), 5, 2, spacing)
    goals_B = formation((0.4, 3.6), 5, 2, spacing)
    NA, NB = len(pos_A), len(pos_B)
    N = NA + NB

    dt = 0.02
    u_max = 1.5
    k_goal = 2.0
    gamma = 6.0           # Bernstein shape CBF gain
    shape_margin = 0.05   # required gap between the two team shapes [m]
    gamma_in = 6.0        # intra-team circle CBF gain
    margin_in = 0.02
    eps_t = 1e-4          # objective weight on the dl/dt epigraph variables t_A, t_B

    G_spd, h_spd = speed_polygon(N, u_max)

    vis_res = 220
    x_vis = np.linspace(-domain, domain, vis_res)
    X, Y = np.meshgrid(x_vis, x_vis, indexing='ij')
    fine = np.linspace(-domain, domain, 440)

    if not args.headless:
        plt.ion()
    fig, ax = plt.subplots(figsize=(8, 8))
    try:
        fig.canvas.manager.set_window_title("B-spline SDF + Bernstein-coefficient CBF (10 vs 10)")
    except Exception:
        pass
    gif_frames = []

    worst = dict(min_robot_gap=np.inf, min_shape_sep=np.inf, min_beta=np.inf, max_slack=0.0)
    print(f"cells {model.K}x{model.K}, ctrl {model.n_ctrl}x{model.n_ctrl}, h = {model.h:.3f} m")

    z_prev = None
    for frame in range(args.frames):
        t0 = time.perf_counter()

        W_A, dW_A = model.fit_team(shapes_A, pos_A)
        W_B, dW_B = model.fit_team(shapes_B, pos_B)
        l_A, dv_A = adaptive_level(model, W_A, dW_A, shapes_A, pos_A)
        l_B, dv_B = adaptive_level(model, W_B, dW_B, shapes_B, pos_B)

        G_bern, h_bern, G_epi, h_epi, beta, n_stage1 = bernstein_shape_rows(
            model, W_A, dW_A, l_A, dv_A, W_B, dW_B, l_B, dv_B, shape_margin, gamma, u_max)
        G_inA, h_inA = intra_team_rows(pos_A, 0, N, shapes_A.bound_r, gamma_in, margin_in)
        G_inB, h_inB = intra_team_rows(pos_B, NA, N, shapes_B.bound_r, gamma_in, margin_in)

        pos = np.vstack([pos_A, pos_B])
        goals = np.vstack([goals_A, goals_B])
        d = goals - pos
        u_nom = k_goal * d
        sp = np.linalg.norm(u_nom, axis=1, keepdims=True)
        u_nom = np.where(sp > u_max, u_nom / np.maximum(sp, 1e-9) * u_max, u_nom)

        # z = [u, s_A, s_B] with t = s / sqrt(eps_t): tiny weight so t is ~free in the objective
        pad = lambda M: np.hstack([M, np.zeros((M.shape[0], 2))])
        G = np.vstack([G_bern, pad(G_inA), pad(G_inB), G_epi, pad(G_spd)])
        G[:, -2:] /= np.sqrt(eps_t)
        h = np.concatenate([h_bern, h_inA, h_inB, h_epi, h_spd])
        n_soft = len(h_bern) + len(h_inA) + len(h_inB)
        z, slack = solve_ldp_qp(np.concatenate([u_nom.ravel(), [0.0, 0.0]]), G, h, n_soft, z_prev)
        z_prev = z
        u = z[:2 * N].reshape(N, 2)

        pos_A = pos_A + u[:NA] * dt
        pos_B = pos_B + u[NA:] * dt
        t_ctrl = time.perf_counter() - t0

        # ---------- safety diagnostics (on the state *before* the step) ----------
        gap = robot_gap(shapes_A, pos[:NA], shapes_B, pos[NA:])  # exact polygon-polygon distance
        zA = model.predict_grid(fine, W_A) - l_A
        zB = model.predict_grid(fine, W_B) - l_B
        shape_sep = np.maximum(zA, zB).min()      # > 0  <=>  shapes do not touch (on a fine grid)
        worst['min_robot_gap'] = min(worst['min_robot_gap'], gap)
        worst['min_shape_sep'] = min(worst['min_shape_sep'], shape_sep)
        worst['min_beta'] = min(worst['min_beta'], beta.min())
        worst['max_slack'] = max(worst['max_slack'], slack)

        if frame % 25 == 0 or frame == args.frames - 1:
            print(f"frame {frame:4d} | cand {n_stage1:5d} rows {len(h_bern):4d} | min beta {beta.min():+.3f} | "
                  f"shape sep {shape_sep:+.3f} | robot gap {gap:+.3f} | slack {slack:.1e} | "
                  f"l_A {l_A:+.3f} l_B {l_B:+.3f} | ctrl {1e3 * t_ctrl:.1f} ms")

        want_gif = frame < args.gif
        if (not args.headless and frame % args.plot_every == 0) or want_gif:
            ax.clear()
            pA = model.predict_grid(x_vis, W_A)
            pB = model.predict_grid(x_vis, W_B)
            ax.contourf(X, Y, pA, levels=[pA.min() - 1, l_A], colors=['royalblue'], alpha=0.3)
            ax.contour(X, Y, pA, levels=[l_A], colors='blue', linewidths=1.5)
            ax.contourf(X, Y, pB, levels=[pB.min() - 1, l_B], colors=['crimson'], alpha=0.3)
            ax.contour(X, Y, pB, levels=[l_B], colors='red', linewidths=1.5)
            if args.red_levels:
                ax.contour(X, Y, pB, levels=0.2 * np.arange(1, 21), cmap='autumn',
                           linewidths=0.8, alpha=0.8)
            for cx, cy in contact_points(fine, zA, zB, args.contact_dist):
                ax.plot(cx, cy, marker='*', ms=16, color='gold', mec='black', mew=1.0, zorder=10)
            for poly, p in zip(shapes_A.polys, pos_A):
                ax.add_patch(Polygon(poly + p, closed=True, color='midnightblue', alpha=0.9))
            for poly, p in zip(shapes_B.polys, pos_B):
                ax.add_patch(Polygon(poly + p, closed=True, color='darkred', alpha=0.9))
            ax.plot(goals_A[:, 0], goals_A[:, 1], 'x', color='blue', ms=7)
            ax.plot(goals_B[:, 0], goals_B[:, 1], 'x', color='red', ms=7)
            ax.set_title(f"B-spline SDF + Bernstein CBF | frame {frame} | "
                         f"min beta {beta.min():+.3f} | shape gap {shape_sep:+.3f}")
            ax.set_xlim(-domain, domain)
            ax.set_ylim(-domain, domain)
            ax.set_aspect('equal')
            if not args.headless:
                plt.pause(0.001)
            if want_gif:
                from PIL import Image
                fig.canvas.draw()
                w, hh = fig.canvas.get_width_height()
                buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(hh, w, 4)
                gif_frames.append(Image.fromarray(buf, 'RGBA').convert('RGB'))
                if frame + 1 == args.gif:
                    gif_frames[0].save(args.out, save_all=True,
                                       append_images=gif_frames[1:], duration=40, loop=0)
                    print(f">>> GIF saved: {args.out} ({args.gif} frames)")

    err_A = np.linalg.norm(pos_A - goals_A, axis=1).max()
    err_B = np.linalg.norm(pos_B - goals_B, axis=1).max()
    print("\n=== summary ===")
    print(f"min robot-robot gap (A vs B)   : {worst['min_robot_gap']:+.4f} m")
    print(f"min shape separation (grid)    : {worst['min_shape_sep']:+.4f}  (>0: shapes never touch)")
    print(f"min Bernstein coefficient beta : {worst['min_beta']:+.4f}")
    print(f"max slack used                 : {worst['max_slack']:.2e}")
    print(f"max distance to goal  A / B    : {err_A:.3f} / {err_B:.3f} m")
    if not args.headless:
        plt.ioff()
        plt.show()


if __name__ == '__main__':
    main()
