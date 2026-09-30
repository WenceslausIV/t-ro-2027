"""
Summed-field barriers with Bernstein coefficient constraints (Sec. V of tro/main.tex).

Safe set (continuous, no boxes):  H = { poses : phi_B(y(x)) >= l_B for every x on S_A = {phi_A = l_A} },
where y(x) maps a point of body A (its frame) into B's frame. Enforcement: on every box Q of the Bernstein cover
of S_A, the lifted barrier
        h~(x) = phi_A(x) - l_A + phi_B(y(x)) - l_B            (= phi_B - l_B on S_A)
must satisfy   Psi~(x, u) = d/dt phi_B(y(x)) + gamma h~(x) >= 0  for all x in Q.
phi_A and phi_B are replaced on Q by their second-order Taylor models at the box center with certified
third-derivative remainders (C^2 majorants M3), so Psi~ is bounded below by a quadratic polynomial in x plus
  - gamma (M3_A + M3_B) r^3 / 6                     (value remainder)
  - (M3_B r^2 / 2) sum_m L_m |u_m|                    (gradient remainder times the point speed; epigraph t_m >= |u_m|)
The 3^n tensor Bernstein coefficients (degree 2) of the quadratic are affine in u; all of them >= 0 is sufficient.
Rows are returned as  A u + T t + C >= 0  (T <= 0).

Frames: every input m moves the point x of Q, in B's frame, with velocity  ydot(x) = Jc[m] + W[m] d  (d = y - y_c),
W[m] the rotation generator (3D: [omega_m]_x, 2D: omega_m J).
"""
import itertools

import numpy as np

from proto3d import Spline3

GAMMA_DEFAULT = 5.
EPS_T = 1e-3                       # QP weight of the epigraph variables t (z = (u, sqrt(EPS_T) t))

_T2 = np.array([-.5, 0., .5])      # degree-2 Bernstein coefficients of (xi - 1/2) on [0, 1]
_S2 = np.array([.25, -.25, .25])   # ... of (xi - 1/2)^2


def bern_tables(dim):
    K = np.array(list(itertools.product(range(3), repeat=dim)))
    TK = _T2[K]                                            # (3^dim, dim)
    QT = TK[:, :, None] * TK[:, None, :]
    for a in range(dim):
        QT[:, a, a] = _S2[K[:, a]]
    return TK, QT


TABLES = {2: bern_tables(2), 3: bern_tables(3)}


def third_cellwise(f, C=None):
    """Certified cellwise bound on the Frobenius norm of the third-derivative tensor of f, from the Bernstein
    coefficients of every third partial derivative (C: cell coefficients, e.g. of a slab of f)."""
    C = f.cell_coeffs() if C is None else C
    d = lambda X, a: (X.shape[3 + a] - 1) * np.diff(X, axis=3 + a) / f.h
    s = 0.
    for a in range(3):
        for b in range(a, 3):
            for c in range(b, 3):
                m = np.abs(d(d(d(C, a), b), c)).max(axis=(3, 4, 5))
                mult = 6 / (1 + (a == b) + (b == c) + (a == c) + 2 * (a == b == c))
                s = s + mult * m ** 2
    return np.sqrt(s)


def third_majorant(f, cell=None):
    """C^2 cubic-spline majorant of the cellwise third-derivative bound: control value m is the maximum over
    cells m-4..m+1 per axis, so M3(c) bounds the tensor on every cell adjacent to the cell of c."""
    W = np.pad(third_cellwise(f) if cell is None else cell, 4, mode='edge')
    for a in range(3):
        W = np.lib.stride_tricks.sliding_window_view(W, 6, axis=a).max(axis=-1)
    g = Spline3.__new__(Spline3)
    g.lo, g.h, g.K, g.n, g.W = f.lo, f.h, f.K, f.n, W
    return g


class LinearMajorant:
    """Lipschitz (piecewise-trilinear) majorant of a cellwise bound on a node grid of spacing h / refine:
    the value at node j is the maximum of the cellwise bound over all cells that meet the cube of half-width
    h / refine + reach around the node. For any point c, trilinear interpolation is a convex combination of
    the corner nodes of c's fine cell, each of whose cubes contains the ball of radius `reach` around c, so
    M(c) bounds the cellwise quantity on that ball."""

    def __init__(self, f, cell, reach, refine=2):
        self.lo, self.hf = np.asarray(f.lo, float), f.h / refine
        fine = cell
        for a in range(3):
            fine = np.repeat(fine, refine, axis=a)
        m = int(np.ceil((self.hf + reach) / self.hf - 1e-12))
        W = np.pad(fine, m, mode='edge')
        for a in range(3):
            W = np.lib.stride_tricks.sliding_window_view(W, 2 * m, axis=a).max(axis=-1)
        self.W = W                                           # node values, shape (Kf + 1,) per axis
        self.K = np.array(W.shape) - 1

    def eval(self, P):
        P = np.atleast_2d(P)
        s = (P - self.lo) / self.hf
        c = np.clip(np.floor(s).astype(int), 0, self.K - 1)
        t = np.clip(s - c, 0, 1)
        out = 0.
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = ((t[:, 0] if dx else 1 - t[:, 0]) * (t[:, 1] if dy else 1 - t[:, 1])
                         * (t[:, 2] if dz else 1 - t[:, 2]))
                    out = out + w * self.W[c[:, 0] + dx, c[:, 1] + dy, c[:, 2] + dz]
        return out


def field_bound(f, reach, cell=None):
    """Online third-derivative bound of a field (B role): Lipschitz majorant with window of about two cells."""
    return LinearMajorant(f, third_cellwise(f) if cell is None else cell, reach)


def box_bound(f, cell, P, r):
    """Offline third-derivative bound on each box (A role): maximum of the cellwise bound over the (at most
    2^3) cells that meet the cube of half-width r around each center P."""
    lo = np.clip(np.floor((P - r - f.lo) / f.h).astype(int), 0, f.K - 1)
    hi = np.clip(np.floor((P + r - f.lo) / f.h).astype(int), 0, f.K - 1)
    assert np.all(hi - lo <= 1), 'boxes must be smaller than the cells'
    out = np.zeros(len(P))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                i = np.minimum(lo + np.array([dx, dy, dz]), hi)
                out = np.maximum(out, cell[i[:, 0], i[:, 1], i[:, 2]])
    return out


def prepare_body(f, level, PC, side, dim, M3=None):
    """Offline Taylor data of phi_A - l_A at the box centers PC (A's frame) and certified bounds on each box."""
    P = PC if PC.shape[1] == 3 else np.c_[PC, np.zeros(len(PC))]
    v, g, H = f.eval(P, order=2)
    r = side * np.sqrt(dim) / 2
    m3 = box_bound(f, third_cellwise(f), P, r)
    g, H = g[:, :dim], H[:, :dim, :dim]
    TK, QT = TABLES[dim]
    # certified lower bound of phi_A - l_A on each box (Bernstein coefficients of its quadratic model - remainder)
    cA = (v - level)[:, None] + side * g @ TK.T + .5 * side ** 2 * np.einsum('kab,nab->nk', QT, H)
    lowA = cA.min(axis=1) - m3 * r ** 3 / 6
    return dict(PC=PC[:, :dim], side=side, r=r, dim=dim, vA=v - level, gA=g, HA=H, m3A=m3, lowA=lowA)


def rows(bodyA, idx, RbA, yc, fB, lB, M3B, Jc, W, Sb, eta, umax, auxmax, gamma=GAMMA_DEFAULT, prune_rows=True):
    """Coefficient rows for the boxes idx of A (centers yc in B's frame, rotation RbA: A frame -> B frame).
    Jc: (n, m, dim) center velocity per input (B frame); W: (m, dim, dim) rotation generators (B frame).
    Sb: (n, k) speed coefficients, |ydot(x)| <= Sb . a on the box for auxiliary variables a >= |E u|
    (a relative twist of the two bodies, componentwise); umax (m,), auxmax (k,): bounds of |u| and a over the
    admissible inputs (used only to drop rows that no admissible input can violate).
    Returns A (N, m), T (N, k), C (N,) with rows A u + T a + C >= 0, boxes kept, and the minimum over the
    kept boxes of the certified lower bound of h~."""
    dim, s, r = bodyA['dim'], bodyA['side'], bodyA['r']
    TK, QT = TABLES[dim]
    P3 = yc if dim == 3 else np.c_[yc, np.zeros(len(yc))]
    v, g, H = fB.eval(P3, order=2)
    m3 = M3B.eval(P3)
    g, H = g[:, :dim], H[:, :dim, :dim]
    R = RbA[:dim, :dim]
    # value part, in A's box coordinates e = x - c (so d = R e)
    b = bodyA['gA'][idx] + g @ R                                         # (n, dim)
    Q = .5 * (bodyA['HA'][idx] + np.einsum('ai,nab,bj->nij', R, H, R))
    val = ((bodyA['vA'][idx] + v - lB - (bodyA['m3A'][idx] + m3) * r ** 3 / 6)[:, None]
           + s * b @ TK.T + s ** 2 * np.einsum('kab,nab->nk', QT, Q))    # (n, K): lower bound of h~
    keep = val.min(axis=1) < eta
    if not keep.any():
        return None
    val, g, H, Jc, m3, Sb = val[keep], g[keep], H[keep], Jc[keep], m3[keep], Sb[keep]
    # velocity part: (g + H d)^T (Jc + W d), d = R e
    lin = np.einsum('nab,nmb->nma', H, Jc) + np.einsum('mba,nb->nma', W, g)            # (n, m, dim), in d
    quad = np.einsum('nab,mbc->nmac', H, W)
    quad = .5 * (quad + quad.transpose(0, 1, 3, 2))
    lin = lin @ R                                                                       # -> e coordinates
    quad = np.einsum('ai,nmab,bj->nmij', R, quad, R)
    Acoef = (np.einsum('nd,nmd->nm', g, Jc)[:, None, :] + s * np.einsum('kd,nmd->nkm', TK, lin)
             + s ** 2 * np.einsum('kab,nmab->nkm', QT, quad))                            # (n, K, m)
    # gradient remainder (M3 r^2 / 2) |ydot(x)| <= (M3 r^2 / 2) Sb . a
    n, K, m = Acoef.shape
    Tcoef = np.repeat(-.5 * (m3 * r ** 2)[:, None] * Sb, K, axis=0)                   # (n*K, k)
    C = (gamma * val).ravel()
    Acoef = Acoef.reshape(-1, m)
    if prune_rows:                                        # rows that some admissible input can violate
        need = C < np.abs(Acoef) @ umax - Tcoef @ auxmax
    else:
        need = np.ones(len(C), bool)
    return Acoef[need], Tcoef[need], C[need], np.flatnonzero(keep), float(val.min())


def solve(u_nom, A, T, C, E, G_in, h_in, z_prev=None):
    """min |u - u_nom|^2 + EPS_T |a|^2  s.t.  A u + T a + C >= 0,  a >= |E u|,  G_in u >= h_in.
    The feasible set in u equals that of A u + T |E u| + C >= 0 (T <= 0). Returns u, slack, z (warm start)."""
    from sdf_cbf_utils import solve_ldp_qp
    m, k = len(u_nom), len(E)
    se = np.sqrt(EPS_T)
    I = np.eye(k)
    G = np.vstack([np.c_[A, T / se], np.c_[-E, I / se], np.c_[E, I / se], np.c_[G_in, np.zeros((len(G_in), k))]])
    h = np.r_[-C, np.zeros(2 * k), h_in]
    zp = z_prev if z_prev is not None and len(z_prev) == m + k else None
    z, s = solve_ldp_qp(np.r_[u_nom, np.zeros(k)], G, h, len(C), zp)
    return z[:m], s, z
