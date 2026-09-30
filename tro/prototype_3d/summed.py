"""
Summed-field barriers with Bernstein coefficient constraints (Sec. V of paper/main.tex).

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
import os

import numpy as np
import scipy.sparse as sp

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
    """Offline third-derivative bound (A role): maximum over all cells meeting
    the cube of half-width r around each center P, including patch-size covers."""
    lo = np.clip(np.floor((P - r - f.lo) / f.h).astype(int), 0, f.K - 1)
    hi = np.clip(np.floor((P + r - f.lo) / f.h).astype(int), 0, f.K - 1)
    out = np.zeros(len(P))
    if not len(P):
        return out
    counts = (hi - lo).max(axis=0) + 1
    for offset in itertools.product(*(range(int(n)) for n in counts)):
        i = np.minimum(lo + np.asarray(offset), hi)
        out = np.maximum(out, cell[i[:, 0], i[:, 1], i[:, 2]])
    return out


def _quad_model(f, level, P, dim, r, cell3):
    """phi - level at the points P (A's frame): value, gradient, Hessian (restricted to dim) and the
    third-derivative bound on the boxes of half-diagonal r around them."""
    P3 = P if P.shape[1] == 3 else np.c_[P, np.zeros(len(P))]
    v, g, H = f.eval(P3, order=2)
    return v - level, g[:, :dim], H[:, :dim, :dim], box_bound(f, cell3, P3, r)


def prepare_body(f, level, PC, side, dim, M3=None):
    """Offline Taylor data of phi_A - l_A at the box centers PC (A's frame) and certified bounds on each box.
    The field, level, and cellwise third-derivative bound are kept for the online refinement of boxes."""
    r = side * np.sqrt(dim) / 2
    cell3 = third_cellwise(f)
    v, g, H, m3 = _quad_model(f, level, PC, dim, r, cell3)
    TK, QT = TABLES[dim]
    # certified lower bound of phi_A - l_A on each box (Bernstein coefficients of its quadratic model - remainder)
    cA = v[:, None] + side * g @ TK.T + .5 * side ** 2 * np.einsum('kab,nab->nk', QT, H)
    lowA = cA.min(axis=1) - m3 * r ** 3 / 6
    return dict(PC=PC[:, :dim], side=side, r=r, dim=dim, vA=v, gA=g, HA=H, m3A=m3, lowA=lowA,
                f=f, level=level, cell3=cell3)


REFINE_THETA = .1        # refine an active box while its gradient remainder M3_B r^2 / 2 exceeds this
REFINE_DEPTH = int(os.environ.get('SUMMED_REFINE_DEPTH', 0))  # native patches by default; subdivision is opt-in
REFINE_NEAR = float(os.environ.get('SUMMED_REFINE_NEAR', .01))   # ... only if its lower bound of h~ is below this [m]


def _layer(D, s, r, R, fB, lB, M3B, dim, gamma):
    """Lower bound of h~ (Bernstein coefficients of the quadratic model minus the value remainders) on the
    boxes of D, and the B-side quantities needed for their velocity rows."""
    TK, QT = TABLES[dim]
    P3 = D['yc'] if dim == 3 else np.c_[D['yc'], np.zeros(len(D['yc']))]
    v, g, H = fB.eval(P3, order=2)
    m3 = M3B.eval(P3)
    g, H = g[:, :dim], H[:, :dim, :dim]
    b = D['gA'] + g @ R                                                 # value part in A's box coordinates
    Q = .5 * (D['HA'] + np.einsum('ai,nab,bj->nij', R, H, R))
    val = ((D['vA'] + v - lB - (D['m3A'] + m3) * r ** 3 / 6)[:, None]
           + s * b @ TK.T + s ** 2 * np.einsum('kab,nab->nk', QT, Q))    # (n, K)
    return val, g, H, m3


def _velocity_rows(val, g, H, m3, Jc, W, Sb, s, r, R, dim, gamma):
    TK, QT = TABLES[dim]
    lin = np.einsum('nab,nmb->nma', H, Jc) + np.einsum('mba,nb->nma', W, g)            # (n, m, dim), in d
    quad = np.einsum('nab,mbc->nmac', H, W)
    quad = .5 * (quad + quad.transpose(0, 1, 3, 2))
    lin = lin @ R                                                                       # -> e coordinates
    quad = np.einsum('ai,nmab,bj->nmij', R, quad, R)
    Acoef = (np.einsum('nd,nmd->nm', g, Jc)[:, None, :] + s * np.einsum('kd,nmd->nkm', TK, lin)
             + s ** 2 * np.einsum('kab,nmab->nkm', QT, quad))                            # (n, K, m)
    n, K, m = Acoef.shape
    Tcoef = np.repeat(-.5 * (m3 * r ** 2)[:, None] * Sb, K, axis=0)                   # (n*K, k)
    return Acoef.reshape(-1, m), Tcoef, (gamma * val).ravel()


ROW_MODE = os.environ.get('SUMMED_ROW', 'vertex')
MULT_MODE = os.environ.get('SUMMED_MULT', 'one')
ROW_MODES = ('vertex', 'bernstein', 'single')
MULT_MODES = ('one', 'norm', 'proj', 'free')
W_MAX = 3.                         # 'free' multipliers: one variable w in [0, W_MAX] per active box
EPS_W = 1e-3                       # QP weight of (w - 1)


def split_values(D, val, s, r, dim):
    """Split the Bernstein coefficients of the unit-lift value bound into the part of phi_A - l_A (with its
    remainder, which scales with the multiplier) and the part of phi_B - l_B (with its remainder)."""
    TK, QT = TABLES[dim]
    valA = ((D['vA'] - D['m3A'] * r ** 3 / 6)[:, None] + s * D['gA'] @ TK.T
            + .5 * s ** 2 * np.einsum('kab,nab->nk', QT, D['HA']))
    return valA, val - valA


def multiplier_weights(gA, gB_local, mode):
    """Cell-constant w=lambda/gamma. This modifies a certificate, not the safe set.

    A small squared-gradient regularizer keeps the weights continuous at gA=0.
    Projection weights may be negative; their Taylor error MUST use abs(w).
    No cancellation or second-order claim is made for the full CBF residual.
    """
    if mode not in MULT_MODES or mode == 'free':
        raise ValueError(f'Multiplier mode {mode} has no fixed weights')
    if mode == 'one':
        return np.ones(len(gA))
    aa = np.einsum('ni,ni->n', gA, gA) + 1e-12
    if mode == 'norm':
        return np.sqrt((np.einsum('ni,ni->n', gB_local, gB_local) + 1e-12) / aa)
    return -np.einsum('ni,ni->n', gA, gB_local) / aa


def weighted_values(val, D, g, R, s, r, dim, mode):
    """Change the unit-lift value bound to w*g_A+h with certified signed weights.

    Screening and refinement continue to use the unit lift, so every variant has
    the same candidate cells at the same pose. On g_A=0 both lifts equal h.
    """
    if mode == 'one':
        return val
    w = multiplier_weights(D['gA'], g @ R, mode)
    TK, QT = TABLES[dim]
    qA = (D['vA'][:, None] + s * D['gA'] @ TK.T
          + .5 * s ** 2 * np.einsum('kab,nab->nk', QT, D['HA']))
    errA = D['m3A'] * r ** 3 / 6
    return val + (w - 1)[:, None] * qA - (np.abs(w) - 1)[:, None] * errA[:, None]


def _single_rows(val, g, H, m3, Jc, Sb, r, dim, gamma):
    """One constraint per box. With d = y - y_c (|d| <= r), the quadratic model of Psi~ reads
        g.ydot_c + d.H ydot_c + g.(W d) + d.H W d + gamma q~(x),
    and ydot_c = V + Omega x (c - p_A), W = [Omega]_x (relative twist of A with respect to B), so on the box
        Psi~ >= g.ydot_c(u) + gamma min_k beta_k(q~) - sigma_V |V|_1 - sigma_O |Omega|_1   with
        sigma_V = r |H| + M3 r^2 / 2,   sigma_O = (rho - r) r |H| + r |g| + r^2 |H| + rho M3 r^2 / 2,
    rho = |c - p_A| + r (the last column of Sb); |H| is the Frobenius norm (>= spectral)."""
    Hn = np.linalg.norm(H, axis=(1, 2))
    gn = np.linalg.norm(g, axis=1)
    rho = Sb[:, -1]
    sV = r * Hn + .5 * m3 * r ** 2
    sO = (rho - r) * r * Hn + r * gn + r ** 2 * Hn + rho * .5 * m3 * r ** 2
    k = Sb.shape[1]
    T = -np.c_[np.repeat(sV[:, None], dim, axis=1), np.repeat(sO[:, None], k - dim, axis=1)]
    A = np.einsum('nd,nmd->nm', g, Jc)
    return A, T, gamma * val.min(axis=1)


def _vertex_rows(val, g, H, m3, Jc, W, Sb, s, r, R, dim, gamma):
    """2^n constraints per box, one per box vertex. The part of the quadratic model of Psi~ that is linear in
    d = y - y_c, (g + H d).ydot_c + g.(W d) = g.ydot_c + d.(H ydot_c) + (W^T g).d, attains its minimum over the
    (rotated) box at a vertex d_k; the value part is bounded by its minimum Bernstein coefficient (constant),
    and the remaining terms by
        |d.H W d| <= r^2 |H| |Omega|,   |grad phi_B - grad q_B| |ydot| <= M3 r^2 / 2 (|V| + rho |Omega|).
    Rows: g.J u + d_k.(H J u) + g.(W(u) d_k) + gamma min beta - sigma_V |V|_1 - sigma_O |Omega|_1 >= 0."""
    off = np.array(list(itertools.product((-.5, .5), repeat=dim))) * s @ R.T        # (2^n, dim) in B's frame
    Hd = np.einsum('nab,kb->nka', H, off)                                         # (n, 2^n, dim)
    A = (np.einsum('nd,nmd->nm', g, Jc)[:, None, :] + np.einsum('nka,nma->nkm', Hd, Jc)
         + np.einsum('nd,mdb,kb->nkm', g, W, off))                                # (n, 2^n, m)
    Hn = np.linalg.norm(H, axis=(1, 2))
    rho = Sb[:, -1]
    sV = .5 * m3 * r ** 2
    sO = rho * .5 * m3 * r ** 2 + r ** 2 * Hn
    k = Sb.shape[1]
    T = -np.c_[np.repeat(sV[:, None], dim, axis=1), np.repeat(sO[:, None], k - dim, axis=1)]
    n, K, m = A.shape
    return A.reshape(-1, m), np.repeat(T, K, axis=0), np.repeat(gamma * val.min(axis=1), K)


def _children(D, bodyA, s, r, R, W, dim):
    """The 2^dim half-size boxes of every box in D, without those that contain no point of S_A."""
    TK, QT = TABLES[dim]
    off = np.array(list(itertools.product((-1., 1.), repeat=dim))) * s / 4          # (2^dim, dim), A frame
    n, c = len(off), D['c']
    cc = (c[:, None] + off[None]).reshape(-1, dim)
    rc, sc = r / 2, s / 2
    v, g, H, m3 = _quad_model(bodyA['f'], bodyA['level'], cc, dim, rc, bodyA['cell3'])
    cA = v[:, None] + sc * g @ TK.T + .5 * sc ** 2 * np.einsum('kab,nab->nk', QT, H)
    has = (cA.min(axis=1) - m3 * rc ** 3 / 6 <= 0) & (cA.max(axis=1) + m3 * rc ** 3 / 6 >= 0)
    dB = off @ R.T                                                                   # offsets in B's frame
    rep = lambda X: np.repeat(X, n, axis=0)
    yc = (D['yc'][:, None] + dB[None]).reshape(-1, dim)
    Jc = rep(D['Jc']) + np.tile(np.einsum('mab,kb->kma', W, dB), (len(c), 1, 1))
    C = dict(c=cc, yc=yc, vA=v, gA=g, HA=H, m3A=m3, Jc=Jc, Sb=rep(D['Sb']))           # rho of the parent still bounds
    return {k: x[has] for k, x in C.items()}


def sampled_tightening(yc, g, rbar, r, sd):
    """Per-box coefficients of the sampled-data tightening (Theorem sampled of the paper) on the auxiliary
    variables (a_V, a_Omega, s): with the input held for dt, every box point satisfies
        |ydot| <= sigma = 1.a_V + rbar 1.a_Omega + dt c2 nu 1.s <= sigma_bar,   |yddot| <= c2 nu 1.s,
    where c2 = 1.5 D (D bounds the lever arms of the link), nu >= |u|_1 and the twist caps are enforced in the
    QP, so |d^2/dt^2 phi_B(y)| <= H sigma_bar sigma + G c2 nu 1.s with H, G bounding |Hess phi_B|, |grad phi_B|
    on the swept ball of radius r + travel. The row is tightened by dt/2 times this bound."""
    dt, V, Om, nu, c2, travel, m = sd['dt'], sd['V'], sd['Om'], sd['nu'], sd['c2'], sd['travel'], sd['m']
    P3 = yc if yc.shape[1] == 3 else np.c_[yc, np.zeros(len(yc))]
    Hq = sd['M2'].eval(P3)
    Gq = np.linalg.norm(g, axis=1) + Hq * (r + travel)
    sbar = V + rbar * Om + dt * c2 * nu * nu
    k = Hq * sbar
    dim = yc.shape[1]
    return .5 * dt * np.c_[np.repeat(k[:, None], dim, axis=1), np.repeat((k * rbar)[:, None], dim, axis=1),
                           np.repeat((k * dt * c2 * nu + Gq * c2 * nu)[:, None], m, axis=1)]


def rows(bodyA, idx, RbA, yc, fB, lB, M3B, Jc, W, Sb, eta, umax, auxmax, gamma=GAMMA_DEFAULT, prune_rows=True,
         theta=None, depth=None, near=None, row_mode=None, mult_mode=None, sd=None):
    """Coefficient rows for the boxes idx of A (centers yc in B's frame, rotation RbA: A frame -> B frame).
    Jc: (n, m, dim) center velocity per input (B frame); W: (m, dim, dim) rotation generators (B frame).
    Sb: (n, k) speed coefficients, |ydot(x)| <= Sb . a on the box for auxiliary variables a >= |E u|
    (a relative twist of the two bodies, componentwise); umax (m,), auxmax (k,): bounds of |u| and a over the
    admissible inputs (used only to drop rows that no admissible input can violate).
    An active box whose gradient remainder M3_B r^2 / 2 exceeds theta is replaced by its half-size boxes that
    meet S_A (at most `depth` times); they cover the same part of the surface.
    Returns A (N, m), T (N, k), C (N,) with rows A u + T a + C >= 0, the number of boxes with rows, the
    minimum over those boxes of the certified lower bound of h~, and the indices (into idx) of the active boxes.
    sd (optional): sampled-data data (see sampled_tightening); the auxiliary vector is then (a, s) with s >= |u|."""
    theta = REFINE_THETA if theta is None else theta
    depth = REFINE_DEPTH if depth is None else depth
    near = REFINE_NEAR if near is None else near
    row_mode = ROW_MODE if row_mode is None else row_mode
    mult_mode = MULT_MODE if mult_mode is None else mult_mode
    if row_mode not in ROW_MODES or mult_mode not in MULT_MODES:
        raise ValueError(f'Unknown certificate modes: {row_mode}, {mult_mode}')
    dim, s, r = bodyA['dim'], bodyA['side'], bodyA['r']
    R = RbA[:dim, :dim]
    D = dict(c=bodyA['PC'][idx], yc=yc, vA=bodyA['vA'][idx], gA=bodyA['gA'][idx], HA=bodyA['HA'][idx],
             m3A=bodyA['m3A'][idx], Jc=Jc, Sb=Sb)
    out, n_box, h_low, first = [], 0, np.inf, None
    for level_ in range(depth + 1):
        if not len(D['c']):
            break
        val, g, H, m3 = _layer(D, s, r, R, fB, lB, M3B, dim, gamma)
        active = val.min(axis=1) < eta
        if first is None:
            first = np.flatnonzero(active)                # indices (into idx) of the active boxes
        refine = (active & (.5 * m3 * r ** 2 > theta) & (val.min(axis=1) < near) if level_ < depth
                  else np.zeros_like(active))
        final = active & ~refine
        if final.any() and mult_mode == 'free':
            # joint certificate with a per-box multiplier w (a QP variable): the rows
            #   b_k(velocity part)(u) + gamma w beta_k^A + gamma beta_k^B - remainders(a) >= 0
            # are affine in (u, a, w); beta^A carries -M3_A r^3/6, so the remainder scales with w >= 0.
            valA, valB = split_values({k: x[final] for k, x in D.items()}, val[final], s, r, dim)
            Ar_, Tr_, Cr_ = _velocity_rows(valB, g[final], H[final], m3[final], D['Jc'][final], W,
                                           D['Sb'][final], s, r, R, dim, gamma)
            K = valA.shape[1]
            if sd is not None:
                Tr_ = np.c_[Tr_, np.zeros((len(Tr_), sd['m']))] - np.repeat(
                    sampled_tightening(D['yc'][final], g[final], D['Sb'][final][:, -1], r, sd), K, axis=0)
            out.append((Ar_, Tr_, Cr_, gamma * valA.ravel(), np.repeat(n_box + np.arange(len(valA)), K)))
            n_box += int(final.sum())
            h_low = min(h_low, float(val[final].min()))
        elif final.any():
            cert_val = weighted_values(val[final], {k: x[final] for k, x in D.items()},
                                       g[final], R, s, r, dim, mult_mode)
            if row_mode == 'vertex':
                out.append(_vertex_rows(cert_val, g[final], H[final], m3[final], D['Jc'][final], W, D['Sb'][final],
                                        s, r, R, dim, gamma))
            elif row_mode == 'single':
                out.append(_single_rows(cert_val, g[final], H[final], m3[final], D['Jc'][final], D['Sb'][final],
                                        r, dim, gamma))
            else:
                out.append(_velocity_rows(cert_val, g[final], H[final], m3[final], D['Jc'][final], W,
                                          D['Sb'][final], s, r, R, dim, gamma))
            if sd is not None:
                Ar_, Tr_, Cr_ = out[-1]
                per = len(Cr_) // int(final.sum())
                Tr_ = np.c_[Tr_, np.zeros((len(Tr_), sd['m']))] - np.repeat(
                    sampled_tightening(D['yc'][final], g[final], D['Sb'][final][:, -1], r, sd), per, axis=0)
                out[-1] = (Ar_, Tr_, Cr_)
            n_box += int(final.sum())
            h_low = min(h_low, float(val[final].min()))
        if not refine.any():
            break
        D = _children({k: x[refine] for k, x in D.items()}, bodyA, s, r, R, W, dim)
        s, r = s / 2, r / 2
    if not out:
        return None
    Acoef, Tcoef, C = (np.vstack([o[0] for o in out]), np.vstack([o[1] for o in out]),
                       np.concatenate([o[2] for o in out]))
    if mult_mode == 'free':
        wcol, bid = np.concatenate([o[3] for o in out]), np.concatenate([o[4] for o in out])
        worst = C + np.minimum(0., wcol * W_MAX)           # the row's minimum over w in [0, W_MAX]
    else:
        worst = C
    if prune_rows:                                        # rows that some admissible input can violate
        need = worst < np.abs(Acoef) @ umax - Tcoef @ auxmax
    else:
        need = np.ones(len(C), bool)
    if mult_mode == 'free':                          # one nonzero per row: its box's multiplier
        nr = int(need.sum())
        Wc = sp.csr_matrix((wcol[need], (np.arange(nr), bid[need])), shape=(nr, n_box))
        return Acoef[need], Tcoef[need], C[need], n_box, h_low, first, Wc
    return Acoef[need], Tcoef[need], C[need], n_box, h_low, first


QP_SOLVER = os.environ.get('SUMMED_QP', 'nnls')  # 'nnls': dense least-distance program; 'clarabel': sparse IPM
QP_MARGIN = 1e-7                  # rows and input bounds are tightened by this much (conservative only)
QP_SPARSE_MIN_ROWS = int(os.environ.get('SUMMED_QP_MIN_ROWS', 2000))   # smaller QPs stay with the dense solver
QP_STATS = dict(sparse=0, fallback=0, infeasible=0, rejected=0, rescued=0)


def _dense(X):
    return X.toarray() if sp.issparse(X) else X


def violation(u, z, A, T, C, E, G_in, h_in, Wc=None):
    """Exact check of an input to be applied: with a = |E u| (the smallest admissible epigraph values) and the
    multipliers of z clipped to [0, W_MAX], returns max(0, -min row, -min input-bound row) in double precision."""
    m, k = len(u), len(E)
    a = np.abs(E @ u)
    r = A @ u + T @ a + C
    if Wc is not None and Wc.shape[1]:
        w = np.clip(1 + z[m + k:m + k + Wc.shape[1]] / np.sqrt(EPS_W), 0., W_MAX)
        r = r + Wc @ w
    v = max(0., -float(np.min(r))) if len(r) else 0.
    if len(G_in):
        v = max(v, -float(np.min(G_in @ u - h_in)))
    return v


def solve(u_nom, A, T, C, E, G_in, h_in, z_prev=None, Wc=None):
    """min |u - u_nom|^2 + EPS_T |a|^2 (+ EPS_W |w - 1|^2)  s.t.  A u + T a (+ Wc w) + C >= 0,  a >= |E u|,
    G_in u >= h_in (and 0 <= w <= W_MAX for 'free' multipliers, one per column of Wc).
    The feasible set in u equals that of A u + T |E u| + C >= 0 (T <= 0). Returns u, slack, z (warm start).
    No solver is trusted: the returned slack is the exact violation of the returned input (function violation),
    zero only if every row and input bound holds in double precision. A rejected solution of the dense
    least-distance solver (its NNLS can return an infeasible point without notice) is retried with the sparse
    interior-point solver. With QP_SOLVER = 'clarabel', large QPs go to the sparse solvers first."""
    tried_sparse = False
    if QP_SOLVER == 'clarabel' and len(C) >= QP_SPARSE_MIN_ROWS:
        res = _solve_working_set(u_nom, A, T, C, E, G_in, h_in, Wc, z_prev)   # warm start first
        if res is None:
            res, tried_sparse = _solve_sparse(u_nom, A, T, C, E, G_in, h_in, Wc), True
        if res is not None and violation(res[0], res[2], A, T, C, E, G_in, h_in, Wc) == 0.:
            return res
    u, s, z = _solve_dense(u_nom, A, T, C, E, G_in, h_in, z_prev, Wc)
    v = violation(u, z, A, T, C, E, G_in, h_in, Wc)
    if v == 0.:
        return u, 0., z
    if s == 0.:
        QP_STATS['rejected'] += 1
    if not tried_sparse:
        res = _solve_sparse(u_nom, A, T, C, E, G_in, h_in, Wc)
        if res is not None and violation(res[0], res[2], A, T, C, E, G_in, h_in, Wc) == 0.:
            QP_STATS['rescued'] += 1
            return res
    return u, max(v, s), z


def _solve_dense(u_nom, A, T, C, E, G_in, h_in, z_prev=None, Wc=None):
    T, Wc = _dense(T), (None if Wc is None else _dense(Wc))
    from sdf_cbf_utils import solve_ldp_qp
    m, k = len(u_nom), len(E)
    nw = 0 if Wc is None else Wc.shape[1]
    se, sw = np.sqrt(EPS_T), np.sqrt(EPS_W)
    I = np.eye(k)
    pad = lambda X, c: np.c_[X, np.zeros((len(X), c))]
    G = np.vstack([np.c_[A, T / se] if nw == 0 else np.c_[A, T / se, Wc / sw],
                   pad(np.c_[-E, I / se], nw), pad(np.c_[E, I / se], nw), pad(G_in, k + nw)])
    h = np.r_[(-C if nw == 0 else -C - Wc.sum(axis=1)) + QP_MARGIN, np.zeros(2 * k), h_in + QP_MARGIN]
    if nw:                                  # scaled variable sw (w - 1): w >= 0 and w <= W_MAX
        Iw = np.eye(nw)
        G = np.vstack([G, np.c_[np.zeros((nw, m + k)), Iw], np.c_[np.zeros((nw, m + k)), -Iw]])
        h = np.r_[h, np.full(nw, -sw), np.full(nw, -sw * (W_MAX - 1))]
    if z_prev is not None and len(z_prev) == m + k + nw:
        zp = z_prev
    elif nw and z_prev is not None and len(z_prev) >= m:   # the numbers of pairs and boxes changed
        zp = np.r_[z_prev[:m], np.zeros(k + nw)]
    else:
        zp = None
    z, s = solve_ldp_qp(np.r_[u_nom, np.zeros(k + nw)], G, h, len(C), zp)
    return z[:m], s, z


def _scaled_sparse(u_nom, A, T, C, E, G_in, h_in, Wc):
    """Sparse G z >= h of the dense path, in its scaled variables z = (u, sqrt(EPS_T) a, sqrt(EPS_W) (w - 1))."""
    m, k = len(u_nom), len(E)
    nw = 0 if Wc is None else Wc.shape[1]
    se, sw = np.sqrt(EPS_T), np.sqrt(EPS_W)
    Z = lambda r, c: sp.csr_matrix((r, c))
    Ik = sp.identity(k, format='csr') / se
    rows_ = [sp.hstack([sp.csr_matrix(A), sp.csr_matrix(T) / se] + ([sp.csr_matrix(Wc) / sw] if nw else [])),
             sp.hstack([sp.csr_matrix(-E), Ik] + ([Z(k, nw)] if nw else [])),
             sp.hstack([sp.csr_matrix(E), Ik] + ([Z(k, nw)] if nw else [])),
             sp.hstack([sp.csr_matrix(G_in), Z(len(G_in), k + nw)])]
    h = [-C - (np.asarray(sp.csr_matrix(Wc).sum(axis=1)).ravel() if nw else 0.) + QP_MARGIN, np.zeros(2 * k),
         h_in + QP_MARGIN]
    if nw:
        Iw = sp.identity(nw, format='csr')
        rows_ += [sp.hstack([Z(nw, m + k), Iw]), sp.hstack([Z(nw, m + k), -Iw])]
        h += [np.full(nw, -sw), np.full(nw, -sw * (W_MAX - 1))]
    return sp.vstack(rows_, format='csr'), np.concatenate(h), m, k, nw


QP_WS = dict(size=200, rounds=10, max_rows=3000, tol=1e-9)
QP_STATS['working_set'] = 0


def _solve_working_set(u_nom, A, T, C, E, G_in, h_in, Wc, z_prev):
    """The dense path's warm-started constraint generation, on sparse rows and with a bounded working set.
    A relaxation's minimizer that satisfies all rows (same tolerance as the dense path) is the QP minimizer.
    Returns None when there is no warm start, the working set does not converge or grows too large, or a
    relaxation is infeasible; the caller then uses the interior-point or dense solver."""
    from sdf_cbf_utils import _ldp
    G, h, m, k, nw = _scaled_sparse(u_nom, A, T, C, E, G_in, h_in, Wc)
    n = m + k + nw
    if z_prev is None or len(z_prev) < m:
        return None
    zp = z_prev if len(z_prev) == n else np.r_[z_prev[:m], np.zeros(k + nw)]
    z0 = np.r_[u_nom, np.zeros(k + nw)]
    work = np.zeros(len(h), dtype=bool)
    work[np.argsort(G @ zp - h)[:QP_WS['size']]] = True
    for _ in range(QP_WS['rounds']):
        z = _ldp(z0, G[work].toarray(), h[work])
        if z is None:
            return None
        viol = G @ z - h < -QP_WS['tol'] * np.maximum(1.0, np.abs(h))
        if not viol.any():
            QP_STATS['working_set'] += 1
            return z[:m], 0., z
        work |= viol
        if work.sum() > QP_WS['max_rows']:
            return None
    return None


def _solve_sparse(u_nom, A, T, C, E, G_in, h_in, Wc):
    """Sparse interior-point solution of the same QP (Clarabel). Returns (u, 0, z) if the applied variables
    satisfy every row exactly, (zeros, inf, None) if the solver certifies primal infeasibility (the caller's
    dense path then handles slack), and None to defer to the dense solver."""
    import clarabel
    m, k = len(u_nom), len(E)
    nw = 0 if Wc is None else Wc.shape[1]
    n = m + k + nw
    P = sp.diags(np.r_[np.full(m, 2.), np.full(k, 2 * EPS_T), np.full(nw, 2 * EPS_W)]).tocsc()
    q = np.r_[-2 * u_nom, np.zeros(k), np.full(nw, -2 * EPS_W)]
    Tm = sp.csr_matrix(T)
    blocks = [[sp.csr_matrix(A), Tm] + ([sp.csr_matrix(Wc)] if nw else [])]
    Ik = sp.identity(k, format='csr')
    Z = lambda r, c: sp.csr_matrix((r, c))
    rows_ = [sp.hstack(blocks[0]),
             sp.hstack([sp.csr_matrix(-E), Ik] + ([Z(k, nw)] if nw else [])),
             sp.hstack([sp.csr_matrix(E), Ik] + ([Z(k, nw)] if nw else [])),
             sp.hstack([sp.csr_matrix(G_in), Z(len(G_in), k + nw)])]
    h = [-C + QP_MARGIN, np.zeros(2 * k), h_in + QP_MARGIN]
    if nw:
        Iw = sp.identity(nw, format='csr')
        rows_ += [sp.hstack([Z(nw, m + k), Iw]), sp.hstack([Z(nw, m + k), -Iw])]
        h += [np.zeros(nw), np.full(nw, -W_MAX)]
    G = sp.vstack(rows_).tocsc()
    b = -np.concatenate(h)
    st = clarabel.DefaultSettings()
    st.verbose = False
    sol = clarabel.DefaultSolver(P, q, -G, b, [clarabel.NonnegativeConeT(G.shape[0])], st).solve()
    status = str(sol.status)
    if 'Infeasible' in status and 'Almost' not in status:
        QP_STATS['infeasible'] += 1
        return None
    if 'Solved' not in status:
        QP_STATS['fallback'] += 1
        return None
    z = np.asarray(sol.x)
    u = z[:m]
    if len(G_in) and np.any(G_in @ u < h_in):
        QP_STATS['fallback'] += 1
        return None
    a = np.abs(E @ u)                                # the smallest admissible epigraph values
    w = np.clip(z[m + k:], 0., W_MAX) if nw else None
    r = A @ u + (Tm @ a) + C + ((sp.csr_matrix(Wc) @ w) if nw else 0.)
    if r.min() < 0.:                                  # exact acceptance test of the applied variables
        QP_STATS['fallback'] += 1
        return None
    QP_STATS['sparse'] += 1
    return u, 0., np.r_[u, np.sqrt(EPS_T) * a, (np.sqrt(EPS_W) * (w - 1)) if nw else np.zeros(0)]
