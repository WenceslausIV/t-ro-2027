"""
3D prototype: collision avoidance between two non-convex rigid links, each represented by a smooth
tensor cubic B-spline SDF, WITHOUT configuration space and WITHOUT closest points.

  ours : the boundary of link A is covered offline by patches (triangle clusters of its mesh). For a
         patch with center c, radius r and hull vertices v_i, a second-order Taylor bound of phi_B with
         a certified Hessian bound M and the convex-hull property give, for every point x of the patch,
             phi_B(x) >= phi_B(c) + grad phi_B(c)^T (v_i - c) - M r^2 / 2   (minimized over hull vertices),
         so every hull vertex yields one smooth barrier  h_i = phi_B(c) + g^T(v_i - c) - M r^2/2 - l_B.
         No optimization online; patches far from B are pruned with a first-order bound.
  kkt  : ICRA'26-style closest point: x* = argmin phi_B(x) s.t. phi_A(x) = 0 (SLSQP, warm-started),
         h = phi_B(x*) - l_B, derivative through the closest point (envelope theorem).

Link A (a J-shaped hook) moves with SE(3) single-integrator dynamics (twist in the world frame) from
above into a U-shaped bracket B while rotating; the nominal goal is in collision, so the barrier
decides where A stops. Ground-truth distance: exact SDFs of the analytic (rounded-box union) shapes.

    python prototype_3d/proto3d.py
"""
import json
import os
import sys
import time

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')   # torch (via sdf_cbf_utils) + scipy OpenMP runtimes

import numpy as np
from scipy import ndimage
from scipy.cluster.vq import kmeans2
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from sdf_cbf_utils import Q_B2B, solve_ldp_qp  # noqa: E402

# ---------------------------------------------------------------------------------------------
# analytic ground truth: unions of rounded boxes
# ---------------------------------------------------------------------------------------------
ROUND = .01
A_BOXES = [((0., 0., 0.), (.25, .035, .035)),          # main bar along x
           ((.25, .09, 0.), (.035, .125, .035)),        # hook post
           ((.17, .20, 0.), (.10, .03, .035))]          # hook tip returning over the bar
B_BOXES = [((0., 0., 0.), (.32, .04, .25)),             # base plate
           ((-.28, .30, 0.), (.04, .30, .25)),          # left arm
           ((.28, .30, 0.), (.04, .30, .25))]           # right arm


def sdf_boxes(P, boxes):
    d = np.full(len(P), np.inf)
    for c, b in boxes:
        q = np.abs(P - np.array(c)) - (np.array(b) - ROUND)
        d = np.minimum(d, np.linalg.norm(np.maximum(q, 0), axis=1) + np.minimum(q.max(axis=1), 0) - ROUND)
    return d


def rotm(axis, ang):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def logm_so3(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    th = np.arccos(c)
    if th < 1e-9:
        return np.zeros(3)
    return th / (2 * np.sin(th)) * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])


# ---------------------------------------------------------------------------------------------
# 3D tensor cubic B-spline field with value / gradient / Hessian and a certified Hessian bound
# ---------------------------------------------------------------------------------------------
def _b(t):
    return np.stack([(1 - t) ** 3, 3 * t ** 3 - 6 * t ** 2 + 4, -3 * t ** 3 + 3 * t ** 2 + 3 * t + 1, t ** 3], -1) / 6


def _db(t):
    return np.stack([-(1 - t) ** 2, 3 * t ** 2 - 4 * t, -3 * t ** 2 + 2 * t + 1, t ** 2], -1) / 2


def _ddb(t):
    return np.stack([1 - t, 3 * t - 2, -3 * t + 1, t], -1)


class Spline3:
    def __init__(self, lo, hi, h):
        self.lo, self.h = np.asarray(lo, float), h
        self.K = np.ceil((np.asarray(hi, float) - self.lo) / h).astype(int)
        self.n = self.K + 3

    def _design(self, ax, x):
        s = (x - self.lo[ax]) / self.h
        c = np.clip(np.floor(s).astype(int), 0, self.K[ax] - 1)
        B = _b(np.clip(s - c, 0, 1))
        D = np.zeros((len(x), self.n[ax]))
        np.add.at(D, (np.repeat(np.arange(len(x)), 4), (c[:, None] + np.arange(4)).ravel()), B.ravel())
        return D

    def fit(self, fun, spacing, lam=1e-4):
        axes = [np.arange(self.lo[a], self.lo[a] + self.K[a] * self.h + 1e-9, spacing) for a in range(3)]
        X, Y, Z = np.meshgrid(*axes, indexing='ij')
        T = fun(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)
        maps = []
        for a in range(3):
            D = self._design(a, axes[a])
            maps.append(np.linalg.solve(D.T @ D + lam * np.eye(self.n[a]), D.T))
        self.W = np.einsum('ai,bj,ck,ijk->abc', maps[0], maps[1], maps[2], T, optimize=True)
        return self

    def hessian_bound(self, cellwise=False):
        """Certified bound on the Frobenius norm of the Hessian from the Bernstein coefficients of every
        second derivative: over the whole domain, or (cellwise=True) one bound per cell."""
        W = self.W
        idx = [np.arange(self.K[a])[:, None] + np.arange(4) for a in range(3)]
        Wg = W[idx[0][:, None, None, :, None, None], idx[1][None, :, None, None, :, None],
               idx[2][None, None, :, None, None, :]]
        C = np.einsum('ip,jq,kr,abcpqr->abcijk', Q_B2B, Q_B2B, Q_B2B, Wg, optimize=True)
        h2 = self.h ** 2
        red = (lambda x: np.abs(x).max(axis=(3, 4, 5))) if cellwise else (lambda x: np.abs(x).max())
        m = [[None] * 3 for _ in range(3)]
        for a in range(3):
            m[a][a] = red(6 * np.diff(C, n=2, axis=3 + a)) / h2
        for a in range(3):
            for b in range(a + 1, 3):
                m[a][b] = m[b][a] = red(9 * np.diff(np.diff(C, axis=3 + a), axis=3 + b)) / h2
        fro = np.sqrt(sum(m[a][b] ** 2 for a in range(3) for b in range(3)))
        return fro if cellwise else float(fro)

    def cell_coeffs(self):
        idx = [np.arange(self.K[a])[:, None] + np.arange(4) for a in range(3)]
        Wg = self.W[idx[0][:, None, None, :, None, None], idx[1][None, :, None, None, :, None],
                    idx[2][None, None, :, None, None, :]]
        return np.einsum('ip,jq,kr,abcpqr->abcijk', Q_B2B, Q_B2B, Q_B2B, Wg, optimize=True)

    def cell_of(self, P):
        return np.clip(np.floor((np.atleast_2d(P) - self.lo) / self.h).astype(int), 0, self.K - 1)

    def eval(self, P, order=0):
        P = np.atleast_2d(P)
        s = (P - self.lo) / self.h
        c = np.clip(np.floor(s).astype(int), 0, self.K - 1)
        t = np.clip(s - c, 0, 1)
        n1, n2 = self.W.shape[1], self.W.shape[2]
        if not hasattr(self, '_off') or self._off_shape != self.W.shape:
            g = np.arange(4)
            self._off = (g[:, None, None] * n1 * n2 + g[None, :, None] * n2 + g[None, None, :]).ravel()
            self._off_shape = self.W.shape
        base = (c[:, 0] * n1 + c[:, 1]) * n2 + c[:, 2]
        Wl = self.W.ravel()[base[:, None] + self._off].reshape(-1, 4, 4, 4)
        bas = [[_b(t[:, a])] for a in range(3)]                     # per axis: value, 1st, 2nd derivative
        if order >= 1:
            for a in range(3):
                bas[a].append(_db(t[:, a]) / self.h)
        if order >= 2:
            for a in range(3):
                bas[a].append(_ddb(t[:, a]) / self.h ** 2)
        Az = [np.einsum('nijk,nk->nij', Wl, bz) for bz in bas[2]]
        f = lambda mx, my, mz: np.einsum('ni,ni->n', np.einsum('nij,nj->ni', Az[mz], bas[1][my]), bas[0][mx])
        val = f(0, 0, 0)
        if order == 0:
            return val
        grad = np.stack([f(1, 0, 0), f(0, 1, 0), f(0, 0, 1)], 1)
        if order == 1:
            return val, grad
        H = np.empty((len(P), 3, 3))
        H[:, 0, 0] = f(2, 0, 0); H[:, 1, 1] = f(0, 2, 0); H[:, 2, 2] = f(0, 0, 2)
        H[:, 0, 1] = H[:, 1, 0] = f(1, 1, 0)
        H[:, 0, 2] = H[:, 2, 0] = f(1, 0, 1)
        H[:, 1, 2] = H[:, 2, 1] = f(0, 1, 1)
        return val, grad, H


# ---------------------------------------------------------------------------------------------
# offline cover of link A's boundary
# ---------------------------------------------------------------------------------------------
def mesh(boxes, lo, hi, res):
    axes = [np.arange(lo[a], hi[a] + 1e-9, res) for a in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing='ij')
    V = sdf_boxes(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1), boxes).reshape(X.shape)
    verts, faces, _, _ = measure.marching_cubes(V, 0., spacing=(res, res, res))
    return verts + np.asarray(lo), faces


def cover(verts, faces, n_patches, seed=0):
    cent = verts[faces].mean(axis=1)
    _, lab = kmeans2(cent, n_patches, seed=seed, minit='++')
    patches = []
    for k in range(n_patches):
        vid = np.unique(faces[lab == k])
        if len(vid) == 0:
            continue
        pts = verts[vid]
        c = pts.mean(axis=0)
        try:
            hv = pts[ConvexHull(pts).vertices]
        except Exception:
            hv = pts
        patches.append(dict(c=c, r=float(np.linalg.norm(pts - c, axis=1).max()), hull=hv))
    return patches


def surface_boxes(f, side=None, level=0.):
    """Certified cover of the continuous zero level set {f = 0}: a box can contain the level set only
    if its restricted Bernstein coefficients change sign (convex-hull property); such boxes are halved
    in all three axes until their side is <= side. With side=None, keep native patches
    without subdivision. Returns lower corners (n,3) and the side."""
    from cspace_cbf_5robots import sub_matrix_np
    side = f.h if side is None else side
    C = f.cell_coeffs() - level                    # Bernstein coefficients of f - level
    cells = np.argwhere((C.min(axis=(3, 4, 5)) <= 0) & (C.max(axis=(3, 4, 5)) >= 0))
    co = C[cells[:, 0], cells[:, 1], cells[:, 2]]
    lo = f.lo + cells * f.h
    size = f.h
    S = [sub_matrix_np(np.array([0.]), np.array([.5]))[0], sub_matrix_np(np.array([.5]), np.array([1.]))[0]]
    while size > side:
        los, cos = [], []
        for bx in (0, 1):
            for by in (0, 1):
                for bz in (0, 1):
                    cc = np.einsum('ip,jq,kr,npqr->nijk', S[bx], S[by], S[bz], co)
                    keep = (cc.min(axis=(1, 2, 3)) <= 0) & (cc.max(axis=(1, 2, 3)) >= 0)
                    los.append(lo[keep] + size / 2 * np.array([bx, by, bz])); cos.append(cc[keep])
        lo, co, size = np.concatenate(los), np.concatenate(cos), size / 2
    return lo, size


def box_patches(lo, size):
    corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float) * size
    return [dict(c=l + size / 2, r=size * np.sqrt(3) / 2, hull=l + corners) for l in lo]


# ---------------------------------------------------------------------------------------------
# closed loop
# ---------------------------------------------------------------------------------------------
def box_rows(vmax, wmax):
    G, h = [], []
    for a in range(6):
        e = np.zeros(6); e[a] = 1.
        lim = vmax if a < 3 else wmax
        G += [e, -e]; h += [-lim, -lim]
    return np.array(G), np.array(h)


def run(method, fB, fA, patches, Mcell, lB, gtA_pts, start, goal, steps=1000, dt=.01, gamma=5., vmax=.5, wmax=1.,
        lsA_pts=None,
        act=.05):
    p, R = np.array(start[0], float), np.array(start[1], float)
    pg, Rg = np.array(goal[0], float), np.array(goal[1], float)
    Gb, hb = box_rows(vmax, wmax)
    log = dict(t_ctrl=[], h=[], d=[], dsdf=[], u=[], rows=[], p=[], R=[])
    xstar = None
    PC = np.array([pt['c'] for pt in patches]); PR = np.array([pt['r'] for pt in patches])
    HV = np.concatenate([pt['hull'] for pt in patches])
    PID = np.concatenate([np.full(len(pt['hull']), i) for i, pt in enumerate(patches)])
    for k in range(steps):
        t0 = time.perf_counter()
        V = np.clip(1.5 * (pg - p), -vmax, vmax)
        W = np.clip(2. * (R @ logm_so3(R.T @ Rg)), -wmax, wmax)
        u_nom = np.r_[V, W]
        rows, hs = [], []
        if method == 'ours':
            C = PC @ R.T + p                                    # patch centers (world)
            val, g = fB.eval(C, order=1)
            ci = fB.cell_of(C)
            Mloc = Mcell[ci[:, 0], ci[:, 1], ci[:, 2]]           # certified local Hessian bound
            marg = .5 * Mloc * PR ** 2
            near = np.flatnonzero(val - np.linalg.norm(g, axis=1) * PR - marg - lB < act)
            if len(near):
                v2, g2, H2 = fB.eval(C[near], order=2)
                pos = np.full(len(PC), -1); pos[near] = np.arange(len(near))
                sel = np.flatnonzero(pos[PID] >= 0)              # hull vertices of the near patches
                q = pos[PID[sel]]
                hv = HV[sel] @ R.T + p
                dv = hv - C[near][q]
                hh = v2[q] + np.einsum('nd,nd->n', dv, g2[q]) - marg[near][q] - lB
                keep = hh < act
                dv, hv, hh, q = dv[keep], hv[keep], hh[keep], q[keep]
                Hd = np.einsum('nd,nde->ne', dv, H2[q])
                rowV = Hd + g2[q]
                rowW = np.cross(C[near][q] - p, Hd) + np.cross(hv - p, g2[q])
                rows, hs = np.hstack([rowV, rowW]), hh
            else:
                rows, hs = np.zeros((0, 6)), np.zeros(0)
        else:
            if xstar is None:                         # initialize at A's surface point with smallest phi_B
                pts = gtA_pts @ R.T + p
                xstar = pts[np.argmin(fB.eval(pts))]
            loc = lambda x: R.T @ (x - p)
            cons = dict(type='eq', fun=lambda x: fA.eval(loc(x)[None])[0],
                        jac=lambda x: fA.eval(loc(x)[None], order=1)[1][0] @ R.T)
            res = minimize(lambda x: fB.eval(x[None])[0], xstar,
                           jac=lambda x: fB.eval(x[None], order=1)[1][0], constraints=[cons],
                           method='SLSQP', options=dict(maxiter=50, ftol=1e-9))
            xstar = res.x
            val, g = fB.eval(xstar[None], order=1)
            rr = xstar - p
            rows = np.hstack([g[0], np.cross(rr, g[0])])[None]
            hs = np.array([val[0] - lB])
        G = np.vstack([rows, Gb]) if len(hs) else Gb
        hh = np.r_[-gamma * hs, hb]
        u, _ = solve_ldp_qp(u_nom, G, hh, len(hs), None)
        log['t_ctrl'].append(time.perf_counter() - t0)
        log['rows'].append(len(hs))
        log['h'].append(float(hs.min()) if len(hs) else np.nan)
        pts = gtA_pts @ R.T + p                        # ground truth: A's dense surface vs exact SDF of B
        log['d'].append(float(sdf_boxes(pts, B_BOXES).min()))
        if lsA_pts is not None:                        # SDF-defined bodies: phi_B on A's level set
            log['dsdf'].append(float(fB.eval(lsA_pts @ R.T + p).min()))
        log['u'].append(u.copy()); log['p'].append(p.copy()); log['R'].append(R.copy())
        p = p + dt * u[:3]
        R = rotm(u[3:], dt * np.linalg.norm(u[3:])) @ R if np.linalg.norm(u[3:]) > 1e-12 else R
    return log


def main():
    t0 = time.perf_counter()
    fB = Spline3((-.7, -.4, -.6), (.7, 1.4, .6), .03).fit(lambda P: sdf_boxes(P, B_BOXES), .02)
    fA = Spline3((-.35, -.15, -.15), (.4, .35, .15), .02).fit(lambda P: sdf_boxes(P, A_BOXES), .01)
    M = fB.hessian_bound()
    from scipy.ndimage import maximum_filter
    Mcell = maximum_filter(fB.hessian_bound(cellwise=True), size=5)   # covers patch balls up to 2 cells
    lB = 0.                                                # bodies are {phi_A <= 0} and {phi_B <= 0}
    blo, bside = surface_boxes(fA, .01)                    # certified cover of {phi_A = 0}
    patches = box_patches(blo, bside)
    ax_ = [np.arange(fA.lo[a], fA.lo[a] + fA.K[a] * fA.h, .004) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax_, indexing='ij')
    VA = fA.eval(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)
    lsA = measure.marching_cubes(VA, 0., spacing=(.004,) * 3)[0] + fA.lo    # samples of {phi_A = 0}
    gtA, _ = mesh(A_BOXES, (-.3, -.1, -.1), (.35, .3, .1), .004)
    setup = dict(t_setup=time.perf_counter() - t0, M=M, M_cell_median=float(np.median(Mcell)), lB=lB,
                 n_boxes=len(patches), box_side=bside,
                 hull_vertices=int(sum(len(pt['hull']) for pt in patches)),
                 r_max=max(pt['r'] for pt in patches), coef_B=int(np.prod(fB.n)))
    print('setup', setup, flush=True)
    start = ((0., 1.0, 0.), rotm((1, 0, 0), .35))
    goal = ((0., .22, 0.), rotm((0, 0, 1), np.pi / 2) @ rotm((1, 0, 0), .2))
    out = dict(setup=setup)
    for method in ('ours', 'kkt'):
        L = run(method, fB, fA, patches, Mcell, lB, gtA, start, goal, lsA_pts=lsA)
        u = np.array(L['u'])
        out[method] = dict(t_med_ms=1e3 * float(np.median(L['t_ctrl'])), t_max_ms=1e3 * float(np.max(L['t_ctrl'])),
                           min_true_dist_mm=1e3 * float(np.min(L['d'])), final_true_dist_mm=1e3 * L['d'][-1],
                           min_sdf_gap_mm=1e3 * float(np.min(L['dsdf'])),
                           min_h_mm=1e3 * float(np.nanmin(L['h'])), tv=float(np.abs(np.diff(u, axis=0)).sum()),
                           rows_med=float(np.median(L['rows'])), rows_max=int(np.max(L['rows'])),
                           final_p=L['p'][-1].tolist())
        np.savez_compressed(os.path.join(HERE, f'traj_{method}.npz'), p=np.array(L['p']), R=np.array(L['R']),
                            d=np.array(L['d']), h=np.array(L['h']), u=u, t=np.array(L['t_ctrl']))
        print(method, out[method], flush=True)
    json.dump(out, open(os.path.join(HERE, 'results.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
