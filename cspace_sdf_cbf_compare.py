"""
Position switching of two robots with smooth, NON-CONVEX boundaries, WITH rotation, using a
C-space SDF CBF with a certified enclosure.  Compares two representations of the same 3-D SDF:

    bspline : tensor-product cubic B-spline (C^2), theta periodic
    bezier  : tensor-product piecewise cubic Bezier (C^0, one Bernstein patch per cell)
              with the same number of coefficients (~48^3)

Relative configuration of A w.r.t. B:  q = (t, theta),  t = R(-theta_B)(p_A - p_B),
theta = theta_A - theta_B.  C-obstacle CO = {q : (R(theta)A + t) meets B}.

Fit (offline): phi(q) ~ signed distance to CO, from rasterised C-space slices (training data
only; the guarantee does not rely on it).

Certified enclosure (offline): every configuration where the boundaries touch satisfies
    t = b(tau) - R(theta) a(sigma)      (dA (+) dB contains d(A (+) B)),
so dCO lies in M = {(b(tau) - R(theta)a(sigma), theta)}.  Branch-and-bound over (sigma, tau,
theta) with Bernstein upper bounds of phi on boxes containing each piece of M gives a certified
    l* >= max_M phi.
Then h(q) = phi(q) - l* >= 0 excludes dCO, and a trajectory starting outside CO can never
reach CO (it would have to cross dCO, where h < 0).  No Minkowski sum is ever formed.

Online: one scalar CBF constraint for the pair, joint QP over both robots' (vx, vy, omega).
"""
import argparse
import itertools
import math
import time

import numpy as np
import torch
from matplotlib.path import Path
from scipy import ndimage
from scipy.signal import fftconvolve

from sdf_cbf_utils import (Q_B2B, cubic_bspline_local, cubic_bspline_local_deriv,
                                             polygon_sdf, solve_ldp_qp)

DEV = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
TWO_PI = 2 * np.pi
ACT_DELTA = 0.125   # CBF of a pair is active for |t|_inf < D - ACT_DELTA (Assumption 1 band width) [m]

# cubic Bernstein -> power basis:  p(s) = [1, s, s^2, s^3] MB P
MB = np.array([[1, 0, 0, 0], [-3, 3, 0, 0], [3, -6, 3, 0], [-1, 3, -3, 1]], dtype=np.float64)
MB_INV = np.linalg.inv(MB)


def bern3(s):
    return np.stack([(1 - s) ** 3, 3 * s * (1 - s) ** 2, 3 * s ** 2 * (1 - s), s ** 3], axis=-1)


def bern3_d(s):
    return np.stack([-3 * (1 - s) ** 2, 3 * (1 - s) ** 2 - 6 * s * (1 - s),
                     6 * s * (1 - s) - 3 * s ** 2, 3 * s ** 2], axis=-1)


def rot(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s], [s, c]])


# ---------------------------------------------------------------------------------------------
# Robot shapes: closed periodic cubic B-spline curves (smooth, non-convex)
# ---------------------------------------------------------------------------------------------
class SmoothShape:
    def __init__(self, ctrl, n_dense=480):
        ctrl = np.asarray(ctrl, dtype=np.float64)
        self.ctrl = ctrl
        K = len(ctrl)
        idx = (np.arange(K)[:, None] + np.arange(4)[None]) % K
        self.bez = np.einsum('ip,kpd->kid', Q_B2B, ctrl[idx])            # (K, 4, 2) per span
        s = np.arange(n_dense // K) / (n_dense // K)
        self.dense = np.einsum('si,kid->ksd', bern3(s), self.bez).reshape(-1, 2)
        self.rho = np.linalg.norm(self.dense, axis=1).max()
        self.edges_t = None

    def world(self, pose):
        return self.dense @ rot(pose[2]).T + pose[:2]


def shape_A():
    """C-shaped (horseshoe) body: outer arc CCW, inner arc back."""
    out = [0.56 * np.array([np.cos(a), np.sin(a)]) for a in np.deg2rad(np.linspace(45, 315, 9))]
    inn = [0.26 * np.array([np.cos(a), np.sin(a)]) for a in np.deg2rad(np.linspace(300, 60, 6))]
    return SmoothShape(out + inn)


def shape_B():
    """Three-lobed flower."""
    a = np.linspace(0, TWO_PI, 15, endpoint=False)
    r = 0.44 + 0.2 * np.cos(3 * a)
    return SmoothShape(np.stack([r * np.cos(a), r * np.sin(a)], axis=1))


# ---------------------------------------------------------------------------------------------
# Training data: rasterised C-space slices (approximate; certification is separate)
# ---------------------------------------------------------------------------------------------
def cspace_sdf_slices(A, B, thetas, D, res=0.01):
    n = 2 * int(round(D / res)) + 1
    g = (np.arange(n) - n // 2) * res
    GX, GY = np.meshgrid(g, g, indexing='ij')
    P = np.stack([GX.ravel(), GY.ravel()], axis=1)
    mB = Path(B.dense).contains_points(P).reshape(n, n).astype(float)
    out = np.empty((len(thetas), n, n))
    for k, th in enumerate(thetas):
        mnegA = Path(-(A.dense @ rot(th).T)).contains_points(P).reshape(n, n).astype(float)
        co = fftconvolve(mB, mnegA, mode='same') > 0.5                   # (R(th)A + t) meets B
        out[k] = (ndimage.distance_transform_edt(~co) - ndimage.distance_transform_edt(co)) * res
    return g, out


def _raster_polygon(P, n, res):
    """
    Occupancy of polygon P at the grid points g = (i - n//2) res, as an (n, n) array indexed
    [ix, iy]; identical to Path.contains_points: PIL fills the interior, and only the one-pixel
    ring on which PIL's fill and outline disagree is decided by the exact point-in-polygon test.
    """
    from PIL import Image, ImageDraw
    Q = [tuple(p) for p in (P / res + n // 2)]
    masks = []
    for outline in (0, 1):
        img = Image.new('1', (n, n), 0)
        ImageDraw.Draw(img).polygon(Q, fill=1, outline=outline)
        masks.append(np.array(img, dtype=bool).T)
    ring = np.argwhere(masks[0] != masks[1])
    ring = np.unique(np.concatenate([ring + d for d in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1))]), axis=0)
    ring = ring[(ring >= 0).all(1) & (ring < n).all(1)]
    m = masks[0].copy()
    m[ring[:, 0], ring[:, 1]] = Path(P).contains_points((ring - n // 2) * res)
    return m.astype(np.float64)


def cspace_sdf_samples(A, B, thetas, D, res, pts_x, pts_y):
    """
    Same targets as cspace_sdf_slices + sample_slices, streamed slice by slice with the occupancy
    convolution on the GPU (FFT) and polygon masks rasterized with PIL; returns (nx, ny, n_theta).
    """
    n = 2 * int(round(D / res)) + 1
    g = (np.arange(n) - n // 2) * res
    X, Y = np.meshgrid((pts_x - g[0]) / res, (pts_y - g[0]) / res, indexing='ij')
    m = 2 * n - 1
    FB = torch.fft.rfft2(torch.as_tensor(_raster_polygon(B.dense, n, res), device=DEV), s=(m, m))
    out = np.empty((len(pts_x), len(pts_y), len(thetas)))
    for k, th in enumerate(thetas):
        mA = torch.as_tensor(_raster_polygon(-(A.dense @ rot(th).T), n, res), device=DEV)
        full = torch.fft.irfft2(FB * torch.fft.rfft2(mA, s=(m, m)), s=(m, m))
        co = (full[n // 2:n // 2 + n, n // 2:n // 2 + n] > 0.5).cpu().numpy()   # mode='same'
        s = (ndimage.distance_transform_edt(~co) - ndimage.distance_transform_edt(co)) * res
        out[:, :, k] = ndimage.map_coordinates(s, [X, Y], order=1, mode='nearest')
    return out


def sample_slices(g, slices, pts_x, pts_y):
    """Bilinear sample of every slice at the grid pts_x x pts_y -> (n_theta, nx, ny)."""
    res = g[1] - g[0]
    X, Y = np.meshgrid((pts_x - g[0]) / res, (pts_y - g[0]) / res, indexing='ij')
    return np.stack([ndimage.map_coordinates(s, [X, Y], order=1, mode='nearest') for s in slices])


# ---------------------------------------------------------------------------------------------
# 1-D bases and the 3-D field (both kinds share cell-wise Bernstein coefficients)
# ---------------------------------------------------------------------------------------------
class Basis1D:
    def __init__(self, kind, lo, hi, K, periodic):
        self.kind, self.lo, self.hi, self.K, self.periodic = kind, lo, hi, K, periodic
        self.h = (hi - lo) / K
        if kind == 'bspline':
            self.n = K if periodic else K + 3
            self.Q = Q_B2B
        else:
            self.n = 3 * K if periodic else 3 * K + 1
            self.Q = np.eye(4)

    def cell(self, x):
        t = (np.asarray(x, dtype=np.float64) - self.lo) / self.h
        c = np.floor(t).astype(int)
        if self.periodic:
            return c % self.K, t - np.floor(t)
        c = np.clip(c, 0, self.K - 1)
        return c, np.clip(t - c, 0.0, 1.0)

    def ctrl_idx(self, c):
        step = 1 if self.kind == 'bspline' else 3
        idx = step * np.asarray(c)[..., None] + np.arange(4)
        return idx % self.n if self.periodic else idx

    def design(self, x):
        c, xi = self.cell(x)
        B = cubic_bspline_local(xi) if self.kind == 'bspline' else bern3(xi)
        out = np.zeros((len(c), self.n))
        np.add.at(out, (np.repeat(np.arange(len(c)), 4), self.ctrl_idx(c).ravel()), B.ravel())
        return out


class Field3D:
    def __init__(self, kind, D, K_xy, K_th):
        self.kind = kind
        self.bx = Basis1D(kind, -D, D, K_xy, False)
        self.by = Basis1D(kind, -D, D, K_xy, False)
        self.bt = Basis1D(kind, 0.0, TWO_PI, K_th, True)
        self.n_coef = self.bx.n * self.by.n * self.bt.n

    def fit(self, T, gx, gy, gt, lam=1e-4):
        """T[ix, iy, it] on the tensor grid -> separable ridge fit."""
        maps = []
        for b, g in ((self.bx, gx), (self.by, gy), (self.bt, gt)):
            Phi = b.design(g)
            maps.append(np.linalg.solve(Phi.T @ Phi + lam * np.eye(b.n), Phi.T))
        self.set_W(np.einsum('ai,bj,ck,ijk->abc', maps[0], maps[1], maps[2], T, optimize=True))

    def set_W(self, W):
        """Set the control weights and precompute every cell's Bernstein coefficients."""
        self.W = W
        ix = self.bx.ctrl_idx(np.arange(self.bx.K))
        iy = self.by.ctrl_idx(np.arange(self.by.K))
        it = self.bt.ctrl_idx(np.arange(self.bt.K))
        Wg = self.W[ix[:, None, None, :, None, None], iy[None, :, None, None, :, None],
                    it[None, None, :, None, None, :]]
        self.C = np.einsum('ip,jq,kr,abcpqr->abcijk', self.bx.Q, self.by.Q, self.bt.Q, Wg,
                           optimize=True)                               # cell Bernstein coeffs
        self.C_t = torch.as_tensor(self.C, device=DEV)

    def eval(self, q, grad=False):
        q = np.atleast_2d(q)
        a, u = self.bx.cell(q[:, 0])
        b, v = self.by.cell(q[:, 1])
        c, w = self.bt.cell(q[:, 2])
        Cc = self.C[a, b, c]
        Bu, Bv, Bw = bern3(u), bern3(v), bern3(w)
        phi = np.einsum('ni,nj,nk,nijk->n', Bu, Bv, Bw, Cc)
        if not grad:
            return phi
        g = np.stack([np.einsum('ni,nj,nk,nijk->n', bern3_d(u) / self.bx.h, Bv, Bw, Cc),
                      np.einsum('ni,nj,nk,nijk->n', Bu, bern3_d(v) / self.by.h, Bw, Cc),
                      np.einsum('ni,nj,nk,nijk->n', Bu, Bv, bern3_d(w) / self.bt.h, Cc)], axis=1)
        return phi, g


# ---------------------------------------------------------------------------------------------
# Certified enclosure level: branch-and-bound for max of phi over M
# ---------------------------------------------------------------------------------------------
def _t(x):
    return torch.as_tensor(x, device=DEV, dtype=torch.float64)


MB_T, MB_INV_T = _t(MB), _t(MB_INV)
BINOM = [[math.comb(k, j) for k in range(4)] for j in range(4)]


def sub_matrix(s0, s1):
    """(R,4,4): cubic Bernstein coefficients on [0,1] -> on the sub-interval [s0, s1]."""
    d = s1 - s0
    T = torch.zeros(len(s0), 4, 4, device=DEV, dtype=torch.float64)
    for j in range(4):
        for k in range(j, 4):
            T[:, j, k] = BINOM[j][k] * s0 ** (k - j) * d ** j
    return MB_INV_T @ T @ MB_T


def bern3_t(s):
    return torch.stack([(1 - s) ** 3, 3 * s * (1 - s) ** 2, 3 * s ** 2 * (1 - s), s ** 3], dim=-1)


def _contains(a0, a1, target):
    return torch.ceil((a0 - target) / TWO_PI) <= torch.floor((a1 - target) / TWO_PI)


class Certifier:
    def __init__(self, field, A, B):
        self.f = field
        self.PA, self.PB = _t(A.bez), _t(B.bez)
        self.rhoA = A.rho

    def phi_t(self, t, th, return_cell=False):
        f = self.f
        vals = []
        for b, x in ((f.bx, t[:, 0]), (f.by, t[:, 1]), (f.bt, th)):
            s = (x - b.lo) / b.h
            c = torch.floor(s)
            if b.periodic:
                vals.append(((c % b.K).long(), s - c))
            else:
                c = c.clamp(0, b.K - 1)
                vals.append((c.long(), (s - c).clamp(0, 1)))
        (a, u), (bb, v), (c, w) = vals
        phi = torch.einsum('ni,nj,nk,nijk->n', bern3_t(u), bern3_t(v), bern3_t(w), f.C_t[a, bb, c])
        return (phi, (a, bb, c)) if return_cell else phi

    def lower_bound(self, R, return_cell=False):
        ia, s0, s1, ib, u0, u1, th0, th1 = R
        a = torch.einsum('ni,nid->nd', bern3_t((s0 + s1) / 2), self.PA[ia])
        b = torch.einsum('ni,nid->nd', bern3_t((u0 + u1) / 2), self.PB[ib])
        th = (th0 + th1) / 2
        c, s = torch.cos(th), torch.sin(th)
        t = b - torch.stack([c * a[:, 0] - s * a[:, 1], s * a[:, 0] + c * a[:, 1]], dim=1)
        return self.phi_t(t, th, return_cell)

    def boxes(self, R):
        ia, s0, s1, ib, u0, u1, th0, th1 = R
        SA = sub_matrix(s0, s1) @ self.PA[ia]                           # (R,4,2) Bezier ctrl
        SB = sub_matrix(u0, u1) @ self.PB[ib]
        r = torch.linalg.norm(SA, dim=-1)
        al = torch.atan2(SA[..., 1], SA[..., 0])
        a0, a1 = al + th0[:, None], al + th1[:, None]
        cx0, cx1, cy0, cy1 = r * torch.cos(a0), r * torch.cos(a1), r * torch.sin(a0), r * torch.sin(a1)
        xmax = torch.where(_contains(a0, a1, 0.0), r, torch.maximum(cx0, cx1))
        xmin = torch.where(_contains(a0, a1, np.pi), -r, torch.minimum(cx0, cx1))
        ymax = torch.where(_contains(a0, a1, np.pi / 2), r, torch.maximum(cy0, cy1))
        ymin = torch.where(_contains(a0, a1, 1.5 * np.pi), -r, torch.minimum(cy0, cy1))
        ra_min = torch.stack([xmin.min(1).values, ymin.min(1).values], dim=1)
        ra_max = torch.stack([xmax.max(1).values, ymax.max(1).values], dim=1)
        tmin = SB.min(1).values - ra_max
        tmax = SB.max(1).values - ra_min
        ext = torch.stack([torch.linalg.norm(SA.max(1).values - SA.min(1).values, dim=1),
                           torch.linalg.norm(SB.max(1).values - SB.min(1).values, dim=1),
                           (th1 - th0) * self.rhoA], dim=1)
        return tmin, tmax, ext

    def upper_bound(self, tmin, tmax, th0, th1, max_cells=3, cell_thr=None):
        """
        Upper bound of phi on the box. With cell_thr (Kx,Ky,Kt), additionally returns whether the
        restricted bound in every covered cell c is below cell_thr[c] (per-cell pruning).
        """
        f = self.f
        axes = []
        for b, lo_v, hi_v in ((f.bx, tmin[:, 0], tmax[:, 0]), (f.by, tmin[:, 1], tmax[:, 1]),
                              (f.bt, th0, th1)):
            f0, f1 = (lo_v - b.lo) / b.h, (hi_v - b.lo) / b.h
            c0 = torch.floor(f0)
            c1 = torch.maximum(torch.ceil(f1) - 1, c0)
            axes.append((b, f0, f1, c0, c1))
        too_big = torch.zeros_like(th0, dtype=torch.bool)
        for (b, f0, f1, c0, c1) in axes:
            too_big |= (c1 - c0 + 1) > max_cells
            if not b.periodic:                          # box leaves the fitted domain: no bound
                too_big |= (f0 < 0) | (f1 > b.K)
        ub = torch.full_like(th0, -np.inf)
        below = torch.ones_like(too_big)
        per_axis = []
        for (b, f0, f1, c0, c1) in axes:
            opts = []
            for o in range(max_cells):
                c = c0 + o
                valid = c <= c1
                S = sub_matrix((f0 - c).clamp(0, 1), (f1 - c).clamp(0, 1))
                ci = (c % b.K) if b.periodic else c.clamp(0, b.K - 1)
                opts.append((ci.long(), valid, S))
            per_axis.append(opts)
        for (cx, vx, Sx) in per_axis[0]:
            for (cy, vy, Sy) in per_axis[1]:
                for (ct, vt, St) in per_axis[2]:
                    valid = vx & vy & vt
                    if not valid.any():
                        continue
                    Cc = f.C_t[cx, cy, ct]
                    Cr = torch.einsum('rip,rjq,rkw,rpqw->rijk', Sx, Sy, St, Cc)
                    m = Cr.amax(dim=(1, 2, 3))
                    ub = torch.where(valid, torch.maximum(ub, m), ub)
                    if cell_thr is not None:
                        below &= ~valid | (m < cell_thr[cx, cy, ct])
        ub = torch.where(too_big, torch.full_like(ub, np.inf), ub)
        if cell_thr is None:
            return ub
        return ub, below & ~too_big

    def certify_cells(self, n_theta0=32, tol=2e-3, chunk=100_000, max_iter=150):
        """
        Per-cell enclosure levels: for every cell c met by M, a certified bound
        l_c >= max_{M in cell c} phi (within tol), by branch-and-bound with per-cell pruning:
        a region is removed once its restricted bound in every covered cell is below that cell's
        best lower bound + tol.  Returns (l_c with -inf where M does not enter, certified, time).
        """
        f = self.f
        K = (f.bx.K, f.by.K, f.bt.K)
        best = torch.full(K, -np.inf, device=DEV, dtype=torch.float64)

        def update(Rc):
            v, (a, b, c) = self.lower_bound(Rc, return_cell=True)
            best.view(-1).scatter_reduce_(0, (a * K[1] + b) * K[2] + c, v, reduce='amax')

        nA, nB = len(self.PA), len(self.PB)
        ia, ib, it = torch.meshgrid(torch.arange(nA, device=DEV), torch.arange(nB, device=DEV),
                                    torch.arange(n_theta0, device=DEV), indexing='ij')
        ia, ib, it = ia.ravel(), ib.ravel(), it.ravel()
        z = torch.zeros(len(ia), device=DEV, dtype=torch.float64)
        o = torch.ones_like(z)
        th0 = it.double() * TWO_PI / n_theta0
        R = [ia, z, o, ib, z.clone(), o.clone(), th0, th0 + TWO_PI / n_theta0]
        t_start = time.perf_counter()
        for s in range(0, len(R[0]), chunk):
            update([x[s:s + chunk] for x in R])
        for _ in range(max_iter):
            if len(R[0]) == 0:
                break
            keep_all, ext_all = [], []
            for s in range(0, len(R[0]), chunk):
                Rc = [x[s:s + chunk] for x in R]
                tmin, tmax, ext = self.boxes(Rc)
                _, below = self.upper_bound(tmin, tmax, Rc[6], Rc[7], cell_thr=best + tol)
                keep_all.append(~below)
                ext_all.append(ext)
            keep = torch.cat(keep_all)
            ext = torch.cat(ext_all)[keep]
            R = [x[keep] for x in R]
            if len(R[0]) == 0:
                break
            dim = ext.argmax(dim=1)
            ia, s0, s1, ib, u0, u1, th0, th1 = R
            sm, um, tm = (s0 + s1) / 2, (u0 + u1) / 2, (th0 + th1) / 2
            c1 = [ia, s0, torch.where(dim == 0, sm, s1), ib, u0, torch.where(dim == 1, um, u1),
                  th0, torch.where(dim == 2, tm, th1)]
            c2 = [ia, torch.where(dim == 0, sm, s0), s1, ib, torch.where(dim == 1, um, u0), u1,
                  torch.where(dim == 2, tm, th0), th1]
            R = [torch.cat([a, b]) for a, b in zip(c1, c2)]
            for s in range(0, len(R[0]), chunk):
                update([x[s:s + chunk] for x in R])
        return (best + tol).cpu().numpy(), len(R[0]) == 0, time.perf_counter() - t_start

    def certify(self, n_theta0=32, tol=2e-3, chunk=200_000, max_iter=80, verbose=True):
        nA, nB = len(self.PA), len(self.PB)
        ia, ib, it = torch.meshgrid(torch.arange(nA, device=DEV), torch.arange(nB, device=DEV),
                                    torch.arange(n_theta0, device=DEV), indexing='ij')
        ia, ib, it = ia.ravel(), ib.ravel(), it.ravel()
        z = torch.zeros(len(ia), device=DEV, dtype=torch.float64)
        o = torch.ones_like(z)
        th0 = it.double() * TWO_PI / n_theta0
        R = [ia, z, o, ib, z.clone(), o.clone(), th0, th0 + TWO_PI / n_theta0]
        best = self.lower_bound(R).max().item()
        n_processed = 0
        t_start = time.perf_counter()
        for k in range(max_iter):
            n = len(R[0])
            if n == 0:
                break
            n_processed += n
            keep_all, ext_all = [], []
            for s in range(0, n, chunk):
                Rc = [x[s:s + chunk] for x in R]
                tmin, tmax, ext = self.boxes(Rc)
                ub = self.upper_bound(tmin, tmax, Rc[6], Rc[7])
                keep_all.append(ub >= best + tol)          # removed only if UB < best + tol (strict)
                ext_all.append(ext)
            keep = torch.cat(keep_all)
            ext = torch.cat(ext_all)[keep]
            R = [x[keep] for x in R]
            if len(R[0]) == 0:
                break
            dim = ext.argmax(dim=1)
            ia, s0, s1, ib, u0, u1, th0, th1 = R
            sm, um, tm = (s0 + s1) / 2, (u0 + u1) / 2, (th0 + th1) / 2
            # child 1 takes the lower half of the split dimension, child 2 the upper half
            c1 = [ia, s0, torch.where(dim == 0, sm, s1), ib, u0, torch.where(dim == 1, um, u1),
                  th0, torch.where(dim == 2, tm, th1)]
            c2 = [ia, torch.where(dim == 0, sm, s0), s1, ib, torch.where(dim == 1, um, u0), u1,
                  torch.where(dim == 2, tm, th0), th1]
            R = [torch.cat([a, b]) for a, b in zip(c1, c2)]
            for s in range(0, len(R[0]), chunk):
                best = max(best, self.lower_bound([x[s:s + chunk] for x in R]).max().item())
            if verbose and k % 5 == 0:
                print(f"    [{self.f.kind}] iter {k:2d}: open regions {len(R[0]):8d} | best {best:+.4f}")
        certified = len(R[0]) == 0
        return best + tol, certified, n_processed, time.perf_counter() - t_start

    # ---------- certificates for the CBF itself ----------
    def restricted(self, cells, lo, hi):
        """Bernstein coefficients (R,4,4,4) of phi on the sub-boxes [lo, hi] (cell-local, in [0,1]^3)."""
        C = self.f.C_t[cells[:, 0], cells[:, 1], cells[:, 2]]
        Sx, Sy, St = (sub_matrix(lo[:, k], hi[:, k]) for k in range(3))
        return torch.einsum('rip,rjq,rkw,rpqw->rijk', Sx, Sy, St, C)

    def _cells_meeting(self, delta):
        """
        (Kx, Ky) masks of the cells that meet the active region {|t|_inf < D - delta} and the
        activation band {D - delta <= |t|_inf <= D}.  Whole cells are used, so both are supersets.
        """
        masks = []
        for b in (self.f.bx, self.f.by):
            e0 = b.lo + b.h * torch.arange(b.K, device=DEV, dtype=torch.float64)
            e1 = e0 + b.h
            r = b.hi - delta
            masks.append(((e1 > -r) & (e0 < r), (e0 <= -r) | (e1 >= r)))
        (ax, bx), (ay, by) = masks
        return ax[:, None] & ay[None, :], bx[:, None] | by[None, :]

    def certify_band(self, level, delta=ACT_DELTA):
        """
        Activation band (Assumption 1): phi > level wherever the CBF is switched off, i.e. on
        D - delta <= |t|_inf <= D, via the minimum Bernstein coefficient of every cell meeting it.
        """
        cmin = self.f.C_t.amin(dim=(3, 4, 5))
        m = cmin[self._cells_meeting(delta)[1]].min().item()
        return m > level, m

    def certify_regular(self, level, delta=ACT_DELTA, max_depth=6, chunk=100_000):
        """
        Regularity: grad phi != 0 on {phi = level} inside the active domain.  A sub-box is cleared
        if level lies outside the range of its Bernstein coefficients (no level set there) or if
        one partial derivative has Bernstein coefficients of one strict sign (never zero there).
        Otherwise it is split into 8. Returns (certified, depth reached, open boxes left).
        """
        C = self.f.C_t
        cand = (C.amin(dim=(3, 4, 5)) <= level) & (C.amax(dim=(3, 4, 5)) >= level)
        cand &= self._cells_meeting(delta)[0][:, :, None]
        cells = torch.nonzero(cand)
        lo = torch.zeros(len(cells), 3, device=DEV, dtype=torch.float64)
        hi = torch.ones_like(lo)
        corners = torch.tensor([[(c >> k) & 1 for k in range(3)] for c in range(8)],
                               device=DEV, dtype=torch.bool)
        for depth in range(max_depth + 1):
            keep = []
            for s in range(0, len(cells), chunk):
                Cr = self.restricted(cells[s:s + chunk], lo[s:s + chunk], hi[s:s + chunk])
                has_level = (Cr.amin(dim=(1, 2, 3)) <= level) & (Cr.amax(dim=(1, 2, 3)) >= level)
                regular = torch.zeros_like(has_level)
                for d in range(3):
                    D = torch.diff(Cr, dim=1 + d)                   # derivative coefficients (up to a
                    regular |= (D > 0).flatten(1).all(1) | (D < 0).flatten(1).all(1)  # positive factor)
                keep.append(has_level & ~regular)
            keep = torch.cat(keep) if keep else torch.zeros(0, dtype=torch.bool, device=DEV)
            cells, lo, hi = cells[keep], lo[keep], hi[keep]
            if len(cells) == 0 or depth == max_depth:
                break
            mid = (lo + hi) / 2
            cells = cells.repeat(8, 1)
            c = corners.repeat_interleave(len(lo), dim=0)
            lo, hi, mid = lo.repeat(8, 1), hi.repeat(8, 1), mid.repeat(8, 1)
            lo, hi = torch.where(c, mid, lo), torch.where(c, hi, mid)
        return len(cells) == 0, depth, len(cells)


# ---------------------------------------------------------------------------------------------
# Geometry for the ground-truth distance
# ---------------------------------------------------------------------------------------------
def min_enclosing_circle(P, seed=0):
    """Smallest circle containing the points P (randomized incremental algorithm, expected linear
    time). Returns (center, radius)."""
    P = np.asarray(P, dtype=np.float64)[np.random.default_rng(seed).permutation(len(P))]

    def circ2(a, b):
        c = (a + b) / 2
        return c, np.linalg.norm(a - c)

    def circ3(a, b, c):
        d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
        if abs(d) < 1e-14:                       # collinear: widest pair
            return max((circ2(a, b), circ2(a, c), circ2(b, c)), key=lambda t: t[1])
        sa, sb, sc = a @ a, b @ b, c @ c
        o = np.array([sa * (b[1] - c[1]) + sb * (c[1] - a[1]) + sc * (a[1] - b[1]),
                      sa * (c[0] - b[0]) + sb * (a[0] - c[0]) + sc * (b[0] - a[0])]) / d
        return o, np.linalg.norm(a - o)

    def inside(c, p):
        return np.linalg.norm(p - c[0]) <= c[1] * (1 + 1e-12) + 1e-12

    c = (P[0], 0.0)
    for i in range(1, len(P)):
        if inside(c, P[i]):
            continue
        c = (P[i], 0.0)
        for j in range(i):
            if inside(c, P[j]):
                continue
            c = circ2(P[i], P[j])
            for k in range(j):
                if not inside(c, P[k]):
                    c = circ3(P[i], P[j], P[k])
    return c[0], float(c[1])


def max_edge(P):
    return float(np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1).max())


def polygons_cross(PA, PB):
    """True if an edge of polygon PA meets an edge of polygon PB (orientation tests; collinear
    overlaps are caught by the bounding-box test)."""
    PA, PB = torch.as_tensor(PA, device=DEV), torch.as_tensor(PB, device=DEV)
    a, b = PA[:, None], torch.roll(PA, -1, 0)[:, None]
    c, d = PB[None], torch.roll(PB, -1, 0)[None]

    def cross(o, p, q):
        return (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])
    o1, o2, o3, o4 = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
    boxes = ((torch.minimum(a, b) <= torch.maximum(c, d)) & (torch.minimum(c, d) <= torch.maximum(a, b))).all(-1)
    return bool(((o1 * o2 <= 0) & (o3 * o4 <= 0) & boxes).any())


def true_distance(PA_world, PB_world):
    """
    Distance between two polygons if they are disjoint; otherwise a value <= 0 (minus the deepest
    vertex penetration, or 0 if the boundaries cross with no vertex inside the other polygon).
    For disjoint polygons the vertex-to-polygon distances give the exact set distance (the closest
    points of two disjoint segments include an endpoint).  A crossing puts a vertex within one edge
    length of the other boundary, so the edge test is run whenever the vertex distance is not larger
    than the longest edge.
    """
    def edges(P):
        P = torch.as_tensor(P, device=DEV)
        return P[None], torch.roll(P, -1, 0)[None]
    Aa, Ab = edges(PA_world)
    Ba, Bb = edges(PB_world)
    dA = polygon_sdf(torch.as_tensor(PA_world, device=DEV)[:, None, :], Ba, Bb)[0].min()
    dB = polygon_sdf(torch.as_tensor(PB_world, device=DEV)[:, None, :], Aa, Ab)[0].min()
    d = min(dA.item(), dB.item())
    if 0 < d <= max(max_edge(PA_world), max_edge(PB_world)) and polygons_cross(PA_world, PB_world):
        return 0.0
    return d


def step_gap(GT, rho_gt, x0, x1, pairs=None, tol=1e-4, max_depth=14):
    """
    Minimum EVALUATED polygon distance over one integration step, poses interpolated linearly between
    x0 and x1 (both endpoints included).  On a sub-interval of relative length s, no point of body i
    moves by more than s * m_i with m_i = |dp_i| + |dtheta_i| rho_i, so a pair with distance
    d > s (m_i + m_j) at the start of the sub-interval cannot touch within it; otherwise the
    sub-interval is bisected (down to a residual motion below tol).
    This return value is not a lower bound on the continuous-time minimum. A
    tolerance/depth stop alone is not a separation certificate. If every evaluated
    distance exceeds max(tol, max_pair_motion / 2**max_depth), all leaves with
    positive endpoints satisfy the motion-based separation test as well.
    """
    dx = x1 - x0
    m = np.linalg.norm(dx[:, :2], axis=1) + np.abs(dx[:, 2]) * rho_gt
    out = np.inf

    def gap(i, j, s):
        x = x0 + s * dx
        dc = np.linalg.norm(x[i, :2] - x[j, :2])
        if dc >= rho_gt[i] + rho_gt[j] + 0.1:
            return dc - rho_gt[i] - rho_gt[j]
        return true_distance(GT[i] @ rot(x[i, 2]).T + x[i, :2], GT[j] @ rot(x[j, 2]).T + x[j, :2])

    for i, j in (pairs if pairs is not None else itertools.combinations(range(len(x0)), 2)):
        mij = m[i] + m[j]
        d0 = gap(i, j, 0.0) if np.linalg.norm(x0[i, :2] - x0[j, :2]) < rho_gt[i] + rho_gt[j] + 0.2 + mij else np.inf
        if not np.isfinite(d0):
            continue
        d1 = gap(i, j, 1.0)
        out = min(out, d0, d1)
        stack = [(0.0, 1.0, d0, d1, 0)]
        while stack:
            a, b, da, db, depth = stack.pop()
            w = (b - a) * mij
            if min(da, db) <= 0 or max(da, db) > w or w < tol or depth >= max_depth:
                continue                    # witness / separation proof / tolerance or work limit
            c = (a + b) / 2
            dc_ = gap(i, j, c)
            out = min(out, dc_)
            stack += [(a, c, da, dc_, depth + 1), (c, b, dc_, db, depth + 1)]
    return out


# ---------------------------------------------------------------------------------------------
# Closed loop: SE(2) single integrators, joint CBF-QP
# ---------------------------------------------------------------------------------------------
J2 = np.array([[0.0, -1.0], [1.0, 0.0]])


def wrap(a):
    return (a + np.pi) % TWO_PI - np.pi


def simulate(field, l_star, A, B, D, steps, dt, start_A, goal_A, start_B, goal_B, a_static=False,
             v_max=1.0, w_max=2.0, gamma=5.0, k_v=1.5, k_w=2.0):
    """Joint CBF-QP over (v_A, omega_A, v_B, omega_B); with a_static, A never moves or rotates."""
    xA, xB = start_A.copy(), start_B.copy()
    free = slice(3, 6) if a_static else slice(0, 6)
    n_sides = 12
    ang = np.linspace(0, TWO_PI, n_sides, endpoint=False)
    nk = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    G_box, h_box = [], []
    for off in (0, 3):                                                   # |v| polygon, |omega|
        for n_ in nk:
            r = np.zeros(6)
            r[off:off + 2] = -n_
            G_box.append(r)
            h_box.append(-v_max * np.cos(np.pi / n_sides))
        for sgn in (1, -1):
            r = np.zeros(6)
            r[off + 2] = -sgn
            G_box.append(r)
            h_box.append(-w_max)
    G_box, h_box = np.array(G_box), np.array(h_box)

    log = dict(xA=[], xB=[], h=[], dist=[], u=[], t_ctrl=[], active=[], step_min=[])
    reached = None
    margin = ACT_DELTA
    for k in range(steps):
        t0 = time.perf_counter()
        uA = np.r_[k_v * (goal_A[:2] - xA[:2]), k_w * wrap(goal_A[2] - xA[2])]
        uB = np.r_[k_v * (goal_B[:2] - xB[:2]), k_w * wrap(goal_B[2] - xB[2])]
        for u_ in (uA, uB):
            s = np.linalg.norm(u_[:2])
            if s > v_max:
                u_[:2] *= v_max / s
            u_[2] = np.clip(u_[2], -w_max, w_max)
        u_nom = np.r_[uA, uB]

        RBt = rot(-xB[2])
        t = RBt @ (xA[:2] - xB[:2])
        th = (xA[2] - xB[2]) % TWO_PI
        near = np.all(np.abs(t) < D - margin)
        if near:
            phi, g = field.eval(np.r_[t, th], grad=True)
            phi, g = phi[0], g[0]
            h = phi - l_star
            gt, gth = g[:2], g[2]
            row = np.zeros(6)
            row[0:2] = rot(xB[2]) @ gt                                   # d/dv_A
            row[3:5] = -rot(xB[2]) @ gt                                  # d/dv_B
            row[2] = gth                                                 # d/domega_A
            row[5] = -gt @ (J2 @ t) - gth                                # d/domega_B
            G = np.vstack([row, G_box])
            hh = np.r_[-gamma * h, h_box]
        else:
            h = np.nan
            G, hh = G_box, h_box
        u = np.zeros(6)
        u[free], _ = solve_ldp_qp(u_nom[free], G[:, free], hh, n_soft=1 if near else 0)
        log['t_ctrl'].append(time.perf_counter() - t0)

        log['xA'].append(xA.copy())
        log['xB'].append(xB.copy())
        log['h'].append(h)
        log['u'].append(u.copy())
        log['active'].append(near)
        log['dist'].append(true_distance(A.world(xA), B.world(xB)))
        x_old = np.stack([xA, xB])
        xA = xA + dt * u[:3]
        xB = xB + dt * u[3:]
        log['step_min'].append(step_gap([A.dense, B.dense], np.array([A.rho, B.rho]), x_old, np.stack([xA, xB])))
        if reached is None and max(np.linalg.norm(xA[:2] - goal_A[:2]), np.linalg.norm(xB[:2] - goal_B[:2])) < 0.03 \
                and max(abs(wrap(xA[2] - goal_A[2])), abs(wrap(xB[2] - goal_B[2]))) < 0.03:
            reached = k * dt
    for key in log:
        log[key] = np.array(log[key])
    log['reached'] = reached
    log['goals'] = (goal_A, goal_B)
    return log


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=900)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--gif-every', type=int, default=5)
    ap.add_argument('--scenario', choices=['switch', 'dock'], default='switch',
                    help='switch: swap positions (both rotate); dock: C-shape A fixed, B pushes '
                         'one lobe into the mouth of A')
    ap.add_argument('--dock-gap', type=float, default=-0.10,
                    help='true gap at the nominal dock goal [m]; negative = the goal itself is in '
                         'collision, so the CBF must stop B at its certified boundary')
    args = ap.parse_args()

    A, B = shape_A(), shape_B()
    D = round(A.rho + B.rho + 0.3, 2)
    print(f"rho_A {A.rho:.3f}  rho_B {B.rho:.3f}  |  C-space domain t in [-{D}, {D}]^2 x S^1")

    # ---------- training / test data (raster C-space slices) ----------
    t0 = time.perf_counter()
    n_th_train, n_xy = 96, 100
    th_train = np.arange(n_th_train) * TWO_PI / n_th_train
    th_test = (np.arange(16) + 0.37) * TWO_PI / 16
    g, S_train = cspace_sdf_slices(A, B, th_train, D)
    _, S_test = cspace_sdf_slices(A, B, th_test, D)
    gxy = np.linspace(-D, D, n_xy)
    T = sample_slices(g, S_train, gxy, gxy).transpose(1, 2, 0)          # [ix, iy, it]
    print(f"raster C-space slices: {time.perf_counter() - t0:.1f} s")

    rng = np.random.default_rng(0)
    Xg, Yg = np.meshgrid(g, g, indexing='ij')
    test_q, test_sdf = [], []
    for k, th in enumerate(th_test):
        sel = rng.choice(Xg.size, 4000, replace=False)
        test_q.append(np.stack([Xg.ravel()[sel], Yg.ravel()[sel], np.full(4000, th)], axis=1))
        test_sdf.append(S_test[k].ravel()[sel])
    test_q, test_sdf = np.concatenate(test_q), np.concatenate(test_sdf)
    near_bd = np.abs(test_sdf) < 0.1

    configs = {'bspline': dict(K_xy=45, K_th=48), 'bezier': dict(K_xy=16, K_th=16)}
    res = {}
    for kind, cfg in configs.items():
        f = Field3D(kind, D, cfg['K_xy'], cfg['K_th'])
        t0 = time.perf_counter()
        f.fit(T, gxy, gxy, th_train)
        t_fit = time.perf_counter() - t0
        err = f.eval(test_q) - test_sdf

        # gradient jump across x-faces of the cells (C^0 vs C^2)
        qs = test_q[near_bd][:3000].copy()
        cx, _ = f.bx.cell(qs[:, 0])
        qs[:, 0] = f.bx.lo + np.clip(cx, 1, f.bx.K - 1) * f.bx.h
        _, gm = f.eval(qs - [1e-7, 0, 0], grad=True)
        _, gp = f.eval(qs + [1e-7, 0, 0], grad=True)
        jump = np.linalg.norm(gp - gm, axis=1).max()

        print(f"\n[{kind}] coefficients {f.n_coef} (cells {f.bx.K}x{f.by.K}x{f.bt.K}) | fit {t_fit:.2f} s")
        print(f"[{kind}] certifying enclosure (branch-and-bound over sigma, tau, theta) ...")
        cert = Certifier(f, A, B)
        l_star, ok, n_reg, t_cert = cert.certify()
        # independent sanity check on random points of M
        sg = rng.uniform(0, 1, (200000, 3))
        ia_, ib_ = rng.integers(0, len(A.bez), 200000), rng.integers(0, len(B.bez), 200000)
        a = np.einsum('ni,nid->nd', bern3(sg[:, 0]), A.bez[ia_])
        b = np.einsum('ni,nid->nd', bern3(sg[:, 1]), B.bez[ib_])
        th = sg[:, 2] * TWO_PI
        tM = b - np.stack([np.cos(th) * a[:, 0] - np.sin(th) * a[:, 1],
                           np.sin(th) * a[:, 0] + np.cos(th) * a[:, 1]], axis=1)
        sampled_max = f.eval(np.column_stack([tM, th])).max()
        print(f"[{kind}] certified l* = {l_star:+.4f} ({'certified' if ok else 'NOT certified'}) | "
              f"{n_reg} regions, {t_cert:.1f} s | sampled max of phi on M = {sampled_max:+.4f}")
        t0 = time.perf_counter()
        band_ok, band_min = cert.certify_band(l_star)
        reg_ok, reg_depth, reg_open = cert.certify_regular(l_star)
        print(f"[{kind}] activation band: min phi = {band_min:+.4f} > l* : {band_ok} | "
              f"regularity grad phi != 0 on {{phi = l*}}: {reg_ok} (depth {reg_depth}, open {reg_open}) "
              f"| {time.perf_counter() - t0:.1f} s")
        res[kind] = dict(field=f, l=l_star, ok=ok, n_reg=n_reg, t_cert=t_cert, t_fit=t_fit,
                         band=band_ok, reg=reg_ok,
                         rmse=np.sqrt(np.mean(err ** 2)), rmse_bd=np.sqrt(np.mean(err[near_bd] ** 2)),
                         maxerr_bd=np.abs(err[near_bd]).max(), jump=jump, sampled_max=sampled_max)

    # ---------- closed loop ----------
    if args.scenario == 'switch':
        poses = (np.array([-2.4, 0.18, 0.0]), np.array([2.4, 0.18, np.pi / 2]),
                 np.array([2.4, -0.18, 0.0]), np.array([-2.4, -0.18, -np.pi / 3]))
    else:
        # A fixed at the origin (mouth facing +x); B turned so one lobe points into the mouth,
        # pushed in as deep as possible while keeping the true gap >= dock_gap
        pA = np.zeros(3)
        th_B = np.pi / 3
        d_goal = next(d for d in np.arange(1.5, 0.0, -0.005)
                      if true_distance(A.world(pA), B.world(np.array([d, 0.0, th_B]))) < args.dock_gap) + 0.005
        tip = d_goal - np.abs(B.world(np.array([0.0, 0.0, th_B]))[:, 0].min())
        print(f"\ndock goal: B at x = {d_goal:.3f}, lobe tip at x = {tip:+.3f} "
              f"(A's inner wall at x = {A.world(pA)[np.abs(A.world(pA)[:, 1]) < 0.03][:, 0].max():+.3f})")
        poses = (pA, pA, np.array([2.6, 0.25, -0.6]), np.array([d_goal, 0.0, th_B]))
    for kind in res:
        t0 = time.perf_counter()
        res[kind]['log'] = simulate(res[kind]['field'], res[kind]['l'], A, B, D, args.steps, args.dt,
                                    *poses, a_static=(args.scenario == 'dock'))
        L = res[kind]['log']
        u = L['u']
        res[kind]['min_dist'] = min(L['dist'].min(), L['step_min'].min())
        res[kind]['min_h'] = np.nanmin(L['h'])
        res[kind]['tv'] = np.abs(np.diff(u, axis=0)).sum()
        res[kind]['tv_w'] = np.abs(np.diff(u[:, [2, 5]], axis=0)).sum()
        res[kind]['path'] = (np.linalg.norm(np.diff(L['xA'][:, :2], axis=0), axis=1).sum()
                             + np.linalg.norm(np.diff(L['xB'][:, :2], axis=0), axis=1).sum())
        res[kind]['t_ctrl'] = np.median(L['t_ctrl'][L['active']]) * 1e3 if L['active'].any() else np.nan
        print(f"[{kind}] simulated {args.steps} steps in {time.perf_counter() - t0:.1f} s")

    # ---------- table ----------
    rows = [('coefficients', 'n_coef', '{:d}'),
            ('fit time [s]', 't_fit', '{:.2f}'),
            ('test RMSE, all [m]', 'rmse', '{:.4f}'),
            ('test RMSE, |sdf|<0.1 [m]', 'rmse_bd', '{:.4f}'),
            ('test max err, |sdf|<0.1 [m]', 'maxerr_bd', '{:.4f}'),
            ('certified l* [m]  (lower = less conservative)', 'l', '{:+.4f}'),
            ('certification time [s]', 't_cert', '{:.1f}'),
            ('activation band certified (phi > l*)', 'band', '{}'),
            ('regularity certified (grad phi != 0 on level set)', 'reg', '{}'),
            ('B&B regions processed', 'n_reg', '{:d}'),
            ('max gradient jump across cell faces', 'jump', '{:.3e}'),
            ('min true distance A-B [m]', 'min_dist', '{:+.4f}'),
            ('min h = phi - l*', 'min_h', '{:+.4f}'),
            ('reached both goals at [s]', 'reached', '{:.2f}'),
            ('final true gap A-B [m]', 'final_gap', '{:+.4f}'),
            ('final x of B (smaller = deeper)', 'final_xB', '{:.4f}'),
            ('total path length [m]', 'path', '{:.2f}'),
            ('input total variation (all)', 'tv', '{:.2f}'),
            ('input total variation (omega)', 'tv_w', '{:.2f}'),
            ('online CBF+QP per step, median [ms]', 't_ctrl', '{:.3f}')]
    for kind in res:
        res[kind]['n_coef'] = res[kind]['field'].n_coef
        res[kind]['reached'] = res[kind]['log']['reached']
        res[kind]['final_gap'] = res[kind]['log']['dist'][-1]
        res[kind]['final_xB'] = res[kind]['log']['xB'][-1, 0]
    print("\n" + "=" * 86)
    print(f"{'metric':50s}{'B-spline (C2)':>18s}{'piecewise Bezier (C0)':>22s}")
    print("-" * 86)
    for name, key, fmt in rows:
        vals = []
        for kind in ('bspline', 'bezier'):
            v = res[kind][key]
            vals.append('n/a' if v is None else fmt.format(v))
        print(f"{name:50s}{vals[0]:>18s}{vals[1]:>22s}")
    print("=" * 86)

    make_figures(res, A, B, D, g, th_test, S_test, args)


def make_figures(res, A, B, D, g, th_test, S_test, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from PIL import Image

    # ---------- C-space slices: true CO vs certified sublevel sets ----------
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.4))
    xs = np.linspace(-D, D, 220)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    for ax, k in zip(axs, (0, 4, 8, 12)):
        th = th_test[k]
        ax.contourf(g, g, (S_test[k] <= 0).T.astype(float), levels=[0.5, 1.5], colors=['0.6'])
        for kind, col, ls in (('bspline', 'tab:blue', '-'), ('bezier', 'tab:orange', '--')):
            f = res[kind]['field']
            Z = f.eval(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, th)])).reshape(X.shape)
            ax.contour(X, Y, Z, levels=[res[kind]['l']], colors=col, linestyles=ls, linewidths=1.6)
        ax.set_title(f"theta = {np.degrees(th):.0f} deg", fontsize=10)
        ax.set_aspect('equal')
        ax.set_xlim(-D, D)
        ax.set_ylim(-D, D)
    axs[0].plot([], [], color='0.6', lw=6, label='true C-obstacle slice')
    axs[0].plot([], [], color='tab:blue', label=f"B-spline  phi = l* ({res['bspline']['l']:+.3f})")
    axs[0].plot([], [], color='tab:orange', ls='--', label=f"Bezier  phi = l* ({res['bezier']['l']:+.3f})")
    axs[0].legend(loc='lower left', fontsize=7)
    fig.suptitle("C-space slices (held-out angles): certified boundaries {phi = l*}", fontsize=11)
    fig.tight_layout()
    fig.savefig('cspace_fit_comparison.png', dpi=120)
    plt.close(fig)
    print(">>> saved cspace_fit_comparison.png")

    # ---------- GIF ----------
    kinds = ('bspline', 'bezier')
    titles = {'bspline': 'B-spline C-space SDF (C2)', 'bezier': 'piecewise Bezier C-space SDF (C0)'}
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1])
    ax_s = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    ax_d = fig.add_subplot(gs[1, 0])
    ax_w = fig.add_subplot(gs[1, 1])
    n = len(res['bspline']['log']['h'])
    w_col, w_name = (2, 'omega_A') if args.scenario == 'switch' else (5, 'omega_B')
    tt = np.arange(n) * args.dt
    frames = []
    for k in list(range(0, n, args.gif_every)) + [n - 1]:
        for ax, kind in zip(ax_s, kinds):
            L = res[kind]['log']
            ax.clear()
            gA, gB = L['goals']
            ax.add_patch(Polygon(A.world(gA), closed=True, fill=False, ec='royalblue', ls=':', lw=1))
            ax.add_patch(Polygon(B.world(gB), closed=True, fill=False, ec='crimson', ls=':', lw=1))
            for j in range(0, k, 60):
                ax.add_patch(Polygon(A.world(L['xA'][j]), closed=True, fill=False, ec='royalblue', alpha=0.18))
                ax.add_patch(Polygon(B.world(L['xB'][j]), closed=True, fill=False, ec='crimson', alpha=0.18))
            ax.plot(L['xA'][:k + 1, 0], L['xA'][:k + 1, 1], color='royalblue', lw=1)
            ax.plot(L['xB'][:k + 1, 0], L['xB'][:k + 1, 1], color='crimson', lw=1)
            ax.add_patch(Polygon(A.world(L['xA'][k]), closed=True, color='midnightblue', alpha=0.9))
            ax.add_patch(Polygon(B.world(L['xB'][k]), closed=True, color='darkred', alpha=0.9))
            h = L['h'][k]
            ax.set_title(f"{titles[kind]}\nt = {k * args.dt:.2f} s | h = "
                         f"{'--' if np.isnan(h) else f'{h:+.3f}'} | true dist = {L['dist'][k]:+.3f} m",
                         fontsize=10)
            if args.scenario == 'switch':
                ax.set_xlim(-3.2, 3.2)
                ax.set_ylim(-2.0, 2.0)
            else:
                ax.set_xlim(-1.0, 3.3)
                ax.set_ylim(-1.4, 1.4)
            ax.set_aspect('equal')
        ax_d.clear()
        ax_w.clear()
        for kind, col in (('bspline', 'tab:blue'), ('bezier', 'tab:orange')):
            L = res[kind]['log']
            ax_d.plot(tt[:k + 1], L['dist'][:k + 1], color=col, label=f"{kind}: true distance")
            ax_d.plot(tt[:k + 1], L['h'][:k + 1], color=col, ls=':', label=f"{kind}: h")
            ax_w.plot(tt[:k + 1], L['u'][:k + 1, w_col], color=col, label=f"{kind}: {w_name}")
        ax_d.axhline(0, color='k', lw=0.8)
        ax_d.set_xlim(0, tt[-1])
        ax_d.set_ylim(-0.1, 1.0)
        ax_d.set_xlabel('t [s]')
        ax_d.legend(fontsize=7, loc='upper right')
        ax_w.set_xlim(0, tt[-1])
        ax_w.set_ylim(-2.2, 2.2)
        ax_w.set_xlabel('t [s]')
        ax_w.set_ylabel(f'{w_name} [rad/s]')
        ax_w.legend(fontsize=7, loc='upper right')
        fig.tight_layout()
        fig.canvas.draw()
        w, hh = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(hh, w, 4)
        frames.append(Image.fromarray(buf, 'RGBA').convert('RGB'))
    out = 'cspace_cbf_switching.gif' if args.scenario == 'switch' else 'cspace_cbf_dock.gif'
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=50, loop=0)
    print(f">>> saved {out} ({len(frames)} frames)")


if __name__ == '__main__':
    main()
