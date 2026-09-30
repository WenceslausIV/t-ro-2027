"""
Experiments for the T-RO draft (tro/main.tex).

    python planar/cspace_experiments.py cache      # 5 random shapes + 10 certified pair fields (seed 3)
    python planar/cspace_experiments.py repr       # representations incl. the global Bernstein polynomial
    python planar/cspace_experiments.py stats      # random start/goal trials: ours vs closest-point vs circle
    python planar/cspace_experiments.py unicycle   # unicycle dynamics, 2 and 4 robots
    python planar/cspace_experiments.py static     # one field for a robot among m static obstacles
    python planar/cspace_experiments.py all

Results go to results/*.json, figures to tro/figs/*.png.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time

import numpy as np
import torch
from matplotlib.path import Path
from scipy import ndimage
from scipy.signal import fftconvolve

import cspace_sdf_cbf_compare as C
from cspace_cbf_5robots import enclose, pair_field, random_shape
from sdf_cbf_utils import polygon_sdf, solve_ldp_qp

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(HERE, 'cache')
RES = os.path.join(HERE, 'results')
FIGS = os.path.join(HERE, 'tro', 'figs')
DEV = C.DEV
TWO_PI = C.TWO_PI
J2 = np.array([[0.0, -1.0], [1.0, 0.0]])


def dump(name, obj):
    os.makedirs(RES, exist_ok=True)

    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.bool_):
            return bool(o)
        return str(o)
    with open(os.path.join(RES, name + '.json'), 'w') as fh:
        json.dump(obj, fh, indent=1, default=conv)


# =============================================================================================
# cache: 5 random shapes (same generator and seed as cspace_cbf_5robots.py) and pair fields
# =============================================================================================
def build_cache(seed=3, N=5):
    os.makedirs(CACHE, exist_ok=True)
    rng = np.random.default_rng(seed)
    GT, ctrl, offs = [], [], []
    for _ in range(N):
        G = random_shape(rng)
        S, off = enclose(G)
        GT.append(G)
        ctrl.append(S.ctrl)
        offs.append(off)
    shapes = [C.SmoothShape(c) for c in ctrl]
    data = {'GT': np.array(GT), 'ctrl': np.array(ctrl)}
    meta = {'offsets': offs, 'pairs': {}}
    for i, j in itertools.combinations(range(N), 2):
        t0 = time.perf_counter()
        f, D, l_star, ok, info = pair_field(shapes[i], shapes[j])
        data[f'W_{i}_{j}'] = f.W
        meta['pairs'][f'{i}_{j}'] = dict(D=D, Kxy=f.bx.K, Kth=f.bt.K, l=l_star, ok=bool(ok),
                                         band_min=info[2], time=time.perf_counter() - t0)
        print(f"cache pair ({i},{j}): l* {l_star:+.4f} ok {ok} ({time.perf_counter() - t0:.1f} s)")
    np.savez(os.path.join(CACHE, 'shapes5.npz'), **data)
    with open(os.path.join(CACHE, 'shapes5.json'), 'w') as fh:
        json.dump(meta, fh, indent=1)


def load_cache():
    d = np.load(os.path.join(CACHE, 'shapes5.npz'))
    meta = json.load(open(os.path.join(CACHE, 'shapes5.json')))
    shapes = [C.SmoothShape(c) for c in d['ctrl']]
    fields = {}
    for key, m in meta['pairs'].items():
        i, j = map(int, key.split('_'))
        f = C.Field3D('bspline', m['D'], m['Kxy'], m['Kth'])
        f.set_W(d[f'W_{key}'])
        fields[(i, j)] = (f, m['D'], m['l'])
    return list(d['GT']), shapes, fields, meta


# =============================================================================================
# representation study: B-spline vs piecewise Bezier vs global Bernstein polynomial
# =============================================================================================
class GlobalBP:
    """Conference representation: one tensor Bernstein polynomial of degree Q-1 per axis."""

    def __init__(self, D, Q):
        self.D, self.Q = D, Q
        self.binom = np.array([math.comb(Q - 1, k) for k in range(Q)], dtype=np.float64)
        self.n_coef = Q ** 3

    def basis(self, x, lo, hi):
        xi = np.clip((np.asarray(x, dtype=np.float64) - lo) / (hi - lo), 0, 1)[:, None]
        k = np.arange(self.Q)
        return self.binom * xi ** k * (1 - xi) ** (self.Q - 1 - k)

    def fit(self, T, gx, gy, gt, lam=1e-4):
        maps = []
        for g, lo, hi in ((gx, -self.D, self.D), (gy, -self.D, self.D), (gt, 0.0, TWO_PI)):
            Phi = self.basis(g, lo, hi)
            maps.append(np.linalg.solve(Phi.T @ Phi + lam * np.eye(self.Q), Phi.T))
        self.W = np.einsum('ai,bj,ck,ijk->abc', maps[0], maps[1], maps[2], T, optimize=True)

    def eval(self, q, chunk=4000):
        out = []
        for s in range(0, len(q), chunk):
            qq = q[s:s + chunk]
            Bx = self.basis(qq[:, 0], -self.D, self.D)
            By = self.basis(qq[:, 1], -self.D, self.D)
            Bt = self.basis(qq[:, 2] % TWO_PI, 0.0, TWO_PI)
            tmp = (Bx @ self.W.reshape(self.Q, -1)).reshape(-1, self.Q, self.Q)
            out.append(np.einsum('nbc,nb,nc->n', tmp, By, Bt))
        return np.concatenate(out)


def sample_M(A, B, n, rng):
    """Random points of M = {(b(tau) - R(theta) a(sigma), theta)}."""
    sg = rng.uniform(0, 1, (n, 3))
    ia, ib = rng.integers(0, len(A.bez), n), rng.integers(0, len(B.bez), n)
    a = np.einsum('ni,nid->nd', C.bern3(sg[:, 0]), A.bez[ia])
    b = np.einsum('ni,nid->nd', C.bern3(sg[:, 1]), B.bez[ib])
    th = sg[:, 2] * TWO_PI
    t = b - np.stack([np.cos(th) * a[:, 0] - np.sin(th) * a[:, 1],
                      np.sin(th) * a[:, 0] + np.cos(th) * a[:, 1]], axis=1)
    return np.column_stack([t, th])


def exp_repr():
    A, B = C.shape_A(), C.shape_B()
    D = round(A.rho + B.rho + 0.3, 2)
    th_train = np.arange(96) * TWO_PI / 96
    th_test = (np.arange(16) + 0.37) * TWO_PI / 16
    g, S_train = C.cspace_sdf_slices(A, B, th_train, D)
    _, S_test = C.cspace_sdf_slices(A, B, th_test, D)
    gxy = np.linspace(-D, D, 100)
    T = C.sample_slices(g, S_train, gxy, gxy).transpose(1, 2, 0)
    rng = np.random.default_rng(0)
    Xg, Yg = np.meshgrid(g, g, indexing='ij')
    tq, ts = [], []
    for k, th in enumerate(th_test):
        sel = rng.choice(Xg.size, 4000, replace=False)
        tq.append(np.stack([Xg.ravel()[sel], Yg.ravel()[sel], np.full(4000, th)], axis=1))
        ts.append(S_test[k].ravel()[sel])
    tq, ts = np.concatenate(tq), np.concatenate(ts)
    near = np.abs(ts) < 0.1
    qM = sample_M(A, B, 200_000, rng)

    out, phiM = {}, {}
    for name, make in (('bspline', lambda: C.Field3D('bspline', D, 45, 48)),
                       ('bezier', lambda: C.Field3D('bezier', D, 16, 16)),
                       ('globalbp48', lambda: GlobalBP(D, 48)),
                       ('globalbp23', lambda: GlobalBP(D, 23))):
        f = make()
        t0 = time.perf_counter()
        f.fit(T, gxy, gxy, th_train)
        t_fit = time.perf_counter() - t0
        err = f.eval(tq) - ts
        pm = f.eval(qM)
        t0 = time.perf_counter()
        for _ in range(200):
            f.eval(qM[:1])
        t_eval = (time.perf_counter() - t0) / 200 * 1e6
        r = dict(n_coef=f.n_coef, t_fit=t_fit, rmse=np.sqrt(np.mean(err ** 2)),
                 rmse_bd=np.sqrt(np.mean(err[near] ** 2)), maxerr_bd=np.abs(err[near]).max(),
                 sampled_max_M=pm.max(), eval_us=t_eval)
        if isinstance(f, C.Field3D):
            cert = C.Certifier(f, A, B)
            l_star, ok, n_reg, t_cert = cert.certify(verbose=False)
            band_ok, band_min = cert.certify_band(l_star)
            reg_ok, _, _ = cert.certify_regular(l_star)
            r.update(l_star=l_star, certified=ok, n_reg=n_reg, t_cert=t_cert, band=band_ok,
                     band_min=band_min, regular=reg_ok)
        phiM[name] = pm
        out[name] = r
        print(name, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    dump('repr', out)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.4))
    bins = np.linspace(-0.9, 0.06, 120)
    style = {'bspline': ('tab:blue', 'B-spline ($C^2$)'), 'bezier': ('tab:orange', 'piecewise Bézier ($C^0$)'),
             'globalbp48': ('tab:green', 'global Bernstein, $Q=48$'), 'globalbp23': ('tab:gray', 'global Bernstein, $Q=23$')}
    for name, (col, lab) in style.items():
        ax.hist(phiM[name], bins=bins, histtype='step', color=col, lw=1.4, label=lab)
        lv = out[name].get('l_star', out[name]['sampled_max_M'])
        ax.axvline(lv, color=col, ls='--' if 'l_star' in out[name] else ':', lw=1.2)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\phi$ on $M$ (200k random points)  [m]')
    ax.set_ylabel('count')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_title(r'dashed: certified $l^\star$;  dotted: sampled $\max_{M}\phi$ (lower bound on any $l$)', fontsize=8)
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, 'M_histogram.png'), dpi=150)
    print('>>> saved tro/figs/M_histogram.png')


# =============================================================================================
# generic multi-robot closed loop (pairwise constraints)
# =============================================================================================
def closest_pair(PA, PB):
    """Closest points between two polygons (torch, world frame). Returns d, pa, pb, n (B -> A)."""
    ea = (PA[None], torch.roll(PA, -1, 0)[None])
    eb = (PB[None], torch.roll(PB, -1, 0)[None])
    sB, gB = polygon_sdf(PA[:, None, :], *eb)
    sA, gA = polygon_sdf(PB[:, None, :], *ea)
    sB, gB, sA, gA = sB[:, 0], gB[:, 0], sA[:, 0], gA[:, 0]
    k1, k2 = int(torch.argmin(sB)), int(torch.argmin(sA))
    if sB[k1] <= sA[k2]:
        d, a, n = sB[k1].item(), PA[k1].cpu().numpy(), gB[k1].cpu().numpy()
        pa, pb = a, a - d * n
    else:
        d, b, g = sA[k2].item(), PB[k2].cpu().numpy(), gA[k2].cpu().numpy()
        pa, pb, n = b - d * g, b, -g
    # boundaries crossing with no vertex inside (see C.true_distance): report contact
    L = max(torch.linalg.norm(torch.roll(P, -1, 0) - P, dim=1).max().item() for P in (PA, PB))
    if 0 < d <= L and C.polygons_cross(PA, PB):
        d = 0.0
    return d, pa, pb, n


def box_rows(N, dyn, v_max, w_max):
    """Input bounds; v_max / w_max may be per-robot arrays (0 fixes a robot)."""
    v_max, w_max = np.broadcast_to(v_max, N), np.broadcast_to(w_max, N)
    rows, rhs = [], []
    nu = 3 if dyn == 'si' else 2
    if dyn == 'si':
        ang = np.linspace(0, TWO_PI, 12, endpoint=False)
        nk = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    for i in range(N):
        if dyn == 'si':
            for n_ in nk:
                r = np.zeros(nu * N)
                r[3 * i:3 * i + 2] = -n_
                rows.append(r)
                rhs.append(-v_max[i] * np.cos(np.pi / 12))
        else:
            for s in (1, -1):
                r = np.zeros(nu * N)
                r[2 * i] = -s
                rows.append(r)
                rhs.append(-v_max[i])
        for s in (1, -1):
            r = np.zeros(nu * N)
            r[nu * i + nu - 1] = -s
            rows.append(r)
            rhs.append(-w_max[i])
    return np.array(rows), np.array(rhs)


def simulate_team(GT, shapes, fields, starts, goals, method, dyn='si', steps=1500, dt=0.01,
                  v_max=1.0, w_max=2.0, gamma=5.0, k_v=1.5, k_w=2.0, margin=0.02, heading=True,
                  record=False, stall_stop=None):
    """
    method: 'cspace' (ours), 'closest' (workspace closest-point barrier on the fitted shapes,
            idealised version of the conference method), 'circle' (bounding circles).
    dyn:    'si' (SE(2) single integrator, u_i = (vx, vy, omega)) or 'uni' (u_i = (v, omega)).
    """
    N = len(shapes)
    nu = 3 if dyn == 'si' else 2
    vm, wm = np.broadcast_to(v_max, N), np.broadcast_to(w_max, N)
    rho = np.array([s.rho for s in shapes])
    rho_gt = np.array([np.linalg.norm(G, axis=1).max() for G in GT])
    mec = [C.min_enclosing_circle(s.dense) for s in shapes] if method == 'circle' else None
    G_box, h_box = box_rows(N, dyn, v_max, w_max)
    x = starts.astype(float).copy()
    z_prev = None
    log = dict(min_gt=np.inf, min_h=np.inf, t_ctrl=[], slack_steps=0, path=np.zeros(N),
               min_row=np.inf, reached=None, traj=[], tv=0.0, pair_min_h=np.full((N, N), np.inf))
    u_prev = None
    for k in range(steps):
        t0 = time.perf_counter()
        # nominal
        u_nom = np.zeros(nu * N)
        for i in range(N):
            d = goals[i, :2] - x[i, :2]
            if dyn == 'si':
                v = k_v * d
                s = np.linalg.norm(v)
                u_nom[3 * i:3 * i + 2] = v if s <= vm[i] else v * vm[i] / max(s, 1e-12)
                u_nom[3 * i + 2] = np.clip(k_w * wrap(goals[i, 2] - x[i, 2]), -wm[i], wm[i]) if heading else 0.0
            else:
                r_ = np.linalg.norm(d)
                if r_ > 0.03:
                    al = wrap(np.arctan2(d[1], d[0]) - x[i, 2])
                    u_nom[2 * i] = np.clip(k_v * r_ * np.cos(al), -vm[i], vm[i])
                    u_nom[2 * i + 1] = np.clip(k_w * al, -wm[i], wm[i])
        # pairwise barrier rows, in (vx, vy, omega) coordinates of every robot
        rows3, hs = [], []
        sA, sT, sC, sE, n_box = [], [], [], [], 0
        for i, j in itertools.combinations(range(N), 2):
            pij = x[i, :2] - x[j, :2]
            r3 = np.zeros(3 * N)
            if method == 'summed':
                # summed-field barriers (prototype_3d/summed.py): Bernstein coefficient rows A u + T a + C >= 0,
                # A as (n, 6) blocks in (v_i, w_i, v_j, w_j), auxiliary a >= |E u| (relative twist of the pair)
                res = fields[(i, j)](x[i], x[j])
                if res is not None:
                    A6, Tk, C6, nb, hl, E6 = res
                    A3, E3 = np.zeros((len(C6), 3 * N)), np.zeros((len(E6), 3 * N))
                    A3[:, 3 * i:3 * i + 3], A3[:, 3 * j:3 * j + 3] = A6[:, :3], A6[:, 3:]
                    E3[:, 3 * i:3 * i + 3], E3[:, 3 * j:3 * j + 3] = E6[:, :3], E6[:, 3:]
                    sA.append(A3); sT.append(Tk); sC.append(C6); sE.append(E3); n_box += nb
                    hs.append(hl)
                    log["pair_min_h"][i, j] = min(log["pair_min_h"][i, j], hl)
                continue
            if method == 'cover':
                # surface-cover barriers (prototype_3d/dock_cover.py): many corner rows per pair, given as
                # (n, 6) rows in (v_i, w_i, v_j, w_j) and values h (n,)
                R6, hv = fields[(i, j)](x[i], x[j])
                if len(hv):
                    R3 = np.zeros((len(hv), 3 * N))
                    R3[:, 3 * i:3 * i + 3], R3[:, 3 * j:3 * j + 3] = R6[:, :3], R6[:, 3:]
                    rows3.extend(list(R3)); hs.extend(list(hv))
                    log["pair_min_h"][i, j] = min(log["pair_min_h"][i, j], float(hv.min()))
                continue
            if method == 'cspace':
                f, D, l_star = fields[(i, j)]
                t = rot(-x[j, 2]) @ pij
                if np.any(np.abs(t) >= D - C.ACT_DELTA):
                    continue
                th = (x[i, 2] - x[j, 2]) % TWO_PI
                phi, g = f.eval(np.r_[t, th], grad=True)
                h = phi[0] - l_star
                log["pair_min_h"][i, j] = min(log["pair_min_h"][i, j], h)
                gt, gth = g[0, :2], g[0, 2]
                r3[3 * i:3 * i + 2] = rot(x[j, 2]) @ gt
                r3[3 * j:3 * j + 2] = -rot(x[j, 2]) @ gt
                r3[3 * i + 2] = gth
                r3[3 * j + 2] = -gt @ (J2 @ t) - gth
            elif method == 'circle':
                # minimum enclosing circles, centres fixed in the body frames (they rotate with the body)
                ci, cj = rot(x[i, 2]) @ mec[i][0], rot(x[j, 2]) @ mec[j][0]
                d = x[i, :2] + ci - x[j, :2] - cj
                R = mec[i][1] + mec[j][1]
                if np.linalg.norm(d) > R + 0.3:
                    continue
                h = d @ d - (R + margin) ** 2
                r3[3 * i:3 * i + 2] = 2 * d
                r3[3 * i + 2] = 2 * d @ (J2 @ ci)
                r3[3 * j:3 * j + 2] = -2 * d
                r3[3 * j + 2] = -2 * d @ (J2 @ cj)
            else:
                if np.linalg.norm(pij) > rho[i] + rho[j] + 0.3:
                    continue
                PA = torch.as_tensor(shapes[i].world(x[i]), device=DEV)
                PB = torch.as_tensor(shapes[j].world(x[j]), device=DEV)
                d, pa, pb, n = closest_pair(PA, PB)
                h = d - margin
                r3[3 * i:3 * i + 2] = n
                r3[3 * i + 2] = n @ (J2 @ (pa - x[i, :2]))
                r3[3 * j:3 * j + 2] = -n
                r3[3 * j + 2] = -n @ (J2 @ (pb - x[j, :2]))
            rows3.append(r3)
            hs.append(h)
        rows = []
        for r3 in rows3:
            if dyn == 'si':
                rows.append(r3)
            else:
                r = np.zeros(2 * N)
                for m in range(N):
                    r[2 * m] = r3[3 * m:3 * m + 2] @ np.array([np.cos(x[m, 2]), np.sin(x[m, 2])])
                    r[2 * m + 1] = r3[3 * m + 2]
                rows.append(r)
        if method == 'summed':
            import summed as SM
            nu_ = 3 * N
            Call = np.concatenate(sC) if sC else np.zeros(0)
            Aall = np.vstack(sA) if sA else np.zeros((0, nu_))
            Eall = np.vstack(sE) if sE else np.zeros((0, nu_))
            Tall = np.zeros((len(Call), len(Eall)))          # one block of auxiliary variables per pair
            i_ = j_ = 0
            for T_, E_ in zip(sT, sE):
                Tall[i_:i_ + len(T_), j_:j_ + len(E_)] = T_
                i_ += len(T_); j_ += len(E_)
            # robots with zero speed bounds are fixed: eliminate their inputs (u = 0) and the rows that bound them
            free = np.repeat((vm > 0) & (wm > 0), 3)
            gfree = ~np.any(G_box[:, ~free] != 0, axis=1)
            uf, slack, z_prev = SM.solve(u_nom[free], Aall[:, free], Tall, Call, Eall[:, free], G_box[gfree][:, free],
                                         h_box[gfree], z_prev)
            u = np.zeros(nu_)
            u[free] = uf
            log.setdefault('n_rows', []).append(len(Call)); log.setdefault('n_box', []).append(n_box)
        else:
            G = np.vstack(rows + [G_box]) if rows else G_box
            hh = np.r_[[-gamma * h for h in hs], h_box]
            u, slack = solve_ldp_qp(u_nom, G, hh, len(rows), z_prev)
            z_prev = u
        log['t_ctrl'].append(time.perf_counter() - t0)
        if slack > 0:
            log['slack_steps'] += 1
        if hs:
            log['min_h'] = min(log['min_h'], min(hs))
            for r, h in zip(rows, hs):
                if h < 0.05:
                    log['min_row'] = min(log['min_row'], np.linalg.norm(r))
        # ground-truth distance (exact) for pairs whose bounding circles overlap
        d_step = np.inf
        for i, j in itertools.combinations(range(N), 2):
            dc = np.linalg.norm(x[i, :2] - x[j, :2])
            if dc < rho_gt[i] + rho_gt[j] + 0.1:          # exact whenever the true distance may be < 0.1 m
                wi = GT[i] @ rot(x[i, 2]).T + x[i, :2]
                wj = GT[j] @ rot(x[j, 2]).T + x[j, :2]
                d_step = min(d_step, C.true_distance(wi, wj))
            else:
                d_step = min(d_step, dc - rho_gt[i] - rho_gt[j])
        log['min_gt'] = min(log['min_gt'], d_step)
        if record:
            log['traj'].append(x.copy())
            log.setdefault('d_t', []).append(d_step)
            log.setdefault('h_t', []).append(min(hs) if hs else np.nan)
            log.setdefault('u_t', []).append(u.copy())
            log.setdefault('u_nom_t', []).append(u_nom.copy())
        if u_prev is not None:
            log['tv'] += np.abs(u - u_prev).sum()
        u_prev = u.copy()
        u = u.reshape(N, nu)
        x_old = x.copy()
        if dyn == 'si':
            x = x + dt * u
        else:
            x[:, 0] += dt * u[:, 0] * np.cos(x[:, 2])
            x[:, 1] += dt * u[:, 0] * np.sin(x[:, 2])
            x[:, 2] += dt * u[:, 1]
        log['path'] += np.linalg.norm(x[:, :2] - x_old[:, :2], axis=1)
        log['min_gt'] = min(log['min_gt'], step_gap(GT, rho_gt, x_old, x))
        pos_ok = np.linalg.norm(x[:, :2] - goals[:, :2], axis=1) < 0.05
        ang_ok = np.abs(wrap(x[:, 2] - goals[:, 2])) < 0.05 if (dyn == 'si' and heading) else True
        if np.all(pos_ok & ang_ok):
            log['reached'] = (k + 1) * dt
            break
        if stall_stop is not None:                  # stop early once every robot has come to rest
            log['_still'] = log.get('_still', 0) + 1 if np.linalg.norm(x - x_old) < 1e-4 else 0
            if log['_still'] >= stall_stop:
                break
    log['final'] = x
    log['t_ctrl'] = np.array(log['t_ctrl'])
    return log


rot, wrap = C.rot, C.wrap


step_gap = C.step_gap


def sample_poses(rng, rho, half, extra):
    N = len(rho)
    while True:
        P = rng.uniform(-half, half, (N, 2))
        if all(np.linalg.norm(P[i] - P[j]) >= rho[i] + rho[j] + extra
               for i, j in itertools.combinations(range(N), 2)):
            return np.column_stack([P, rng.uniform(-np.pi, np.pi, N)])


def exp_stats(n_trials=20, steps=1500):
    GT, shapes, fields, _ = load_cache()
    rho = np.array([s.rho for s in shapes])
    rng = np.random.default_rng(11)
    trials = [(sample_poses(rng, rho, 3.0, 0.5), sample_poses(rng, rho, 3.0, 0.5)) for _ in range(n_trials)]
    out = {}
    for method in ('cspace', 'closest', 'circle'):
        rec = []
        t0 = time.perf_counter()
        for k, (s, g) in enumerate(trials):
            L = simulate_team(GT, shapes, fields, s, g, method, steps=steps)
            rec.append(dict(min_gt=L['min_gt'], min_h=L['min_h'], reached=L['reached'],
                            path=L['path'].sum(), t_med=np.median(L['t_ctrl']) * 1e3,
                            t_max=L['t_ctrl'].max() * 1e3, slack_steps=L['slack_steps']))
        reached = [r['reached'] for r in rec if r['reached'] is not None]
        succ = [i for i, r in enumerate(rec) if r['reached'] is not None]
        out[method] = dict(
            trials=rec,
            collisions=int(sum(r['min_gt'] <= 0 for r in rec)),
            success_collision_free=int(sum(r['reached'] is not None and r['min_gt'] > 0 for r in rec)),
            min_gt=min(r['min_gt'] for r in rec),
            success=len(reached),
            mean_time=float(np.mean(reached)) if reached else None,
            mean_path_success=float(np.mean([rec[i]['path'] for i in succ])) if succ else None,
            t_med=float(np.median([r['t_med'] for r in rec])),
            t_max=float(max(r['t_max'] for r in rec)),
            slack_steps=int(sum(r['slack_steps'] for r in rec)))
        print(f"[{method}] collisions {out[method]['collisions']}/{n_trials} | min GT dist {out[method]['min_gt']:+.4f} | "
              f"reached {len(reached)}/{n_trials} | mean time {out[method]['mean_time']} | "
              f"ctrl median {out[method]['t_med']:.2f} ms max {out[method]['t_max']:.1f} ms | "
              f"infeasible steps {out[method]['slack_steps']} | {time.perf_counter() - t0:.0f} s")
    # paired comparison on trials that every method solved
    common = [k for k in range(n_trials) if all(out[m]['trials'][k]['reached'] is not None
                                                and out[m]['trials'][k]['min_gt'] > 0 for m in out)]
    out['common_trials'] = common
    for m in ('cspace', 'closest', 'circle'):
        out[m]['common_mean_time'] = float(np.mean([out[m]['trials'][k]['reached'] for k in common])) if common else None
        out[m]['common_mean_path'] = float(np.mean([out[m]['trials'][k]['path'] for k in common])) if common else None
    print('common trials:', len(common), {m: (out[m]['common_mean_time'], out[m]['common_mean_path']) for m in ('cspace', 'closest', 'circle')})
    dump('stats', out)


# =============================================================================================
# unicycle
# =============================================================================================
def plot_traj(GT, shapes, traj, goals, fname, title, lim, obstacles=None, every=60):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    N = len(shapes)
    cols = plt.cm.tab10(np.arange(N))
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    if obstacles is not None:
        for Gw, Fw in obstacles:
            ax.add_patch(Polygon(Gw, closed=True, color='0.55'))
            ax.add_patch(Polygon(Fw, closed=True, fill=False, ec='k', lw=1.2))
    traj = np.array(traj)
    for i in range(N):
        ax.plot(traj[:, i, 0], traj[:, i, 1], color=cols[i], lw=0.9)
        for k in list(range(0, len(traj), every)) + [len(traj) - 1]:
            ax.add_patch(Polygon(shapes[i].world(traj[k, i]), closed=True, fill=False, ec=cols[i],
                                 lw=0.8, alpha=0.35 if k < len(traj) - 1 else 1.0))
        xf = traj[-1, i]
        ax.add_patch(Polygon(GT[i] @ rot(xf[2]).T + xf[:2], closed=True, color='0.55'))
        ax.plot(goals[i, 0], goals[i, 1], 'x', color=cols[i], ms=7)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, fname), dpi=150)
    plt.close(fig)


def ab_field():
    """Certified field for the C-shape / three-lobed pair (cached)."""
    A, B = C.shape_A(), C.shape_B()
    fn = os.path.join(CACHE, 'AB.npz')
    if os.path.exists(fn):
        d = np.load(fn)
        f = C.Field3D('bspline', float(d['D']), int(d['Kxy']), int(d['Kth']))
        f.set_W(d['W'])
        return A, B, f, float(d['D']), float(d['l'])
    f, D, l_star, ok, _ = pair_field(A, B)
    assert ok
    os.makedirs(CACHE, exist_ok=True)
    np.savez(fn, W=f.W, D=D, Kxy=f.bx.K, Kth=f.bt.K, l=l_star)
    return A, B, f, D, l_star


def exp_dock(steps=700):
    """Docking into the C-shape's concavity: ours vs closest-point vs bounding circles."""
    A, B, f, D, l_star = ab_field()
    GT2 = [A.dense, B.dense]
    starts = np.array([[0.0, 0.0, 0.0], [2.6, 0.25, -0.6]])
    goals = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, np.pi / 3]])        # goal itself is in collision
    out = {}
    for method in ('cspace', 'closest', 'circle'):
        L = simulate_team(GT2, [A, B], {(0, 1): (f, D, l_star)}, starts, goals, method, steps=steps,
                          v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        xB = L['final'][1]
        out[method] = dict(final_xB=xB[0], final_gt=C.true_distance(A.world(L['final'][0]), B.world(xB)),
                           min_gt=L['min_gt'], min_h=L['min_h'], tv=L['tv'],
                           t_med=np.median(L['t_ctrl']) * 1e3)
        print(f"dock [{method}]: {out[method]}", flush=True)
    dump('dock', out)


def exp_unicycle():
    out = {}
    A, B, f, D, l_star = ab_field()
    GT2 = [A.dense, B.dense]
    fig_two = fig_four = None

    def summary(L, goals):
        return dict(min_gt=L['min_gt'], min_h=L['min_h'], reached=L['reached'], min_row=L['min_row'],
                    interacted=bool(L['min_h'] < 0.05),
                    final_err=np.linalg.norm(L['final'][:, :2] - goals[:, :2], axis=1),
                    t_med=np.median(L['t_ctrl']) * 1e3, slack_steps=L['slack_steps'])

    # (a) two robots swapping along a line, sweep of the lateral offset
    out['two'] = {}
    for dy in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        starts = np.array([[-2.4, dy, 0.0], [2.4, -dy, np.pi]])
        goals = np.array([[2.4, dy, 0.0], [-2.4, -dy, np.pi]])
        L = simulate_team(GT2, [A, B], {(0, 1): (f, D, l_star)}, starts, goals, 'cspace', dyn='uni',
                          steps=2000, w_max=np.pi / 2, k_v=1.0, k_w=2.0, record=True)
        out['two'][dy] = summary(L, goals)
        print(f"unicycle two dy={dy}: {out['two'][dy]}", flush=True)
        if fig_two is None and L['reached'] and out['two'][dy]['interacted']:
            fig_two = dy
            plot_traj(GT2, [A, B], L['traj'], goals, 'unicycle_two.png',
                      f"unicycle, offset {dy} m: min dist {L['min_gt']:.3f} m, goals at {L['reached']:.1f} s", 3.0)
    # (b) four heterogeneous robots (cached shapes 0-3) moving to the next corner, sweep of the square
    GT, shapes, fields, _ = load_cache()
    GT4, sh4 = GT[:4], shapes[:4]
    f4 = {k: v for k, v in fields.items() if max(k) < 4}
    out['four'] = {}
    for mode, perm in (('cyclic', [1, 2, 3, 0]), ('diagonal', [2, 3, 0, 1])):
        for c in (0.9, 1.2, 1.6, 2.2):
            corners = np.array([[-c, -c], [c, -c], [c, c], [-c, c]])
            head = [np.arctan2(*(corners[p] - corners[i])[::-1]) for i, p in enumerate(perm)]
            starts = np.column_stack([corners, head])
            goals = np.column_stack([corners[perm], np.zeros(4)])
            if min(np.linalg.norm(corners[i] - corners[j]) - shapes[i].rho - shapes[j].rho
                   for i, j in itertools.combinations(range(4), 2)) < 0.15:
                continue
            L = simulate_team(GT4, sh4, f4, starts, goals, 'cspace', dyn='uni', steps=3000,
                              w_max=np.pi / 2, k_v=1.0, k_w=2.0, record=True)
            out['four'][f'{mode}_{c}'] = summary(L, goals)
            print(f"unicycle four {mode} c={c}: {out['four'][f'{mode}_{c}']}", flush=True)
            if fig_four is None and L['reached'] and out['four'][f'{mode}_{c}']['interacted']:
                fig_four = (mode, c)
                plot_traj(GT4, sh4, L['traj'], goals, 'unicycle_four.png',
                          f"unicycle, 4 robots ({mode}, half-side {c} m): min dist {L['min_gt']:.3f} m", c + 1.0)
    out['figures'] = dict(two=fig_two, four=fig_four)
    dump('unicycle', out)


# =============================================================================================
# static obstacles: one field for a robot among m obstacles
# =============================================================================================
class Union:
    def __init__(self, curves):
        self.bez = np.concatenate([c.bez for c in curves])
        self.dense_list = [c.dense for c in curves]


def union_cspace_slices(A, polys, thetas, D, res=0.02):
    n = 2 * int(round(D / res)) + 1
    g = (np.arange(n) - n // 2) * res
    GX, GY = np.meshgrid(g, g, indexing='ij')
    P = np.stack([GX.ravel(), GY.ravel()], axis=1)
    mB = np.zeros(n * n, dtype=bool)
    for poly in polys:
        mB |= Path(poly).contains_points(P)
    mB = mB.reshape(n, n).astype(float)
    out = np.empty((len(thetas), n, n), dtype=np.float32)
    for k, th in enumerate(thetas):
        mA = Path(-(A.dense @ rot(th).T)).contains_points(P).reshape(n, n).astype(float)
        co = fftconvolve(mB, mA, mode='same') > 0.5
        out[k] = (ndimage.distance_transform_edt(~co) - ndimage.distance_transform_edt(co)) * res
    return g, out


def exp_static(ms=(1, 3, 5, 7), W=3.8):
    GT, shapes, _, _ = load_cache()
    robot, robot_gt = shapes[0], GT[0]
    rng = np.random.default_rng(7)
    cells = [(-2.2, 0), (-2.2, 2.2), (0, -2.2), (0, 0), (0, 2.2), (2.2, -2.2), (2.2, 0)]
    obs_all = []
    for cxy in cells:
        G = random_shape(rng)
        S, _ = enclose(G)
        c = np.array(cxy) + rng.uniform(-0.15, 0.15, 2)
        obs_all.append((G + c, S, c))
    order = rng.permutation(len(cells))
    # waypoints: horizontal legs pass the middle-row obstacles at 0.6 m (closer than the robot's
    # clearance needs, so the barrier is active), vertical legs run between columns; rotating
    wps = np.array([[-3.3, -0.6, 0.0], [-1.1, -0.6, np.pi / 2], [-1.1, 0.6, np.pi],
                    [1.1, 0.6, 3 * np.pi / 2], [1.1, -0.6, 2 * np.pi], [3.3, -0.6, 5 * np.pi / 2]])
    out = {}
    for m in ms:
        obs = [obs_all[k] for k in order[:m]]
        fitted_world = [C.SmoothShape(S.ctrl + c) for _, S, c in obs]           # world-frame curves
        t0 = time.perf_counter()
        th = np.arange(96) * TWO_PI / 96
        g, Sl = union_cspace_slices(robot, [F.dense for F in fitted_world], th, W + 0.1)
        n_xy = int(round(2 * W / 0.037))
        gxy = np.linspace(-W, W, n_xy)
        f = C.Field3D('bspline', W, int(round(2 * W / 0.0625)), 48)
        f.fit(C.sample_slices(g, Sl.astype(np.float64), gxy, gxy).transpose(1, 2, 0), gxy, gxy, th)
        t_fit = time.perf_counter() - t0
        cert = C.Certifier(f, robot, Union(fitted_world))
        l_star, ok, n_reg, t_cert = cert.certify(verbose=False)
        band_ok, band_min = cert.certify_band(l_star)
        reg_ok, _, _ = cert.certify_regular(l_star)
        print(f"[m={m}] {f.bx.K}x{f.by.K}x{f.bt.K} cells | fit {t_fit:.1f} s | l* {l_star:+.4f} cert {ok} "
              f"({t_cert:.1f} s) band {band_ok} ({band_min:+.3f}) regular {reg_ok}", flush=True)
        res = dict(cells=[f.bx.K, f.by.K, f.bt.K], t_fit=t_fit, l_star=l_star, certified=ok,
                   t_cert=t_cert, band=band_ok, regular=reg_ok)
        for method in ('cspace', 'closest'):
            L = simulate_static(robot, robot_gt, f, l_star, fitted_world, [o[0] for o in obs], wps,
                                method, W, record=(m == ms[-1]))
            res[method] = dict(min_gt=L['min_gt'], min_h=L['min_h'], reached=L['reached'],
                               t_med=np.median(L['t_ctrl']) * 1e3, t_max=L['t_ctrl'].max() * 1e3,
                               n_constraints=L['n_con'], slack_steps=L['slack_steps'],
                               final_err=float(np.linalg.norm(L['final'][:2] - wps[-1, :2])))
            print(f"   {method}: {res[method]}", flush=True)
            if method == 'cspace' and m == ms[-1]:
                plot_traj([robot_gt], [robot], [[p] for p in L['traj']], wps[-1:], f'static_m{m}.png',
                          f"one field for {m} obstacles: min dist {L['min_gt']:.3f} m", W,
                          obstacles=[(o[0], F.dense) for o, F in zip(obs, fitted_world)], every=60)
        out[m] = res
        del f, cert
        torch.cuda.empty_cache()
    dump('static', out)


def simulate_static(robot, robot_gt, f, l_star, fitted_world, gt_world, wps, method, W,
                    steps=4000, dt=0.01, v_max=1.0, w_max=2.0, gamma=5.0, k_v=1.5, k_w=2.0,
                    margin=0.02, r_switch=0.3, record=False):
    G_box, h_box = box_rows(1, 'si', v_max, w_max)
    x = wps[0].copy()
    wi = 1
    log = dict(min_gt=np.inf, min_h=np.inf, t_ctrl=[], slack_steps=0, reached=None, traj=[], n_con=0)
    obs_t = [torch.as_tensor(F.dense, device=DEV) for F in fitted_world]
    obs_c = [F.dense.mean(0) for F in fitted_world]
    obs_r = [np.linalg.norm(F.dense - F.dense.mean(0), axis=1).max() for F in fitted_world]
    # ground truth as bodies for step_gap: robot (index 0) and the static obstacles about their centers
    gc = [Gw.mean(0) for Gw in gt_world]
    GT_all = [robot_gt] + [Gw - c for Gw, c in zip(gt_world, gc)]
    rho_all = np.array([np.linalg.norm(G, axis=1).max() for G in GT_all])
    X_obs = np.array([[c[0], c[1], 0.0] for c in gc]).reshape(-1, 3)
    pairs = [(0, k + 1) for k in range(len(gt_world))]
    z_prev = None
    for k in range(steps):
        t0 = time.perf_counter()
        goal = wps[wi]
        if wi < len(wps) - 1 and np.linalg.norm(goal[:2] - x[:2]) < r_switch:
            wi += 1
            goal = wps[wi]
        d = goal[:2] - x[:2]
        v = k_v * d if wi == len(wps) - 1 else d / max(np.linalg.norm(d), 1e-9) * v_max
        s = np.linalg.norm(v)
        u_nom = np.r_[v if s <= v_max else v * v_max / s, np.clip(k_w * wrap(goal[2] - x[2]), -w_max, w_max)]
        rows, hs = [], []
        if method == 'cspace':
            if np.all(np.abs(x[:2]) < W - C.ACT_DELTA):
                phi, g = f.eval(np.r_[x[:2], x[2] % TWO_PI], grad=True)
                rows.append(g[0])
                hs.append(phi[0] - l_star)
        else:
            PA = torch.as_tensor(robot.world(x), device=DEV)
            for Ot, oc, orr in zip(obs_t, obs_c, obs_r):
                if np.linalg.norm(x[:2] - oc) > robot.rho + orr + 0.3:
                    continue
                dd, pa, _, n = closest_pair(PA, Ot)
                rows.append(np.r_[n, n @ (J2 @ (pa - x[:2]))])
                hs.append(dd - margin)
        log['n_con'] = max(log['n_con'], len(rows))
        G = np.vstack(rows + [G_box]) if rows else G_box
        hh = np.r_[[-gamma * h for h in hs], h_box]
        u, slack = solve_ldp_qp(u_nom, G, hh, len(rows), z_prev)
        z_prev = u
        log['t_ctrl'].append(time.perf_counter() - t0)
        log['slack_steps'] += int(slack > 0)
        if hs:
            log['min_h'] = min(log['min_h'], min(hs))
        rw = robot_gt @ rot(x[2]).T + x[:2]
        for Gw in gt_world:
            if np.linalg.norm(x[:2] - Gw.mean(0)) < 2.0:
                log['min_gt'] = min(log['min_gt'], C.true_distance(rw, Gw))
        if record:
            log['traj'].append(x.copy())
        x_old = x
        x = x + dt * u
        log['min_gt'] = min(log['min_gt'], step_gap(GT_all, rho_all, np.vstack([x_old, X_obs]),
                                                    np.vstack([x, X_obs]), pairs))
        if wi == len(wps) - 1 and np.linalg.norm(x[:2] - goal[:2]) < 0.05 and abs(wrap(x[2] - goal[2])) < 0.05:
            log['reached'] = (k + 1) * dt
            break
    log['final'] = x
    log['t_ctrl'] = np.array(log['t_ctrl'])
    return log


# =============================================================================================
# sensitivity: resolution, representation, training data, tolerance; per-cell levels
# =============================================================================================
def exp_sens():
    A, B = C.shape_A(), C.shape_B()
    D = round(A.rho + B.rho + 0.3, 2)
    rng = np.random.default_rng(0)
    th_test = (np.arange(16) + 0.37) * TWO_PI / 16
    g1, S_test = C.cspace_sdf_slices(A, B, th_test, D, res=0.01)
    Xg, Yg = np.meshgrid(g1, g1, indexing='ij')
    tq, ts = [], []
    for k, th in enumerate(th_test):
        sel = rng.choice(Xg.size, 4000, replace=False)
        tq.append(np.stack([Xg.ravel()[sel], Yg.ravel()[sel], np.full(4000, th)], axis=1))
        ts.append(S_test[k].ravel()[sel])
    tq, ts = np.concatenate(tq), np.concatenate(ts)
    near = np.abs(ts) < 0.1

    def training(res, n_th, n_xy):
        th = np.arange(n_th) * TWO_PI / n_th
        g, S = C.cspace_sdf_slices(A, B, th, D, res=res)
        gxy = np.linspace(-D, D, n_xy)
        return C.sample_slices(g, S, gxy, gxy).transpose(1, 2, 0), gxy, th

    def run(tag, kind, Kxy, Kth, data, tol=2e-3):
        T, gxy, th = data
        f = C.Field3D(kind, D, Kxy, Kth)
        f.fit(T, gxy, gxy, th)
        err = f.eval(tq) - ts
        l_star, ok, n_reg, t_cert = C.Certifier(f, A, B).certify(tol=tol, verbose=False)
        r = dict(tag=tag, kind=kind, Kxy=Kxy, Kth=Kth, h=2 * D / Kxy, n_coef=f.n_coef, tol=tol,
                 rmse_bd=float(np.sqrt(np.mean(err[near] ** 2))), max_bd=float(np.abs(err[near]).max()),
                 l_star=l_star, certified=ok, t_cert=t_cert, n_reg=n_reg)
        print(f"  {tag:10s} {kind:8s} {Kxy:3d}x{Kxy:3d}x{Kth:3d} coef {f.n_coef:7d} | rmse_bd {r['rmse_bd']:.4f} "
              f"max_bd {r['max_bd']:.4f} | l* {l_star:+.4f} ({ok}) {t_cert:5.1f} s", flush=True)
        return r, f

    rows = []
    fine = training(0.01, 192, 120)                      # rich data for the model-size sweeps
    print('resolution sweeps (training: 1 cm raster, 192 angles, 120^2 grid)', flush=True)
    for h in (0.125, 0.0833, 0.0625, 0.047):
        rows.append(run('space', 'bspline', int(round(2 * D / h)), 48, fine)[0])
    for Kth in (24, 96):
        rows.append(run('theta', 'bspline', 45, Kth, fine)[0])
    for K in (8, 12, 16, 20):
        rows.append(run('bezier', 'bezier', K, K, fine)[0])
    print('training-data resolution (B-spline 45x45x48, 96 angles, 100^2 grid)', flush=True)
    base = None
    for res in (0.02, 0.01, 0.005):
        r, f = run(f'raster{res * 100:g}cm', 'bspline', 45, 48, training(res, 96, 100))
        rows.append(r)
        if res == 0.01:
            base = f
    print('certificate tolerance (base field)', flush=True)
    for tol in (5e-4, 1e-3, 4e-3):
        l_star, ok, n_reg, t_cert = C.Certifier(base, A, B).certify(tol=tol, verbose=False)
        rows.append(dict(tag='tol', kind='bspline', Kxy=45, Kth=48, tol=tol, l_star=l_star,
                         certified=ok, t_cert=t_cert, n_reg=n_reg, n_coef=base.n_coef))
        print(f"  tol {tol * 1e3:.1f} mm: l* {l_star:+.4f} ({ok}) {t_cert:.1f} s, {n_reg} regions", flush=True)

    print('per-cell levels (base field)', flush=True)
    lc, ok, t_cells = C.Certifier(base, A, B).certify_cells(tol=2e-3)
    has = np.isfinite(lc)
    vals = lc[has]
    l_glob = vals.max()
    percell = dict(certified=ok, time=t_cells, n_cells_total=int(lc.size), n_cells_M=int(has.sum()),
                   l_global=float(l_glob), quantiles={q: float(np.quantile(vals, q)) for q in (0.5, 0.9, 0.99)},
                   frac_above_half=float(np.mean(vals > l_glob / 2)), frac_below_zero=float(np.mean(vals < 0)),
                   mean_slack=float(np.mean(l_glob - vals)))
    print('  ', percell, flush=True)
    dump('sens', dict(rows=rows, percell=percell))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for tag, col, mk, lab in (('space', 'tab:blue', 'o', 'B-spline, cell size'),
                              ('theta', 'tab:cyan', 's', r'B-spline, $\theta$ cells'),
                              ('bezier', 'tab:orange', '^', r'piecewise B$\'e$zier')):
        rr = sorted([r for r in rows if r['tag'] == tag] +
                    ([r for r in rows if r['tag'] == 'space' and r['Kxy'] == 45] if tag == 'theta' else []),
                    key=lambda r: r['n_coef'])
        axs[0].plot([r['n_coef'] for r in rr], [r['l_star'] for r in rr], mk + '-', color=col, label=lab)
        axs[1].plot([r['n_coef'] for r in rr], [r['max_bd'] for r in rr], mk + '-', color=col, label=lab)
    for ax, yl in zip(axs, (r'certified $l^\star$ [m]', 'max. error near boundary [m]')):
        ax.set_xscale('log')
        ax.set_xlabel('number of coefficients')
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'sens_resolution.png'), dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    axs[0].hist(vals, bins=80, color='tab:blue')
    axs[0].axvline(l_glob, color='k', ls='--', lw=1)
    axs[0].set_yscale('log')
    axs[0].set_xlabel(r'per-cell certified level $\ell_c$ [m]')
    axs[0].set_ylabel('cells')
    mp = np.where(has, lc, -np.inf).max(axis=2)
    xs = base.bx.lo + (np.arange(base.bx.K) + 0.5) * base.bx.h
    im = axs[1].pcolormesh(xs, xs, np.where(np.isfinite(mp), mp, np.nan).T, cmap='viridis', shading='auto')
    fig.colorbar(im, ax=axs[1], label=r'$\max_\theta \ell_c$ [m]')
    axs[1].set_aspect('equal')
    axs[1].set_xlabel(r'$t_x$ [m]')
    axs[1].set_ylabel(r'$t_y$ [m]')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'percell.png'), dpi=150)
    plt.close(fig)
    print('>>> saved tro/figs/sens_resolution.png, tro/figs/percell.png')


# =============================================================================================
# docking with a long probe into a deep slot (replaces the C-shape / three-lobed docking)
# =============================================================================================
DOCK2_K_CURVE = 192      # control points of each certified boundary curve
DOCK2_DENSE = 1920       # boundary samples (polygon) of each fitted body
DOCK2_RASTER = 0.0025    # raster resolution of the training targets [m]
DOCK2_N_TRAIN = 576      # training orientations
DOCK2_K_TH = 384         # angular cells of the B-spline field


def dock2_setup():
    """Certified boundaries and fields (B-spline and piecewise Bezier) for the station/craft pair."""
    from dock_shapes import craft_body, station_body
    from cspace_cbf_5robots import certify_curve_encloses, fit_curve
    fn = os.path.join(CACHE, 'dock_craft.npz')
    GA, GB = station_body(), craft_body()
    if os.path.exists(fn):
        d = np.load(fn)
        A, B = C.SmoothShape(d['ctrlA'], n_dense=DOCK2_DENSE), C.SmoothShape(d['ctrlB'], n_dense=DOCK2_DENSE)
        fields = {}
        for kind in ('bspline', 'bezier'):
            f = C.Field3D(kind, float(d['D']), int(d[f'Kxy_{kind}']), int(d[f'Kth_{kind}']))
            f.set_W(d[f'W_{kind}'])
            fields[kind] = (f, float(d['D']), float(d[f'l_{kind}']))
        return GA, GB, A, B, fields, json.loads(str(d['meta']))
    meta = {}
    ctrl = []
    for name, G in (('A', GA), ('B', GB)):
        for off in np.arange(0.001, 0.06, 0.001):
            S, ok = certify_curve_encloses(fit_curve(G, DOCK2_K_CURVE, off, n_fit=4 * DOCK2_K_CURVE), G)
            if ok:
                break
        assert ok
        ctrl.append(S.ctrl)
        meta[f'offset_{name}'] = float(off)
    A, B = C.SmoothShape(ctrl[0], n_dense=DOCK2_DENSE), C.SmoothShape(ctrl[1], n_dense=DOCK2_DENSE)
    D = round(A.rho + B.rho + 0.3, 2)
    t0 = time.perf_counter()
    th = np.arange(DOCK2_N_TRAIN) * TWO_PI / DOCK2_N_TRAIN
    gxy = np.linspace(-D, D, int(round(2 * D / 0.02)))
    fn_T = os.path.join(CACHE, 'dock_craft_T.npz')          # training data (raster), reusable
    if os.path.exists(fn_T):
        T = np.load(fn_T)['T']
    else:
        T = C.cspace_sdf_samples(A, B, th, D, DOCK2_RASTER, gxy, gxy)
        np.savez(fn_T, T=T)
    meta['t_raster'] = time.perf_counter() - t0
    data = dict(ctrlA=ctrl[0], ctrlB=ctrl[1], D=D)
    fields = {}
    for kind, Kxy, Kth in (('bspline', int(round(2 * D / 0.0625)), DOCK2_K_TH), ('bezier', 30, 30)):
        f = C.Field3D(kind, D, Kxy, Kth)
        t0 = time.perf_counter()
        f.fit(T, gxy, gxy, th)
        t_fit = time.perf_counter() - t0
        cert = C.Certifier(f, A, B)
        l_star, ok, n_reg, t_cert = cert.certify(verbose=False)
        band_ok, band_min = cert.certify_band(l_star)
        reg_ok, _, _ = cert.certify_regular(l_star)
        meta[kind] = dict(n_coef=f.n_coef, cells=[f.bx.K, f.by.K, f.bt.K], l_star=l_star, certified=ok,
                          t_fit=t_fit, t_cert=t_cert, band=band_ok, band_min=band_min, regular=reg_ok)
        print(f"  dock2 {kind}: {meta[kind]}", flush=True)
        data.update({f'W_{kind}': f.W, f'Kxy_{kind}': Kxy, f'Kth_{kind}': Kth, f'l_{kind}': l_star})
        fields[kind] = (f, D, l_star)
    data['meta'] = json.dumps(meta)
    os.makedirs(CACHE, exist_ok=True)
    np.savez(fn, **data)
    return GA, GB, A, B, fields, meta


from dock_shapes import SEAT_X as _SEAT_X  # noqa: E402

DOCK2_START = np.array([[0.0, 0.0, 0.0], [3.4, 0.35, np.pi - 0.7]])
DOCK2_GOAL = np.array([[0.0, 0.0, 0.0], [_SEAT_X - 0.10, 0.0, np.pi]])   # 10 cm past full seating: infeasible


def exp_dock2(steps=900):
    GA, GB, A, B, fields, meta = dock2_setup()
    print('dock2 setup:', meta, flush=True)
    out = dict(setup=meta)
    # surface cover (prototype_3d/dock_cover.py): fields of the ground-truth polygons, tightest levels,
    # 5-mm boxes on the craft's certified curve
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'prototype_3d'))
    import franka3d as F3
    from dock_cover import body_field, tight_level_2d, cover_2d, make_cover_fn
    fS, fC = body_field(GA, .0125, .15), body_field(GB, .01, .06)
    MS, MC = F3.hessian_majorant(fS), F3.hessian_majorant(fC)
    lS, lC = tight_level_2d(fS, MS, GA), tight_level_2d(fC, MC, GB)
    PC, sz = cover_2d(fC, lC, .005)
    fields = dict(fields, cover=make_cover_fn(fS, MS, lS, PC, sz))
    out['cover_setup'] = dict(l_station_mm=1e3 * lS, l_craft_mm=1e3 * lC, n_boxes=len(PC), side_mm=1e3 * sz)
    runs = (('cover', 'cover', 'cover'), ('ours_bspline', 'cspace', 'bspline'), ('ours_bezier', 'cspace', 'bezier'),
            ('closest', 'closest', 'bspline'), ('circle', 'circle', 'bspline'))
    for name, method, kind in runs:
        L = simulate_team([GA, GB], [A, B], {(0, 1): fields[kind]}, DOCK2_START, DOCK2_GOAL, method,
                          steps=steps, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        xB = L['final'][1]
        u = np.array(L['u_t'])
        from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH
        tip = xB[:2] + NOSE_TIP_X * np.array([np.cos(xB[2]), np.sin(xB[2])])
        out[name] = dict(final_xB=xB[0], final_thB=xB[2], tip_depth=MOUTH_X - tip[0],
                         seat_shortfall=SEAT_DEPTH - (MOUTH_X - tip[0]),
                         final_gt=C.true_distance(GA, GB @ C.rot(xB[2]).T + xB[:2]), min_gt=L['min_gt'],
                         min_h=float(np.nanmin(L['h_t'])), tv=L['tv'],
                         tv_w=float(np.abs(np.diff(u[:, 5])).sum()),
                         t_med=np.median(L['t_ctrl']) * 1e3, t_p95=np.percentile(L['t_ctrl'], 95) * 1e3,
                         path=float(L['path'][1]))
        print(f"dock2 [{name}]: {out[name]}", flush=True)
    dump('dock2', out)


def exp_dock2_dt(T=9.0):
    """Docking with the same initial conditions and time steps of 10, 5, 2, and 1 ms."""
    GA, GB, A, B, fields, meta = dock2_setup()
    from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH
    out = {}
    for dt in (0.01, 0.005, 0.002, 0.001):
        for name, method, kind in (('ours_bspline', 'cspace', 'bspline'), ('ours_bezier', 'cspace', 'bezier'),
                                   ('closest', 'closest', 'bspline')):
            L = simulate_team([GA, GB], [A, B], {(0, 1): fields[kind]}, DOCK2_START, DOCK2_GOAL, method,
                              steps=int(round(T / dt)), dt=dt, v_max=np.array([0.0, 1.0]),
                              w_max=np.array([0.0, 2.0]), record=True)
            xB = L['final'][1]
            tip = xB[:2] + NOSE_TIP_X * np.array([np.cos(xB[2]), np.sin(xB[2])])
            u = np.array(L['u_t'])
            out[f'{name}_{dt}'] = dict(dt=dt, method=name, min_h=float(np.nanmin(L['h_t'])), min_gt=L['min_gt'],
                                       depth=MOUTH_X - tip[0], tv=L['tv'],
                                       tv_w=float(np.abs(np.diff(u[:, 5])).sum()))
            print(f"dock2 dt {dt}: [{name}] {out[f'{name}_{dt}']}", flush=True)
    dump('dock2_dt', out)


def exp_dock2_res():
    """Docking pair: certified level and insertion depth versus the angular resolution of the
    B-spline field, and a piecewise Bezier field with a matched number of coefficients."""
    GA, GB, A, B, fields, meta = dock2_setup()
    from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH
    D = fields['bspline'][1]
    T = np.load(os.path.join(CACHE, 'dock_craft_T.npz'))['T']
    th = np.arange(T.shape[2]) * TWO_PI / T.shape[2]
    gxy = np.linspace(-D, D, T.shape[0])
    K0 = int(round(2 * D / 0.0625))
    out = []
    for kind, Kxy, Kth in (('bspline', K0, 96), ('bspline', K0, 192), ('bspline', K0, 384),
                           ('bezier', 46, 46)):
        f = C.Field3D(kind, D, Kxy, Kth)
        t0 = time.perf_counter()
        f.fit(T, gxy, gxy, th)
        t_fit = time.perf_counter() - t0
        cert = C.Certifier(f, A, B)
        l_star, ok, _, t_contact = cert.certify(verbose=False)
        t0 = time.perf_counter()
        band, _ = cert.certify_band(l_star)
        t_band = time.perf_counter() - t0
        t0 = time.perf_counter()
        reg, _, _ = cert.certify_regular(l_star)
        t_reg = time.perf_counter() - t0
        L = simulate_team([GA, GB], [A, B], {(0, 1): (f, D, l_star)}, DOCK2_START, DOCK2_GOAL, 'cspace',
                          steps=900, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        xB = L['final'][1]
        tip = xB[:2] + NOSE_TIP_X * np.array([np.cos(xB[2]), np.sin(xB[2])])
        out.append(dict(kind=kind, Kxy=Kxy, Kth=Kth, n_coef=f.n_coef, l_star=l_star, certified=ok,
                        band=band, regular=reg, t_fit=t_fit, t_contact=t_contact, t_band=t_band,
                        t_regular=t_reg, depth=MOUTH_X - tip[0], shortfall=SEAT_DEPTH - (MOUTH_X - tip[0]),
                        min_h=float(np.nanmin(L['h_t'])), min_gt=L['min_gt'], tv=L['tv']))
        print('dock2 res:', out[-1], flush=True)
        del f, cert
    dump('dock2_res', out)


def exp_percell_fig():
    """Per-cell certified levels of the base B-spline field, shown on theta slices."""
    A, B = C.shape_A(), C.shape_B()
    D = round(A.rho + B.rho + 0.3, 2)
    th = np.arange(96) * TWO_PI / 96
    g, S = C.cspace_sdf_slices(A, B, th, D, res=0.01)
    gxy = np.linspace(-D, D, 100)
    f = C.Field3D('bspline', D, 45, 48)
    f.fit(C.sample_slices(g, S, gxy, gxy).transpose(1, 2, 0), gxy, gxy, th)
    cert = C.Certifier(f, A, B)
    l_star, _, _, _ = cert.certify(verbose=False)
    fn = os.path.join(CACHE, 'percell_lc.npy')
    if os.path.exists(fn):
        lc = np.load(fn)
    else:
        lc, ok, t = cert.certify_cells(tol=2e-3)
        print(f"per-cell certificate: {ok}, {t:.0f} s", flush=True)
        np.save(fn, lc)
    has = np.isfinite(lc)
    vals = lc[has]
    k_star = np.unravel_index(np.argmax(np.where(has, lc, -np.inf)), lc.shape)[2]
    slices = [k_star, (k_star + f.bt.K // 4) % f.bt.K]

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig = plt.figure(figsize=(12, 3.9))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05])
    ax0 = fig.add_subplot(gs[0])
    ax0.hist(vals, bins=80, color='tab:blue')
    ax0.axvline(l_star, color='r', ls='--', lw=1)
    ax0.set_yscale('log')
    ax0.set_xlabel(r'per-cell certified level $\ell_c$ [m]')
    ax0.set_ylabel('cells meeting $M$')
    ax0.set_title('(a) all cells', fontsize=10)
    edges = f.bx.lo + np.arange(f.bx.K + 1) * f.bx.h
    xs = np.linspace(-D, D, 300)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    cmap = plt.cm.viridis
    for n, kt in enumerate(slices):
        ax = fig.add_subplot(gs[1 + n])
        thc = (kt + 0.5) * f.bt.h
        layer = lc[:, :, kt]
        inner = np.where(np.isfinite(layer) & (layer <= 0), 1.0, np.nan)
        ax.pcolormesh(edges, edges, inner.T, cmap=ListedColormap(['0.85']), shading='flat')
        pos = np.where(np.isfinite(layer) & (layer > 0), layer, np.nan)
        im = ax.pcolormesh(edges, edges, pos.T, cmap=cmap, vmin=0, vmax=l_star, shading='flat')
        gi, Si = C.cspace_sdf_slices(A, B, [thc], D, res=0.01)
        ax.contour(gi, gi, Si[0].T, levels=[0], colors='k', linewidths=1.0)
        Z = f.eval(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, thc)])).reshape(X.shape)
        ax.contour(X, Y, Z, levels=[l_star], colors='r', linewidths=1.0)
        ax.set_aspect('equal')
        ax.set_xlim(-D, D)
        ax.set_ylim(-D, D)
        ax.set_xlabel(r'$t_x$ [m]')
        if n == 0:
            ax.set_ylabel(r'$t_y$ [m]')
        tag = ' (cell attaining $l^\\star$)' if n == 0 else ''
        ax.set_title(f"({'bc'[n]}) $\\theta$-cell at {np.degrees(thc):.0f}$^\\circ${tag}", fontsize=10)
    fig.colorbar(im, cax=fig.add_subplot(gs[3]), label=r'$\ell_c$ [m] (cells with $\ell_c>0$)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'percell.png'), dpi=160)
    print(f">>> saved tro/figs/percell.png (slices at theta cells {slices}, l* {l_star:.4f})")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('what', choices=['cache', 'repr', 'stats', 'dock', 'unicycle', 'static', 'sens', 'percell', 'dock2', 'dock2dt', 'dock2res', 'all'])
    ap.add_argument('--trials', type=int, default=20)
    a = ap.parse_args()
    todo = ['cache', 'repr', 'stats', 'dock', 'unicycle', 'static'] if a.what == 'all' else [a.what]
    for w in todo:
        t0 = time.perf_counter()
        print(f"\n######## {w}", flush=True)
        {'cache': build_cache, 'repr': exp_repr, 'stats': lambda: exp_stats(a.trials), 'dock': exp_dock,
         'unicycle': exp_unicycle, 'static': exp_static, 'sens': exp_sens, 'percell': exp_percell_fig, 'dock2': exp_dock2, 'dock2dt': exp_dock2_dt, 'dock2res': exp_dock2_res}[w]()
        print(f"######## {w} done in {time.perf_counter() - t0:.0f} s", flush=True)
