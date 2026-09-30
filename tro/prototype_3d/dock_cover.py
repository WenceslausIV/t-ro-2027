"""
Planar docking (Sec. VIII-C of the paper) with the SURFACE-COVER construction instead of the tabulated
configuration-space field: one workspace SDF per body, certified levels that enclose the ground-truth
polygons, a Bernstein cover of the craft's certified curve, and corner barriers against the station
field with a C^2 curvature majorant. Same simulator, nominal controller, start, and goal as the
tabulated runs (cspace_experiments.simulate_team, DOCK2_START / DOCK2_GOAL).

2D fields are tensor cubic B-splines of proto3d.Spline3 that are constant in a thin z layer and
evaluated at z = 0; the cover keeps the boxes whose z range contains 0 and uses their xy squares.

    python prototype_3d/dock_cover.py [box_side_mm ...]
"""
import json
import os
import sys
import time

import numpy as np
from matplotlib.path import Path
from scipy.spatial import cKDTree

import franka3d as F
from proto3d import Spline3, surface_boxes

sys.path.insert(0, os.path.dirname(F.HERE))
sys.path.insert(0, os.path.join(os.path.dirname(F.HERE), 'planar'))  # planar modules
import cspace_sdf_cbf_compare as C                                      # noqa: E402
from cspace_experiments import DOCK2_GOAL, DOCK2_START, simulate_team   # noqa: E402
from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH, craft_body, station_body   # noqa: E402

J2 = np.array([[0., -1.], [1., 0.]])
CORNERS2 = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]], float) / 2
ACT = .05


def refine(P, ds=.001):
    Q = np.vstack([P, P[:1]])
    out = []
    for a, b in zip(Q[:-1], Q[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / ds)))
        out.append(a + (b - a) * np.arange(n)[:, None] / n)
    return np.vstack(out)


def seg_distance(P, A, B):
    """Exact distance from each point of P to the polygon with edges A[k] -> B[k]."""
    d = np.full(len(P), np.inf)
    for s in range(0, len(A), 256):
        a, b = A[s:s + 256], B[s:s + 256]
        ab = b - a
        t = np.clip(np.einsum('pkd,kd->pk', P[:, None] - a, ab) / np.einsum('kd,kd->k', ab, ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(P[:, None] - a - t[..., None] * ab, axis=2).min(axis=1))
    return d


def body_field(G, h, margin, training_spacing=None):
    lo, hi = G.min(axis=0) - margin, G.max(axis=0) + margin
    tree, path = cKDTree(refine(G)), Path(G)

    def sdf(P):
        d = tree.query(P[:, :2])[0]
        return np.where(path.contains_points(P[:, :2]), -d, d)

    spacing = h / 2 if training_spacing is None else training_spacing
    f = Spline3(np.r_[lo, -1.5 * h], np.r_[hi, 1.5 * h], h).fit(sdf, spacing)
    return f


def level_2d(f, Mf, G, s=.002):
    """Certified l with the ground-truth polygon boundary inside {phi < l} (Proposition 1, n = 2)."""
    delta = s * np.sqrt(2) / 2
    lo, hi = G.min(axis=0) - s, G.max(axis=0) + s
    X, Y = np.meshgrid(np.arange(lo[0], hi[0] + 1e-9, s), np.arange(lo[1], hi[1] + 1e-9, s), indexing='ij')
    P = np.c_[X.ravel(), Y.ravel()]
    near = P[cKDTree(G).query(P)[0] <= delta + np.linalg.norm(np.diff(np.vstack([G, G[:1]]), axis=0), axis=1).max()]
    band = near[seg_distance(near, G, np.roll(G, -1, axis=0)) <= delta]
    P3 = np.c_[band, np.zeros(len(band))]
    v, g = f.eval(P3, order=1)
    return float(np.max(v + np.linalg.norm(g[:, :2], axis=1) * delta + .5 * Mf.eval(P3) * delta ** 2)) + 1e-6


def tight_level_2d(f, Mf, G, eps=1e-4):
    """Tightest certified level for a polygon boundary (branch-and-bound on its edges, cf. tight_level_mesh):
    UB = phi(mid) + |grad phi(mid)| L/2 + kappa L^2/8 per edge of length L; LB = phi at points on the edges."""
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('eps must be finite and positive')
    A, B = G, np.roll(G, -1, axis=0)
    if hasattr(f, 'h'):
        while True:
            split = np.linalg.norm(B - A, axis=1) / 2 > f.h
            if not split.any():
                break
            m = (A[split] + B[split]) / 2
            A, B = (np.vstack([A[~split], A[split], m]),
                    np.vstack([B[~split], m, B[split]]))
    z = lambda P: np.c_[P, np.zeros(len(P))]
    LB = float(f.eval(z(A)).max())
    while len(A):
        m = (A + B) / 2
        r = np.linalg.norm(B - A, axis=1) / 2
        v, g = f.eval(z(m), order=1)
        LB = max(LB, float(v.max()))
        UB = v + np.linalg.norm(g[:, :2], axis=1) * r + .5 * Mf.eval(z(m)) * r ** 2
        keep = UB >= LB + eps
        A, B, m = A[keep], B[keep], m[keep]
        A, B = np.vstack([A, m]), np.vstack([m, B])
    return LB + eps


def cover_2d(f, level, side=None):
    lo, sz = surface_boxes(f, side, level=level)
    keep = (lo[:, 2] <= 0) & (lo[:, 2] + sz > 0)
    return lo[keep, :2] + sz / 2, sz


def make_cover_fn(fS, MS, lS, PC, side, bounds=None):
    """Corner barriers of the craft's cover (robot j) against the station field (robot i), both SE(2).

    Certified cluster pruning (as for the Franka links, franka3d.barrier_rows): the boxes are grouped
    into clusters with center cc and radius R (covering the boxes' balls); a cluster whose ball lies in
    the domain is skipped if phi(cc) - G R - M r^2 / 2 - l >= ACT, with G, M certified bounds of |grad phi|
    and of the majorant on the cells the ball and the majorant's support can reach. This implies the
    per-box test below for all its boxes, so the barrier rows are unchanged. `bounds` = (Gn, Mn) may be
    passed to share them between pairs with the same field."""
    r = side * np.sqrt(2) / 2
    dlo, dhi = fS.lo[:2], (fS.lo + fS.K * fS.h)[:2]
    off = CORNERS2 * side
    groups, cen, crad = F.clusters(PC, r)
    Gn, Mn = bounds if bounds is not None else F.neighborhood_bounds(fS, float(crad.max()))

    def fn(xi, xj):
        Ri, Rj = C.rot(xi[2]), C.rot(xj[2])
        ccb = (cen @ Rj.T + xj[:2] - xi[:2]) @ Ri          # cluster centers in the station frame
        keep = ~np.all((ccb - crad[:, None] > dlo) & (ccb + crad[:, None] < dhi), axis=1)
        ins = np.flatnonzero(~keep)
        if len(ins):
            P3 = np.c_[ccb[ins], np.zeros(len(ins))]
            ci = fS.cell_of(P3)
            low = (fS.eval(P3) - Gn[ci[:, 0], ci[:, 1], ci[:, 2]] * crad[ins]
                   - .5 * Mn[ci[:, 0], ci[:, 1], ci[:, 2]] * r ** 2 - lS)
            keep[ins[low < ACT]] = True
        if not keep.any():
            return np.zeros((0, 6)), np.zeros(0)
        sub = np.concatenate([groups[k] for k in np.flatnonzero(keep)])
        cw = PC[sub] @ Rj.T + xj[:2]                        # box centers (world)
        cb = (cw - xi[:2]) @ Ri                             # ... in the station frame
        ins = np.flatnonzero(np.all((cb - r > dlo) & (cb + r < dhi), axis=1))
        if not len(ins):
            return np.zeros((0, 6)), np.zeros(0)
        P3 = np.c_[cb[ins], np.zeros(len(ins))]
        val, g = fS.eval(P3, order=1)
        near = ins[val - np.linalg.norm(g[:, :2], axis=1) * r - .5 * MS.eval(P3) * r ** 2 - lS < ACT]
        if not len(near):
            return np.zeros((0, 6)), np.zeros(0)
        P3 = np.c_[cb[near], np.zeros(len(near))]
        v2, g2, H2 = fS.eval(P3, order=2)
        Mv, gM = MS.eval(P3, order=1)
        g2, H2, gM = g2[:, :2], H2[:, :2, :2], gM[:, :2]
        d = off @ (Ri.T @ Rj).T                             # corner offsets in the station frame (4, 2)
        hh = v2[:, None] + g2 @ d.T - .5 * Mv[:, None] * r ** 2 - lS
        bi, vi = np.nonzero(hh < ACT)
        if not len(bi):
            return np.zeros((0, 6)), np.zeros(0)
        dd = d[vi]
        a = g2[bi] + np.einsum('nij,nj->ni', H2[bi], dd) - .5 * r ** 2 * gM[bi]
        aw, gw, dw = a @ Ri.T, g2[bi] @ Ri.T, dd @ Ri.T
        c = cw[near][bi]
        R6 = np.zeros((len(bi), 6))
        R6[:, 3:5] = aw
        R6[:, 5] = np.einsum('nd,nd->n', aw, (c - xj[:2]) @ J2.T) + np.einsum('nd,nd->n', gw, dw @ J2.T)
        R6[:, 0:2] = -aw
        R6[:, 2] = -np.einsum('nd,nd->n', aw, (c - xi[:2]) @ J2.T) - np.einsum('nd,nd->n', gw, dw @ J2.T)
        return R6, hh[bi, vi]

    return fn


def make_summed_fn(fS, MS, lS, bodyC, bounds=None, eta=ACT, umax=(1., 1., 2., 1., 1., 2.)):
    """Summed-field barriers (summed.py) of the cover of robot j (prepared by summed.prepare_body, dim 2)
    against the field of robot i, both SE(2) single integrators with inputs (v_i, w_i, v_j, w_j).
    Returns fn(xi, xj) -> (A (n,6), T (n,3), C (n,), boxes, lower bound of h~, E (3,6)) or None; the rows read
    A u + T a + C >= 0 with auxiliary variables a >= |E u| (relative twist of the pair)."""
    import summed as SM
    PC, side, r = bodyC['PC'], bodyC['side'], bodyC['r']
    dlo, dhi = fS.lo[:2], (fS.lo + fS.K * fS.h)[:2]
    groups, cen, crad = F.clusters(PC, r)
    gmin = np.array([bodyC['lowA'][g].min() for g in groups])
    Gn, _ = bounds if bounds is not None else F.neighborhood_bounds(fS, float(crad.max()))
    if not hasattr(fS, '_M3') or getattr(fS, '_M3_reach', -1.) != r:
        fS._M3 = SM.field_bound(fS, r)
        fS._M3_reach = r
    M3S = fS._M3
    # A native patch may be wider than a cell of the other body's field.
    if r > fS.h:
        MS = SM.LinearMajorant(fS, fS.hessian_bound(cellwise=True), r)
    umax = np.asarray(umax, float)
    Wg = np.zeros((6, 2, 2))
    Wg[2], Wg[5] = -J2, J2

    def fn(xi, xj):
        Ri, Rj = C.rot(xi[2]), C.rot(xj[2])
        ccb = (cen @ Rj.T + xj[:2] - xi[:2]) @ Ri
        inside = np.all((ccb - crad[:, None] > dlo) & (ccb + crad[:, None] < dhi), axis=1)
        keep = ~inside
        ins = np.flatnonzero(inside)
        if len(ins):
            P3 = np.c_[ccb[ins], np.zeros(len(ins))]
            ci = fS.cell_of(P3)
            low = fS.eval(P3) - Gn[ci[:, 0], ci[:, 1], ci[:, 2]] * crad[ins] - lS + gmin[ins]
            keep[ins[low < eta]] = True
        if not keep.any():
            return None
        sub = np.concatenate([groups[k] for k in np.flatnonzero(keep)])
        cw = PC[sub] @ Rj.T + xj[:2]
        cb = (cw - xi[:2]) @ Ri
        ok = np.all((cb - r > dlo) & (cb + r < dhi), axis=1)
        sub, cw, cb = sub[ok], cw[ok], cb[ok]
        if not len(sub):
            return None
        P3 = np.c_[cb, np.zeros(len(cb))]
        val, g = fS.eval(P3, order=1)
        pre = val - np.linalg.norm(g[:, :2], axis=1) * r - .5 * MS.eval(P3) * r ** 2 - lS + bodyC['lowA'][sub] < eta
        sub, cw, cb = sub[pre], cw[pre], cb[pre]
        if not len(sub):
            return None
        n = len(sub)
        Jc = np.zeros((n, 6, 2))
        # center velocity in i's frame: R_i^T (v_j + w_j J (c - p_j) - v_i - w_i J (c - p_i))
        Jc[:, 0], Jc[:, 1] = -Ri.T[:, 0], -Ri.T[:, 1]
        Jc[:, 3], Jc[:, 4] = Ri.T[:, 0], Ri.T[:, 1]
        Jc[:, 2] = -((cw - xi[:2]) @ J2.T) @ Ri
        Jc[:, 5] = ((cw - xj[:2]) @ J2.T) @ Ri
        # relative twist of j with respect to i at p_j (world): V = v_j - v_i - w_i J (p_j - p_i), Omega = w_j - w_i;
        # every box point moves (relative to i) with speed <= |V_x| + |V_y| + |Omega| (|c - p_j| + r)
        Jp = J2 @ (xj[:2] - xi[:2])
        E = np.array([[-1., 0., -Jp[0], 1., 0., 0.], [0., -1., -Jp[1], 0., 1., 0.], [0., 0., -1., 0., 0., 1.]])
        Sb = np.c_[np.ones((n, 2)), np.linalg.norm(cw - xj[:2], axis=1) + r]
        res = SM.rows(bodyC, sub, Ri.T @ Rj, cb, fS, lS, M3S, Jc, Wg, Sb, eta, umax, np.abs(E) @ umax)
        if res is None:
            return None
        A6, Tk, C6, nbox, hl, _ = res
        return A6, Tk, C6, nbox, hl, E

    return fn


def main():
    corner = '--corner' in sys.argv                 # previous corner barriers; default: summed-field barriers
    tag = 'cover' if corner else 'summed'
    sides = [float(a) / 1e3 for a in sys.argv[1:] if not a.startswith('--')] or [.01]
    if not corner:
        tag += '_native' if sides == [.01] else '_subdivided'
    GS, GC = station_body(), craft_body()
    t0 = time.perf_counter()
    fS, fC = body_field(GS, .0125, .15), body_field(GC, .01, .06)
    MS, MC = F.hessian_majorant(fS), F.hessian_majorant(fC)
    lS, lC = tight_level_2d(fS, MS, GS), tight_level_2d(fC, MC, GC)
    t_setup = time.perf_counter() - t0
    out = dict(l_station_mm=1e3 * lS, l_craft_mm=1e3 * lC, t_setup=t_setup)
    print(out, flush=True)
    for side in sides:
        t0 = time.perf_counter()
        PC, sz = cover_2d(fC, lC, side)
        t_cover = time.perf_counter() - t0
        if corner:
            fn = make_cover_fn(fS, MS, lS, PC, sz)
        else:
            import summed as SM
            fn = make_summed_fn(fS, MS, lS, SM.prepare_body(fC, lC, PC, sz, 2))
        shapes = [type('S', (), dict(rho=float(np.linalg.norm(G, axis=1).max())))() for G in (GS, GC)]
        L = simulate_team([GS, GC], shapes, {(0, 1): fn}, DOCK2_START, DOCK2_GOAL, 'cover' if corner else 'summed',
                          steps=900, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        xB = L['final'][1]
        tip = xB[:2] + NOSE_TIP_X * np.array([np.cos(xB[2]), np.sin(xB[2])])
        u = np.array(L['u_t'])
        res = dict(n_boxes=len(PC), side_mm=1e3 * sz, t_cover=t_cover, tip_depth=MOUTH_X - tip[0],
                   seat_shortfall=SEAT_DEPTH - (MOUTH_X - tip[0]),
                   final_gt_mm=1e3 * C.true_distance(GS, GC @ C.rot(xB[2]).T + xB[:2]),
                   min_gt_mm=1e3 * L['min_gt'], min_h_mm=1e3 * float(np.nanmin(L['h_t'])), tv=L['tv'],
                   t_med_ms=1e3 * float(np.median(L['t_ctrl'])), t_max_ms=1e3 * float(np.max(L['t_ctrl'])),
                   t_p95_ms=1e3 * float(np.percentile(L['t_ctrl'], 95)), slack=L['slack_steps'],
                   rows_max=int(max(L.get('n_rows', [0]))), boxes_max=int(max(L.get('n_box', [0]))))
        out[f'side_{1e3 * side:g}'] = res
        print(res, flush=True)
        np.savez_compressed(os.path.join(F.HERE, f'dock_{tag}_traj_{1e3 * side:g}.npz'),
                            traj=np.array(L['traj']), h=np.array(L['h_t']), d=np.array(L['d_t']),
                            nb=np.array(L.get('n_box', [])), nr=np.array(L.get('n_rows', [])))
    json.dump(out, open(os.path.join(F.HERE, f'dock_{tag}_results.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
