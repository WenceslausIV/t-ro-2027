"""
Dual-arm collision avoidance: two Franka Pandas facing each other share a workspace. For every pair
(link a of arm 1, link b of arm 2), the certified surface of link a is covered by boxes (franka3d cache)
and bounded against the spline SDF of link b (Sec. VI of the paper); both links move, so each corner
barrier is affine in the joint velocities of both arms:

  h = phi_b(c) + g^T d - M_b(c) r^2/2 - l_b      (c, d = v - c in link b's frame)
  dh/dt = sum_j z1_j.(c x a + d x g) - a.(z1_j x o1_j)  qd1_j   -  [same with arm 2]  qd2_j
  with a = R_b (g + H d - r^2/2 grad M_b), g = R_b grad phi_b (world frame).

Evaluation: approximate mesh distance every step (KD trees of nominally 3-mm mesh samples, queried within 3 cm), then the certified
lower bound (exact point-to-triangle distances from the subdivided triangles of arm 1 to the meshes of
arm 2, minus the triangle radii) at up to 20 steps whose approximate distance is within 5 mm of the
trial's minimum. Sample-based pruning uses an independently bounded surface-cover radius.

    python prototype_3d/dual3d.py [n_trials]
"""
import json
import os
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

import franka3d as F
from sdf_cbf_utils import solve_ldp_qp

BASE2 = np.array([[-1., 0., 0., .9], [0., -1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
DOF2 = 2 * F.DOF


def fk2(q):
    """Frames of both arms and their joint axes/origins (world)."""
    T1, Z1, O1 = F.fk(q[:F.DOF])
    T2, Z2, O2 = F.fk(q[F.DOF:])
    T2 = [BASE2 @ T for T in T2]
    Z2 = Z2 @ BASE2[:3, :3].T
    O2 = O2 @ BASE2[:3, :3].T + BASE2[:3, 3]
    return (T1, Z1, O1), (T2, Z2, O2)


def build(cover_side=None):
    links, _, info = F.build(cover_side=cover_side)
    for L in links:
        f = F.spline_from(L['f'])
        L['fs'], L['Ms'] = f, F.hessian_majorant(f)
        L['dlo'], L['dhi'] = f.lo, f.lo + f.K * f.h
        mesh = trimesh.Trimesh(L['V'], L['F'])
        L['smp'] = np.vstack([trimesh.sample.sample_surface_even(mesh, int(mesh.area / .003 ** 2), seed=4)[0],
                              L['V']])
        L['tree'] = cKDTree(L['smp'])
        L['bb'] = (L['smp'].min(axis=0) - .03, L['smp'].max(axis=0) + .03)
        L['Vs'], L['Fs'] = trimesh.remesh.subdivide_to_size(L['V'], L['F'], max_edge=.004)
        triangles = L['Vs'][L['Fs']]
        centers = triangles.mean(axis=1)
        radii = np.linalg.norm(triangles - centers[:, None], axis=2).max(axis=1)
        # Every mesh point lies within this distance of an actual KD-tree sample.
        L['sample_cover_radius'] = float((L['tree'].query(centers)[0] + radii).max())
    reach = max(L['crad'].max() for L in links)
    for L in links:
        L['Gn'], L['Mn'] = F.neighborhood_bounds(L['fs'], reach)
    return links, info


def barrier_rows(q, links):
    (T1, Z1, O1), (T2, Z2, O2) = fk2(q)
    Z1xO1, Z2xO2 = np.cross(Z1, O1), np.cross(Z2, O2)
    rows, hs = [], []
    for A in links:
        RA, pA = T1[A['frame']][:3, :3], T1[A['frame']][:3, 3]
        for B in links:
            RB, pB = T2[B['frame']][:3, :3], T2[B['frame']][:3, 3]
            if np.linalg.norm(pA - pB) > A['rad'] + B['rad'] + .06 + F.ACT:
                continue
            # arm-1 boxes in link b's frame; clusters first (certified pruning)
            RbA = RB.T @ RA
            tb = RB.T @ (pA - pB)
            cc = A['cc'] @ RbA.T + tb
            keep = ~np.all((cc - A['crad'][:, None] > B['dlo']) & (cc + A['crad'][:, None] < B['dhi']), axis=1)
            ins = np.flatnonzero(~keep)
            if len(ins):
                ci = B['fs'].cell_of(cc[ins])
                low = (B['fs'].eval(cc[ins]) - B['Gn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['crad'][ins]
                       - .5 * B['Mn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['r'] ** 2 - B['l'])
                keep[ins[low < F.ACT]] = True
            if not keep.any():
                continue
            sub = np.concatenate([A['groups'][j] for j in np.flatnonzero(keep)])
            C = A['PC'][sub] @ RbA.T + tb
            idx = np.flatnonzero(np.all((C - A['r'] > B['dlo']) & (C + A['r'] < B['dhi']), axis=1))
            if not len(idx):
                continue
            val, g = B['fs'].eval(C[idx], order=1)
            near = idx[val - np.linalg.norm(g, axis=1) * A['r'] - .5 * B['Ms'].eval(C[idx]) * A['r'] ** 2
                       - B['l'] < F.ACT]
            if not len(near):
                continue
            c = C[near]
            v2, g2, H2 = B['fs'].eval(c, order=2)
            Mv, gM = B['Ms'].eval(c, order=1)
            d = (F.CORNERS * A['side']) @ RbA.T                        # corner offsets in link b's frame
            hh = v2[:, None] + g2 @ d.T - .5 * Mv[:, None] * A['r'] ** 2 - B['l']
            bi, vi = np.nonzero(hh < F.ACT)
            if not len(bi):
                continue
            dd = d[vi]
            a = g2[bi] + np.einsum('nij,nj->ni', H2[bi], dd) - .5 * A['r'] ** 2 * gM[bi]
            aw, gw = a @ RB.T, g2[bi] @ RB.T                             # world frame
            cw, dw = c[bi] @ RB.T + pB, dd @ RB.T
            w = np.cross(cw, aw) + np.cross(dw, gw)
            r1 = w @ Z1.T - aw @ Z1xO1.T
            r1[:, A['n_joints']:] = 0.
            r2 = -(w @ Z2.T - aw @ Z2xO2.T)
            r2[:, B['n_joints']:] = 0.
            rows.append(np.hstack([r1, r2])); hs.append(hh[bi, vi])
    if not rows:
        return np.zeros((0, DOF2)), np.zeros(0)
    return np.vstack(rows), np.concatenate(hs)


def pair_frames(q, links):
    (T1, _, _), (T2, _, _) = fk2(q)
    for A in links:
        for B in links:
            pA, pB = T1[A['frame']][:3, 3], T2[B['frame']][:3, 3]
            if np.linalg.norm(pA - pB) <= A['rad'] + B['rad'] + .02:
                yield A, B, T2[B['frame']][:3, :3].T @ T1[A['frame']][:3, :3], \
                    T2[B['frame']][:3, :3].T @ (pA - pB)


def approx_gap(q, links):
    g = np.inf
    for A, B, Rab, tab in pair_frames(q, links):
        P = A['smp'] @ Rab.T + tab
        P = P[np.all((P > B['bb'][0]) & (P < B['bb'][1]), axis=1)]
        if len(P):
            g = min(g, float(B['tree'].query(P, distance_upper_bound=.03)[0].min()))
    return g


def certified_gap(q, links):
    """Whole-mesh surface distance lower bound; nonpositive does not prove collision.

    Omitted link pairs have separation >2 cm because their enclosing spheres do.
    Omitted centroids retain the AABB/sample-cover lower bound minus triangle radius.
    """
    g = .02
    for A, B, Rab, tab in pair_frames(q, links):
        P = A['gt'] @ Rab.T + tab
        m = np.all((P > B['bb'][0]) & (P < B['bb'][1]), axis=1)
        if (~m).any():
            g = min(g, .03 - float(A['gt_r'][~m].max()))
        idx = np.flatnonzero(m)
        if len(idx):
            sampled_far = B['tree'].query(P[idx], distance_upper_bound=.03)[0] >= .03
            far = idx[sampled_far]
            if len(far):
                floor = .03 - B['sample_cover_radius'] - float(A['gt_r'][far].max())
                if floor > 0:
                    g = min(g, floor)
                    m[far] = False
        if m.any():
            d = F.mesh_distance(B['Vs'], B['Fs'], P[m], .03)
            g = min(g, float((d - A['gt_r'][m]).min()))
    return g


BOX_G = np.vstack([np.eye(DOF2), -np.eye(DOF2)])
BOX_H = -F.QD_MAX * np.ones(2 * DOF2)


def simulate(q0, qg, links, method, steps=1000, record=False):
    q = np.array(q0, float)
    log = dict(t=[], rows=[], h=[], gap=[], slack=0, reached=None, q=[], boxes=[])
    u = z = None
    for k in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        if method == 'ours':                  # summed-field barriers with Bernstein coefficient rows (summed3d.py)
            import summed as SM
            import summed3d as S3
            Ar, Tr, Cr, nb, hl, Er = S3.dual_rows(q, links)
            u, s, z = SM.solve(u_nom, Ar, Tr, Cr, Er, BOX_G, BOX_H, z)
            log['slack'] += s > 0
            log['rows'].append(len(Cr)); log['h'].append(hl); log['boxes'].append(nb)
        elif method == 'corner':              # previous corner barriers
            G, h = barrier_rows(q, links)
            u, s = solve_ldp_qp(u_nom, np.vstack([G, BOX_G]), np.r_[-F.GAMMA * h, BOX_H], len(h), u)
            log['slack'] += s > 0
            log['rows'].append(len(h)); log['h'].append(float(h.min()) if len(h) else np.inf)
        else:
            u = u_nom
        log['t'].append(time.perf_counter() - t0)
        log['gap'].append(approx_gap(q, links))
        log['q'].append(q.copy())
        q = q + F.DT * u
        if np.linalg.norm(qg - q) < .05:
            log['reached'] = (k + 1) * F.DT
            break
    log['q'].append(q.copy())
    return log


def certify_log(log, links, window=.005, n_max=20):
    """Certified distance at the (at most n_max) steps with the smallest sampled distance within `window` of
    the minimum. The nominal sample spacing alone is not an error certificate."""
    g = np.array(log['gap'])
    cand = np.flatnonzero(g <= g.min() + window)
    cand = cand[np.argsort(g[cand])[:n_max]]
    return min(certified_gap(log['q'][k], links) for k in cand)


def sample_config(rng, links, clear=.05):
    mid, half = (F.Q_MIN + F.Q_MAX) / 2, (F.Q_MAX - F.Q_MIN) / 2 * .8
    while True:
        q = np.r_[rng.uniform(mid - half, mid + half), rng.uniform(mid - half, mid + half)]
        (T1, _, _), (T2, _, _) = fk2(q)
        if T1[8][2, 3] < .15 or T2[8][2, 3] < .15:
            continue
        if approx_gap(q, links) > clear and barrier_rows(q, links)[1].min(initial=1.) > 0:
            return q


def main():
    n_trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    links, info = build()
    import summed3d as S3
    S3.prep_all(links)
    rng = np.random.default_rng(0)
    res = []
    tried = 0
    reuse = os.path.join(F.HERE, 'dual_trials.json')   # fixed trial set: [q0, qg] pairs
    fixed = [(np.array(a), np.array(b)) for a, b in json.load(open(reuse))] if os.path.exists(reuse) else []
    while len(res) < n_trials:
        if fixed:
            q0, qg = fixed.pop(0)
        else:
            q0, qg = sample_config(rng, links), sample_config(rng, links)
        tried += 1
        nom = simulate(q0, qg, links, 'nominal')
        if min(nom['gap']) > .002:                     # select sampled near-contact, not a collision witness
            continue
        L = simulate(q0, qg, links, 'ours', record=True)
        r = dict(trial=len(res), q0=q0.tolist(), qg=qg.tolist(), nominal_min_gap_mm=1e3 * min(nom['gap']),
                 approx_min_gap_mm=1e3 * min(L['gap']), certified_min_gap_mm=1e3 * certify_log(L, links),
                 reached=L['reached'], slack_steps=int(L['slack']), min_h_mm=1e3 * min(L['h']),
                 t_med_ms=1e3 * float(np.median(L['t'])), t_max_ms=1e3 * float(np.max(L['t'])),
                 rows_max=int(max(L['rows'])), boxes_max=int(max(L['boxes'])))
        res.append(r)
        print(r['trial'], {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()
                           if k not in ('q0', 'qg')}, flush=True)
    gaps = np.array([r['certified_min_gap_mm'] for r in res])
    summ = dict(n=len(res), sampled_pairs=tried, unresolved_selected_state_trials=int((gaps <= 0).sum()),
                distance_scope='At most 20 selected states per trial; negative lower bounds are not collision witnesses',
                min_gap_mm=float(gaps.min()),
                median_min_gap_mm=float(np.median(gaps)), reached=int(sum(r['reached'] is not None for r in res)),
                slack_trials=int(sum(r['slack_steps'] > 0 for r in res)),
                t_med_ms=float(np.median([r['t_med_ms'] for r in res])),
                t_max_ms=float(max(r['t_max_ms'] for r in res)))
    print('summary', json.dumps(summ), flush=True)
    json.dump(dict(summary=summ, trials=res, base2=BASE2.tolist()),
              open(os.path.join(F.HERE, 'dual_results_native.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
