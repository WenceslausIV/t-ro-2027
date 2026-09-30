"""
Baselines for the Franka experiment on the trials of franka_results.json:
  spheres : the collision spheres of the RDF repository (RDF/panda_layer/franka_sphere.yaml),
            h = d_O(center) - radius - margin, one barrier per sphere and obstacle;
  points  : N surface samples per link mesh (as the 256 points of [li2024representing]),
            h = d_O(point) - margin, one barrier per sample and obstacle.
Both use the EXACT analytic obstacle SDF (gradient by central differences) and a 1-cm margin; safety is
evaluated as in franka3d.py (certified lower bound of the mesh-obstacle distance).

    python prototype_3d/franka_baselines.py [spheres|points ...]
"""
import json
import os
import sys
import time

import numpy as np
import trimesh
import yaml

import franka3d as F
from proto3d import sdf_boxes
from sdf_cbf_utils import solve_ldp_qp

MARGIN, N_POINTS = .01, 256
SPHERE_YAML = os.path.join(F.HERE, 'RDF', 'panda_layer', 'franka_sphere.yaml')


def sdf_grad(P, shape, eps=1e-5):
    """Exact obstacle SDF and its central-difference gradient; shape: rounded boxes or an SDF function."""
    fn = shape if callable(shape) else (lambda X: sdf_boxes(X, shape))
    d = fn(P)
    g = np.stack([(fn(P + eps * e) - fn(P - eps * e)) / (2 * eps) for e in np.eye(3)], 1)
    return d, g


def primitives(links, kind):
    """Per link: points in the link frame and radii."""
    out = []
    if kind == 'spheres':
        conf = yaml.safe_load(open(SPHERE_YAML))['collision_spheres']
        names = {i: f'panda_link{i}' for i in range(1, 8)}
        names[8] = 'panda_hand'
        for L in links:
            sp = conf[names[L['frame']]]
            out.append((np.array([s['center'] for s in sp], float), np.array([s['radius'] for s in sp], float)))
    else:
        for L in links:
            P = trimesh.sample.sample_surface_even(trimesh.Trimesh(L['V'], L['F']), N_POINTS, seed=3)[0]
            out.append((P, np.zeros(len(P))))
    return out


def rows(q, links, obst, prim):
    T, Z, O = F.fk(q)
    ZxO = np.cross(Z, O)
    Gs, hs = [], []
    for L, (P, rad) in zip(links, prim):
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        W = P @ R.T + p
        for o in obst.values():
            m = np.all((W > o['blo'] - F.ACT - rad[:, None] - MARGIN) & (W < o['bhi'] + F.ACT + rad[:, None] + MARGIN),
                       axis=1)
            if not m.any():
                continue
            d, g = sdf_grad(W[m], o['sdf'] if 'sdf' in o else o['boxes'])
            h = d - rad[m] - MARGIN
            keep = h < F.ACT
            if not keep.any():
                continue
            x, g, h = W[m][keep], g[keep], h[keep]
            row = np.cross(x, g) @ Z.T - g @ ZxO.T
            row[:, L['n_joints']:] = 0.
            Gs.append(row); hs.append(h)
    if not Gs:
        return np.zeros((0, F.DOF)), np.zeros(0)
    return np.vstack(Gs), np.concatenate(hs)


def kkt_rows(q, links, obst, state):
    """ICRA'26-style closest point per (link, obstacle) pair: y* = argmin phi_O(R y + p) s.t. phi_A(y) = l_A
    (SLSQP, warm-started from the previous solution, initialized at the best box center), h = phi_O(x*) - l_O,
    row = grad phi_O(x*)^T J(x*) (envelope theorem)."""
    from scipy.optimize import minimize
    T, Z, O = F.fk(q)
    ZxO = np.cross(Z, O)
    Gs, hs = [], []
    for L in links:
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        fA = L.setdefault('fs', F.spline_from(L['f']))
        for name, o in obst.items():
            if np.any(p + L['rad'] < o['blo'] - F.ACT - .05) or np.any(p - L['rad'] > o['bhi'] + F.ACT + .05):
                state.pop((L['frame'], name), None)
                continue
            key = (L['frame'], name)
            if key not in state:
                W = L['PC'] @ R.T + p
                ins = np.all((W > o['dlo']) & (W < o['dhi']), axis=1)
                if not ins.any():
                    continue
                cand = np.flatnonzero(ins)
                state[key] = L['PC'][cand[np.argmin(o['f'].eval(W[cand]))]]
            fO, lA = o['f'], L['l']
            obj = lambda y: fO.eval((R @ y + p)[None])[0]
            jac = lambda y: fO.eval((R @ y + p)[None], order=1)[1][0] @ R
            con = dict(type='eq', fun=lambda y: fA.eval(y[None])[0] - lA,
                       jac=lambda y: fA.eval(y[None], order=1)[1][0])
            res = minimize(obj, state[key], jac=jac, constraints=[con], method='SLSQP',
                           options=dict(maxiter=30, ftol=1e-10))
            state[key] = res.x
            x = R @ res.x + p
            if np.any(x < o['dlo']) or np.any(x > o['dhi']):
                continue
            val, g = fO.eval(x[None], order=1)
            h = val[0] - o['l']
            if h >= F.ACT:
                continue
            row = np.cross(x, g[0]) @ Z.T - g[0] @ ZxO.T
            row[L['n_joints']:] = 0.
            Gs.append(row[None]); hs.append([h])
    if not Gs:
        return np.zeros((0, F.DOF)), np.zeros(0)
    return np.vstack(Gs), np.concatenate(hs)


def simulate(q0, qg, links, obst, prim, steps=1000):
    q = np.array(q0, float)
    log = dict(t=[], rows=[], h=[], gap=[], witness=[], slack=0, reached=None)
    u = None
    state = {}
    for k in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        G, h = kkt_rows(q, links, obst, state) if prim == 'kkt' else rows(q, links, obst, prim)
        u, s = solve_ldp_qp(u_nom, np.vstack([G, F.BOX_G]), np.r_[-F.GAMMA * h, F.BOX_H], len(h), u)
        log['t'].append(time.perf_counter() - t0)
        log['slack'] += s > 0
        log['rows'].append(len(h)); log['h'].append(float(h.min()) if len(h) else np.inf)
        lower, witness = F.real_gap(q, links, obst, return_witness=True)
        log['gap'].append(lower); log['witness'].append(witness)
        q = q + F.DT * u
        if np.linalg.norm(qg - q) < .05:
            log['reached'] = (k + 1) * F.DT
            break
    return log


def main():
    global MARGIN
    args = sys.argv[1:]
    if '--margin' in args:                              # e.g. --margin 0 -> results under 'spheres_m0'
        i = args.index('--margin')
        MARGIN = float(args[i + 1]) / 1e3
        del args[i:i + 2]
    kinds = args or ['spheres', 'points']
    links, obst, _ = F.build()
    trials = json.load(open(os.path.join(F.HERE, 'franka_results.json')))['trials']
    path = os.path.join(F.HERE, 'franka_baselines.json')
    out = json.load(open(path)) if os.path.exists(path) else {}
    for kind in kinds:
        prim = 'kkt' if kind == 'kkt' else primitives(links, kind)
        res = []
        for tr in trials:
            L = simulate(np.array(tr['q0']), np.array(tr['qg']), links, obst, prim)
            r = dict(trial=tr['trial'], min_gap_mm=1e3 * min(L['gap']), reached=L['reached'],
                     min_centroid_sdf_mm=1e3 * min(L['witness']),
                     slack_steps=int(L['slack']), min_h_mm=1e3 * min(L['h']),
                     t_med_ms=1e3 * float(np.median(L['t'])), t_max_ms=1e3 * float(np.max(L['t'])),
                     rows_max=int(max(L['rows'])))
            res.append(r)
            print(kind, MARGIN, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        gaps = np.array([r['min_gap_mm'] for r in res])
        confirmed = np.array([r['min_centroid_sdf_mm'] < 0 for r in res])
        summ = dict(n=len(res), collisions=int(confirmed.sum()),
                    negative_bound_trials=int((gaps < 0).sum()),
                    unresolved_trials=int(((gaps <= 0) & ~confirmed).sum()),
                    collision_definition='negative SDF at an actual mesh centroid; negative bounds alone are unresolved',
                    min_gap_mm=float(gaps.min()),
                    median_min_gap_mm=float(np.median(gaps)),
                    reached=int(sum(r['reached'] is not None for r in res)),
                    slack_trials=int(sum(r['slack_steps'] > 0 for r in res)),
                    t_med_ms=float(np.median([r['t_med_ms'] for r in res])),
                    t_max_ms=float(max(r['t_max_ms'] for r in res)))
        print('summary', kind, MARGIN, json.dumps(summ), flush=True)
        key = kind if abs(MARGIN - .01) < 1e-12 else f'{kind}_m{1e3 * MARGIN:g}'
        out = json.load(open(path)) if os.path.exists(path) else {}
        out[key] = dict(summary=summ, trials=res, margin=MARGIN,
                         n_primitives=(len(links) * len(obst) if kind == 'kkt'
                                       else int(sum(len(p[0]) for p in prim))))
        json.dump(out, open(path, 'w'), indent=1)


if __name__ == '__main__':
    main()
