"""Continuous-boundary baselines on the 30 Franka trials (see results/continuous_baselines/PLAN.md).

    python prototype_3d/continuous_baselines.py --methods points_delta capsule spheres_enclosing [--edge-mm 10]
Every baseline certifies a continuous set that encloses each link mesh, against the EXACT obstacle SDF d_O
(union of boxes, 1-Lipschitz; no fitting error, which favors the baselines). With a primitive point x and a
geometric margin m such that every mesh point lies within m of some x:  d_O >= 0 on the mesh whenever
h = d_O(x) - m >= 0 for all primitives (1-Lipschitz). CBF rows grad d_O(x)^T J(x) u + gamma h >= 0 at the 10-ms
samples, rows only for h < 3 cm (as ours).
  points_delta       mesh remeshed to max edge e; its vertices; m = delta = e / sqrt(3) (any point of a triangle
                     with edges <= e is within e/sqrt(3) of a vertex). 3D Poisson-safety-function style.
  capsule            one capsule per link: segment on the principal axis of the mesh vertices, radius r = max vertex
                     distance (the capsule is convex, so it contains the mesh); segment sampled at spacing e,
                     m = r + e/2. Capsule-link style (e.g. forestry crane, EDF collision checking).
  spheres_enclosing  the RDF collision-sphere centers; every remeshed vertex goes to its nearest center,
                     radius = max assigned distance + delta, so the spheres contain the mesh; m = radius.
  spheres_kmeans     same enclosure with K k-means centers of the remeshed vertices per link (--spheres K).
Same nominal, dt, gamma, input bounds, horizon, and QP (summed.solve, DAQP, exact check) as
certificate_upgrades.py. Audit: franka3d.real_gap at every state (outside the timer).
Writes results/continuous_baselines/<method>/franka_XX.json and .npz.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import time
from pathlib import Path

import numpy as np
import trimesh
import yaml

import franka3d as F
import summed as S
from proto3d import sdf_boxes

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'continuous_baselines'
SPHERE_YAML = Path(F.HERE) / 'RDF' / 'panda_layer' / 'franka_sphere.yaml'


def remeshed_vertices(L, edge):
    V, Fc = trimesh.remesh.subdivide_to_size(np.asarray(L['V']), np.asarray(L['F']), max_edge=edge)
    e = np.linalg.norm(V[Fc] - V[np.roll(Fc, 1, axis=1)], axis=2).max()
    assert e <= edge * (1 + 1e-9)
    return V


def primitives(links, method, edge, n_spheres=16):
    """Per link: points (link frame) and margins m (every mesh point is within m of one of the points)."""
    delta = edge / np.sqrt(3)
    out = []
    if method == 'spheres_kmeans':
        from scipy.cluster.vq import kmeans2
    if method == 'spheres_enclosing':
        conf = yaml.safe_load(open(SPHERE_YAML))['collision_spheres']
        names = {i: f'panda_link{i}' for i in range(1, 8)}
        names[8] = 'panda_hand'
    for L in links:
        if method == 'points_delta':
            P = remeshed_vertices(L, edge)
            out.append((P, np.full(len(P), delta)))
        elif method == 'capsule':
            V = np.asarray(L['V'], float)
            c = V.mean(axis=0)
            axis = np.linalg.svd(V - c, full_matrices=False)[2][0]
            t = (V - c) @ axis
            a, b = c + t.min() * axis, c + t.max() * axis
            s = np.clip(((V - a) @ (b - a)) / ((b - a) @ (b - a)), 0, 1)
            r = np.linalg.norm(V - (a + s[:, None] * (b - a)), axis=1).max()
            n = max(2, int(np.ceil(np.linalg.norm(b - a) / edge)) + 1)
            P = a + np.linspace(0, 1, n)[:, None] * (b - a)
            spacing = np.linalg.norm(b - a) / (n - 1)
            out.append((P, np.full(n, r + spacing / 2)))
        else:
            Vr = remeshed_vertices(L, edge)
            if method == 'spheres_kmeans':
                C = kmeans2(Vr, n_spheres, seed=0, minit='++')[0]
            else:
                C = np.array([s_['center'] for s_ in conf[names[L['frame']]]], float)
            near = np.argmin(np.linalg.norm(Vr[:, None] - C[None], axis=2), axis=1)
            keep = np.unique(near)
            rad = np.array([np.linalg.norm(Vr[near == j] - C[j], axis=1).max() for j in keep]) + delta
            out.append((C[keep], rad))
    return out


def sdf_grad(P, boxes, eps=1e-6):
    d = sdf_boxes(P, boxes)
    g = np.stack([(sdf_boxes(P + eps * e, boxes) - sdf_boxes(P - eps * e, boxes)) / (2 * eps) for e in np.eye(3)], 1)
    return d, g


def rows(q, links, obst, prim):
    T, Z, O = F.fk(q)
    ZxO = np.cross(Z, O)
    A, C = [], []
    for L, (P, m) in zip(links, prim):
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        W = P @ R.T + p
        for o in obst.values():
            near = np.all((W > o['blo'] - F.ACT - m[:, None]) & (W < o['bhi'] + F.ACT + m[:, None]), axis=1)
            if not near.any():
                continue                                  # farther than ACT + m from the obstacle's bounding box
            d, g = sdf_grad(W[near], o['boxes'])
            h = d - m[near]
            act = h < F.ACT
            if not act.any():
                continue
            x, g, h = W[near][act], g[act], h[act]
            row = np.cross(x, g) @ Z.T - g @ ZxO.T
            row[:, L['n_joints']:] = 0.
            A.append(row); C.append(F.GAMMA * h)
    if not A:
        return np.zeros((0, F.DOF)), np.zeros(0)
    return np.vstack(A), np.concatenate(C)


def simulate(q0, qg, links, obst, prim, steps=1000):
    q = np.array(q0, float)
    log = dict(t=[], rows=[], gap=[], slack=0, reached=None, q=[q.copy()])
    z, E = None, np.zeros((1, F.DOF))
    for k in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        A, C = rows(q, links, obst, prim)
        u, s, z = S.solve(u_nom, A, np.zeros((len(C), 1)), C, E, F.BOX_G, F.BOX_H, z)
        log['t'].append(time.perf_counter() - t0)
        log['slack'] += s > 0
        log['rows'].append(len(C))
        log['gap'].append(F.real_gap(q, links, obst))
        q = q + F.DT * u
        log['q'].append(q.copy())
        if np.linalg.norm(qg - q) < .05:
            log['reached'] = (k + 1) * F.DT
            break
    log['gap'].append(F.real_gap(q, links, obst))
    return log


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--methods', nargs='+', default=['points_delta', 'capsule', 'spheres_enclosing'])
    p.add_argument('--edge-mm', type=float, default=10.)
    p.add_argument('--trials', nargs=2, type=int, default=[0, 30])
    p.add_argument('--steps', type=int, default=1000)
    p.add_argument('--spheres', type=int, default=16, help='centers per link for spheres_kmeans')
    p.add_argument('--tag', default='')
    a = p.parse_args()
    S.QP_SOLVER = 'daqp'
    links, obst, _ = F.build()
    trials = json.loads((Path(F.HERE) / 'franka_trials.json').read_text())
    for method in a.methods:
        prim = primitives(links, method, a.edge_mm / 1000, a.spheres)
        out = OUT / (method + (f'{a.spheres}' if method == 'spheres_kmeans' else '') + a.tag)
        out.mkdir(parents=True, exist_ok=True)
        info = dict(method=method, edge_mm=a.edge_mm, spheres_per_link=a.spheres if method == 'spheres_kmeans' else None,
                    primitives=int(sum(len(x[0]) for x in prim)),
                    margin_mm=[1e3 * float(x[1].max()) for x in prim])
        (out / 'setup.json').write_text(json.dumps(info, indent=1))
        print(info, flush=True)
        for i in range(*a.trials):
            path = out / f'franka_{i:02d}.json'
            if path.exists():
                continue
            q0, qg, _ = trials[i]
            log = simulate(q0, np.asarray(qg), links, obst, prim, a.steps)
            t, gaps = 1e3 * np.asarray(log['t']), np.asarray(log['gap'])
            m = dict(steps=len(t), reached_s=log['reached'], slack_steps=int(log['slack']),
                     min_gap_bound_mm=1e3 * float(gaps.min()), final_gap_mm=1e3 * float(gaps[-1]),
                     nonpositive_gap_states=int((gaps <= 0).sum()), t_median_ms=float(np.median(t)),
                     t_p95_ms=float(np.percentile(t, 95)), t_max_ms=float(t.max()), rows_max=int(max(log['rows'])))
            path.write_text(json.dumps(dict(trial=i, **info, metrics=m), indent=1))
            np.savez_compressed(out / f'franka_{i:02d}.npz', q=np.asarray(log['q']), t=t, gap=gaps,
                                rows=np.asarray(log['rows']))
            print(method, i, m, flush=True)


if __name__ == '__main__':
    main()
