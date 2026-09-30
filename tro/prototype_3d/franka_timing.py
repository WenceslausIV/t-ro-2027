"""
Timing of the surface-cover filter for the paper, measured in isolation: replays the Franka trials
(franka_trials.json) and the dual-arm trials (dual_trials.json) with the same controller and records, per
step, the time to build the barrier rows (FK, pruning, spline evaluations) and the time of the QP.
Run with no other heavy process and one thread:

    set OMP_NUM_THREADS=1 & set MKL_NUM_THREADS=1 & python prototype_3d/franka_timing.py
"""
import json
import os
import platform
import sys
import time

for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(v, '1')

import numpy as np

import dual3d as D
import franka3d as F
from sdf_cbf_utils import solve_ldp_qp


def run(q0, qg, rows_fn, n_dof, box_G, box_H, steps=1000):
    q = np.array(q0, float)
    t_rows, t_qp, n_rows = [], [], []
    u = None
    for _ in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        G, h = rows_fn(q)
        t1 = time.perf_counter()
        u, _ = solve_ldp_qp(u_nom, np.vstack([G, box_G]), np.r_[-F.GAMMA * h, box_H], len(h), u)
        t2 = time.perf_counter()
        t_rows.append(t1 - t0); t_qp.append(t2 - t1); n_rows.append(len(h))
        q = q + F.DT * u
        if np.linalg.norm(qg - q) < .05:
            break
    return t_rows, t_qp, n_rows


def run_summed(q0, qg, rows_fn, n_dof, box_G, box_H, steps=1000):
    import summed as SM
    q = np.array(q0, float)
    t_rows, t_qp, n_rows = [], [], []
    z = None
    for _ in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        Ar, Tr, Cr, nb, hl, Er = rows_fn(q)
        t1 = time.perf_counter()
        u, _, z = SM.solve(u_nom, Ar, Tr, Cr, Er, box_G, box_H, z)
        t2 = time.perf_counter()
        t_rows.append(t1 - t0); t_qp.append(t2 - t1); n_rows.append(len(Cr))
        q = q + F.DT * u
        if np.linalg.norm(qg - q) < .05:
            break
    return t_rows, t_qp, n_rows


def stats(x):
    x = 1e3 * np.asarray(x)
    return dict(median=float(np.median(x)), p95=float(np.percentile(x, 95)), max=float(x.max()), mean=float(x.mean()))


def main():
    # python franka_timing.py [ours] [corner] [spheres] [points] [kkt] [dual]   (default: ours dual);
    # results are merged into timing_results.json ('franka' = ours = summed-field barriers)
    todo = sys.argv[1:] or ['ours', 'dual']
    path = os.path.join(F.HERE, 'timing_results_native.json')
    out = json.load(open(path)) if os.path.exists(path) else {}
    out.update(cpu=platform.processor(), python=platform.python_version(), threads=os.environ['OMP_NUM_THREADS'])
    import franka_baselines as Bl
    import summed3d as S3
    links, obst, _ = F.build()
    S3.prep_all(links, obst.values())
    trials = json.load(open(os.path.join(F.HERE, 'franka_trials.json')))
    prims = {k: Bl.primitives(links, k) for k in ('spheres', 'points')}
    methods = {'ours': lambda st: (lambda q: S3.franka_rows(q, links, obst)),
               'corner': lambda st: (lambda q: F.barrier_rows(q, links, obst)),
               'spheres': lambda st: (lambda q: Bl.rows(q, links, obst, prims['spheres'])),
               'points': lambda st: (lambda q: Bl.rows(q, links, obst, prims['points'])),
               'kkt': lambda st: (lambda q: Bl.kkt_rows(q, links, obst, st))}
    for name, make in methods.items():
        if name not in todo:
            continue
        R, Q, N = [], [], []
        for q0, qg, _ in trials:
            r, qp, n = (run_summed if name == 'ours' else run)(q0, qg, make({}), F.DOF, F.BOX_G, F.BOX_H)
            R += r; Q += qp; N += n
        tot = np.array(R) + np.array(Q)
        N = np.array(N)
        key = 'franka' if name == 'ours' else f'franka_{name}'
        out[key] = dict(steps=len(N), rows=stats(R), qp=stats(Q), total=stats(tot),
                        total_active=stats(tot[N > 0]) if (N > 0).any() else None, barriers_max=int(N.max()),
                        barriers_median_active=float(np.median(N[N > 0])) if (N > 0).any() else 0.,
                        frac_steps_active=float(np.mean(N > 0)))
        print(key, json.dumps(out[key]), flush=True)
        json.dump(out, open(path, 'w'), indent=1)
    if 'dual' not in todo:
        return
    dlinks, _ = D.build()
    S3.prep_all(dlinks)
    dtrials = json.load(open(os.path.join(F.HERE, 'dual_trials.json')))
    res = json.load(open(os.path.join(F.HERE, 'dual_results_native.json')))['trials']
    keep = {(tuple(np.round(t['q0'], 9)), tuple(np.round(t['qg'], 9))) for t in res}
    R, Q, N = [], [], []
    for q0, qg in dtrials:
        if (tuple(np.round(q0, 9)), tuple(np.round(qg, 9))) not in keep:
            continue
        r, qp, n = run_summed(q0, qg, lambda q: S3.dual_rows(q, dlinks), D.DOF2, D.BOX_G, D.BOX_H)
        R += r; Q += qp; N += n
    tot = np.array(R) + np.array(Q)
    N = np.array(N)
    out['dual'] = dict(steps=len(N), rows=stats(R), qp=stats(Q), total=stats(tot),
                       total_active=stats(tot[N > 0]), barriers_max=int(N.max()),
                       barriers_median_active=float(np.median(N[N > 0])),
                       frac_steps_active=float(np.mean(N > 0)))
    print('dual', json.dumps(out['dual']), flush=True)
    json.dump(out, open(path, 'w'), indent=1)


if __name__ == '__main__':
    main()
