"""Fairness item 3 of results/continuous_baselines/PLAN.md: does each method start inside its certified set?

    python prototype_3d/continuous_baselines_starts.py
For every baseline configuration, h0 = min over its rows of h at the start pose q0 (rows() returns gamma*h).
A trial with h0 < 0 starts outside the certified set, so the CBF guarantee does not apply to it. For ours, the
certified lower bound of the lifted barrier over the active boxes at q0 (summed3d.franka_rows) is reported.
Writes results/continuous_baselines/start_validity.json.
"""
import json
import os
from pathlib import Path

for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np                                    # noqa: E402

import continuous_baselines as CBL                    # noqa: E402
import franka3d as F                                  # noqa: E402
import summed3d as S3                                 # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'continuous_baselines' / 'start_validity.json'
CONFIGS = [  # (folder stem, method, obstacle, edge_mm, spheres)
    ('points_delta_fitted', 'points_delta', 'fitted', 10., 16),
    ('points_delta_fitted_e5', 'points_delta', 'fitted', 5., 16),
    ('points_delta_fitted_e15', 'points_delta', 'fitted', 15., 16),
    ('points_delta_fitted_e20', 'points_delta', 'fitted', 20., 16),
    ('capsule_fitted', 'capsule', 'fitted', 10., 16),
    ('spheres_enclosing_fitted', 'spheres_enclosing', 'fitted', 10., 16),
    ('spheres_kmeans16_fitted', 'spheres_kmeans', 'fitted', 10., 16),
    ('spheres_kmeans64_fitted', 'spheres_kmeans', 'fitted', 10., 64),
    ('spheres_kmeans128_fitted', 'spheres_kmeans', 'fitted', 10., 128),
    ('points_delta', 'points_delta', 'exact', 10., 16),
    ('capsule', 'capsule', 'exact', 10., 16),
    ('spheres_enclosing', 'spheres_enclosing', 'exact', 10., 16),
    ('spheres_kmeans16', 'spheres_kmeans', 'exact', 10., 16),
    ('spheres_kmeans64', 'spheres_kmeans', 'exact', 10., 64),
]


def main():
    trials = json.loads((Path(F.HERE) / 'franka_trials.json').read_text())
    q0s = [np.asarray(t[0], float) for t in trials]
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    links, obst, _ = F.build()
    for stem, method, obstacle, edge, k in CONFIGS:
        if stem in res:
            continue
        prim = CBL.primitives(links, method, edge / 1000, k)
        if obstacle == 'fitted':
            reach = max(float(x[1].max()) for x in prim)
            for o in obst.values():
                o['G_cb'] = F.neighborhood_bounds(o['f'], reach)[0]
        h0 = []
        for q in q0s:
            _, C = CBL.rows(q, links, obst, prim, obstacle)
            h0.append(float(C.min() / F.GAMMA) if len(C) else float('inf'))
        h0 = np.array(h0)
        res[stem] = dict(h0_mm=[1e3 * x for x in h0], outside_trials=[i for i in range(len(h0)) if h0[i] < 0],
                         n_outside=int((h0 < 0).sum()))
        print(stem, 'outside', res[stem]['n_outside'], res[stem]['outside_trials'], flush=True)
        OUT.write_text(json.dumps(res, indent=1))
    for stem, cache in (('ours_12mm', None),
                        ('ours_6mm', ROOT / 'results' / 'fine_native_6mm_trial' / 'cache_franka_6mm.npz')):
        if stem in res:
            continue
        L2, O2, _ = F.build(cache_path=cache) if cache else F.build()
        S3.prep_all(L2, O2.values())
        h0 = np.array([S3.franka_rows(q, L2, O2)[4] for q in q0s], float)
        res[stem] = dict(h0_lower_bound_mm=[1e3 * x for x in h0],
                         outside_trials=[i for i in range(len(h0)) if h0[i] < 0], n_outside=int((h0 < 0).sum()))
        print(stem, 'certified lower bound < 0 in', res[stem]['n_outside'], flush=True)
        OUT.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
