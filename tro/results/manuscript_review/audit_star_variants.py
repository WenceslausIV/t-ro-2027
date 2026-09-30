"""Audit every saved star certificate-comparison state without fitted-field pruning."""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'prototype_3d'))
import franka3d as F
from certificate_audit_star import link_bound


def main():
    setup = json.loads((ROOT / 'prototype_3d/interactive/setups/setup_20260928_171837.json').read_text())
    obj = setup['objects'][0]
    RB, pB = np.asarray(obj['R']), np.asarray(obj['p'])
    links, _, _ = F.build(cover_side='cached')
    radii = [float(np.linalg.norm(L['V'], axis=1).max()) for L in links]
    results = {}
    for path in sorted((ROOT / 'results/certificate_ablation/pilot').glob('star_*.npz')):
        start = time.monotonic()
        states = np.load(path)['q']
        cached = [None] * len(links)
        bounds = []
        direct = transferred = 0
        for k, q in enumerate(states):
            frames, _, _ = F.fk(q)
            per_link = []
            for j, L in enumerate(links):
                T = frames[L['frame']]
                candidate = -np.inf
                if cached[j] is not None:
                    old_T, old_bound = cached[j]
                    movement = np.linalg.norm(T[:3, :3]-old_T[:3, :3], ord=2)*radii[j] + np.linalg.norm(T[:3, 3]-old_T[:3, 3])
                    candidate = old_bound - movement
                if candidate < .001:
                    candidate = link_bound(T, L, RB, pB)
                    cached[j] = (T.copy(), candidate)
                    direct += 1
                else:
                    transferred += 1
                per_link.append(candidate)
            bounds.append(min(per_link))
            if k % 100 == 0:
                print(path.stem, k, '/', len(states)-1, 'min mm', 1000*min(bounds), flush=True)
        np.save(OUT / (path.stem + '_mesh_bounds.npy'), bounds)
        results[path.stem] = dict(source=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            states=len(states), min_bound_mm=1000*min(bounds), unresolved_states=int(np.count_nonzero(np.asarray(bounds)<=0)),
            direct_link_checks=direct, transferred_bounds=transferred, reuse_threshold_mm=1.,
            scope='All saved states including final state; all modeled mesh triangles; no fitted-field pruning; no between-step certificate; double precision.',
            elapsed_s=time.monotonic()-start)
        (OUT / 'star_comparison_audit.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
        print(results[path.stem], flush=True)
    assert len(results)==6 and all(r['unresolved_states']==0 for r in results.values())


if __name__=='__main__':
    main()
