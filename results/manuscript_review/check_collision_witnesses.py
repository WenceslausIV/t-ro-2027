"""Separate negative distance lower bounds from observed mesh/obstacle intersections.

Replays only historical baseline trials whose negative bound does not already imply a negative
centroid SDF through the global triangle-radius bound. Keeps the historical 6-mm cover used by
the KKT initializer. Writes only review artifacts; does not overwrite baseline results.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.path.join(str(ROOT), 'planar'))  # planar modules
import franka3d as F
import franka_baselines as B
from proto3d import sdf_boxes

OUT = Path(__file__).resolve().parent


def main():
    links, obst, _ = F.build(cover_side='cached')
    radius = max(float(L['gt_r'].max()) for L in links)
    old = json.loads((ROOT / 'prototype_3d/franka_baselines.json').read_text())
    trials = json.loads((ROOT / 'prototype_3d/franka_results.json').read_text())['trials']
    results = dict(max_triangle_radius_mm=1000*radius, methods={})
    nominal = [t['nominal_min_gap_mm'] for t in trials]
    results['nominal_collisions_implied_by_centroid_witness'] = sum(g < -1000*radius for g in nominal)
    for method in ('spheres_m0', 'points_m0', 'kkt'):
        kind = method.split('_')[0]
        B.MARGIN = old[method]['margin']
        prim = 'kkt' if kind == 'kkt' else B.primitives(links, kind)
        output = results['methods'][method] = []
        for previous in old[method]['trials']:
            i, bound = previous['trial'], previous['min_gap_mm']
            rec = dict(trial=i, historical_lower_bound_mm=bound)
            if bound >= 0:
                rec['status'] = 'positive_historical_saved_state_bound'
            elif bound < -1000*radius:
                rec.update(status='collision_implied_by_centroid',
                           witness_sdf_upper_mm=bound+1000*radius)
            else:
                print(f'Replaying ambiguous {method} trial {i}, bound {bound:.6g} mm', flush=True)
                states, sampled, bounds = [], [], []

                def gap_and_witness(q, link_data, objects, return_witness=False):
                    transforms, _, _ = F.fk(q)
                    lower, witness = np.inf, np.inf
                    for L in link_data:
                        T = transforms[L['frame']]
                        p = T[:3, 3]
                        near = [o for o in objects.values() if
                            np.all(p + L['rad'] > o['blo'] - .05) and
                            np.all(p - L['rad'] < o['bhi'] + .05)]
                        if not near:
                            continue
                        W = L['gt'] @ T[:3, :3].T + p
                        for o in near:
                            keep = np.all((W > o['blo']-.05) & (W < o['bhi']+.05), axis=1)
                            if keep.any():
                                d = sdf_boxes(W[keep], o['boxes'])
                                witness = min(witness, float(d.min()))
                                lower = min(lower, float((d-L['gt_r'][keep]).min()))
                    states.append(q.copy()); sampled.append(witness); bounds.append(lower)
                    return (lower, witness) if return_witness else lower

                original = F.real_gap
                F.real_gap = gap_and_witness
                start = time.monotonic()
                try:
                    log = B.simulate(np.array(trials[i]['q0']), np.array(trials[i]['qg']), links, obst, prim)
                finally:
                    F.real_gap = original
                rec.update(status='collision_witness' if min(sampled) < 0 else 'unresolved',
                           replay_lower_bound_mm=1000*min(bounds),
                           bound_difference_mm=1000*min(bounds)-bound,
                           minimum_centroid_sdf_mm=1000*min(sampled),
                           reached=log['reached'], slack_steps=log['slack'], elapsed_s=time.monotonic()-start)
                np.savez_compressed(OUT/f'{method}_trial_{i:02d}_witness.npz',
                                    q=states, lower=bounds, centroid_sdf=sampled)
                print(method, i, rec, flush=True)
            output.append(rec)
            (OUT/'check_collision_witnesses.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
        print(method, 'confirmed', sum(r['status'].startswith('collision') for r in output),
              'unresolved', sum(r['status']=='unresolved' for r in output), flush=True)


if __name__ == '__main__':
    main()
