"""Audit all saved star-tube poses using continuous triangle bounds and geometric pruning.

No fitted-field cutoff is used. For each mesh triangle, the 1-Lipschitz analytic
obstacle field at its centroid minus its radius bounds the entire triangle.
Bounding-box pruning only discards triangles whose distance is already >= cap.
Rigid-motion bounds transfer a previous link certificate to another saved pose;
the full triangle check is repeated whenever that transferred bound is too small.
This audits saved states, NOT motion between states, and assumes exact arithmetic.
"""
import argparse
import json
import os
from pathlib import Path
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import franka3d as F
import dolphin3d as D


def link_bound(T, L, RB, pB, cap=.05):
    ext = np.array([D.HALF_T, D.R_OUT, D.R_OUT])
    P = (L['gt'] @ T[:3, :3].T + T[:3, 3] - pB) @ RB
    # Distance to a containing AABB lower-bounds distance to the solid for
    # outside points. All inside points are kept, including deeply interior ones.
    outside = np.linalg.norm(np.maximum(np.abs(P) - ext, 0.), axis=1)
    use = outside - L['gt_r'] < cap
    if use.any():
        signed = D.hoop_sdf_rounded(P[use] + D.HOOP_C)
        return min(cap, float((signed - L['gt_r'][use]).min()))
    return cap


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('trajectory', type=Path)
    p.add_argument('--reuse-margin-mm', type=float, default=1.,
                   help='Recompute a link if its transferred lower bound is below this margin; inf disables reuse.')
    args = p.parse_args()
    setup = json.loads((Path(F.HERE) / 'interactive/setups/setup_20260928_171837.json').read_text())
    obj = setup['objects'][0]
    RB, pB = np.asarray(obj['R']), np.asarray(obj['p'])
    links, _, _ = F.build()
    qs = np.load(args.trajectory)['q']
    bounds = []
    # Every triangle point has norm at most the largest vertex norm by convexity.
    radii = [float(np.linalg.norm(L['V'], axis=1).max()) for L in links]
    cached = [None] * len(links)
    direct_checks, transferred = 0, 0
    for k, q in enumerate(qs):
        frames, _, _ = F.fk(q)
        per_link = []
        for j, L in enumerate(links):
            T = frames[L['frame']]
            candidate = -np.inf
            if cached[j] is not None:
                old_T, old_bound = cached[j]
                # ||(R-R0)x + t-t0|| <= ||R-R0||_2 max||x|| + ||t-t0||.
                # The analytic obstacle field is 1-Lipschitz, so subtracting
                # this displacement bounds every point of every triangle.
                displacement = (np.linalg.norm(T[:3, :3] - old_T[:3, :3], ord=2) * radii[j]
                                + np.linalg.norm(T[:3, 3] - old_T[:3, 3]))
                candidate = old_bound - displacement
            if candidate < args.reuse_margin_mm * .001:
                candidate = link_bound(T, L, RB, pB)
                cached[j] = (T.copy(), candidate)
                direct_checks += 1
            else:
                transferred += 1
            per_link.append(candidate)
        bounds.append(min(per_link))
        if k % 100 == 0:
            print('audited', k, 'min mm', 1e3 * min(bounds), flush=True)
    result = dict(trajectory=str(args.trajectory), states=len(qs), scope='all saved states, all modeled link mesh triangles; not inter-step motion',
                  triangle_max_edge_mm=4., min_bound_mm=1e3 * min(bounds),
                  unresolved_states=int(np.sum(np.asarray(bounds) <= 0)),
                  uses_fitted_field_pruning=False, direct_link_checks=direct_checks,
                  rigid_motion_transfers=transferred, reuse_margin_mm=args.reuse_margin_mm)
    target = args.trajectory.with_suffix('.audit.json')
    target.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
