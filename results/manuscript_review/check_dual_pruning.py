"""Recover valid lower bounds for historical selected-state dual-arm audits.

Historical sample-distance pruning is justified with an explicit mesh-to-sample
cover radius. Every skipped term has a finite floor. Clamp old minima by that
floor; no new trajectory or whole-trajectory claim is inferred.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
import dual3d as D


def main():
    links, _ = D.build(cover_side='cached')
    rmax = max(float(L['gt_r'].max()) for L in links)
    cover = max(L['sample_cover_radius'] for L in links)
    floor = min(.02, .03 - rmax, .03 - cover - rmax)
    assert floor > 0
    link_data = []
    for L in links:
        historical_radius = float(np.linalg.norm(L['gt'], axis=1).max()) + .02
        mesh_radius = float(np.linalg.norm(L['V'], axis=1).max())
        assert historical_radius >= mesh_radius
        f = L['fs']
        radius = L['rad'] + .06 + D.F.ACT
        index = np.stack(np.meshgrid(*(np.arange(k) for k in f.K), indexing='ij'), axis=-1)
        lo = f.lo + f.h * index
        farthest = np.linalg.norm(np.maximum(np.abs(lo), np.abs(lo + f.h)), axis=-1)
        # Include every cell INTERSECTING the exterior of the sphere, not only
        # cells whose nearest point is outside. The latter misses cut cells.
        exterior_cells = farthest >= radius
        coeff_min = f.cell_coeffs().min(axis=(3,4,5)) - L['l']
        field_floor = float(coeff_min[exterior_cells].min()) if exterior_cells.any() else None
        assert field_floor is None or field_floor > 0
        link_data.append(dict(link=L['frame'], sample_cover_radius_mm=1000*L['sample_cover_radius'],
                              historical_sphere_margin_mm=1000*(historical_radius-mesh_radius),
                              controller_pair_skip_radius_m=radius,
                              controller_pair_skip_field_bound_mm=None if field_floor is None else 1000*field_floor))
    source = ROOT / 'prototype_3d/dual_results.json'
    previous = json.loads(source.read_text())
    revised = [dict(trial=t['trial'], historical_selected_lower_bound_mm=t['certified_min_gap_mm'] if np.isfinite(t['certified_min_gap_mm']) else None,
                    audited_selected_lower_bound_mm=min(1000*floor, t['certified_min_gap_mm']))
               for t in previous['trials']]
    result = dict(source=str(source.relative_to(ROOT)), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  scope='Same selected states as historical records; no all-state or between-step claim.',
                  omitted_term_lower_bound_mm=1000*floor, links=link_data, trials=revised,
                  minimum_selected_bound_mm=min(t['audited_selected_lower_bound_mm'] for t in revised))
    assert all(t['audited_selected_lower_bound_mm'] > 0 for t in revised)
    # Check revised finite-bound routine at a recorded initial pose, without simulating.
    initial = np.asarray(previous['trials'][7]['q0'])
    result['initial_pose_check_mm'] = 1000*D.certified_gap(initial, links)
    assert 0 < result['initial_pose_check_mm'] <= 20
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
