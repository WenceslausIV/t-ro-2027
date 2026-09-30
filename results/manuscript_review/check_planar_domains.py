"""Reconstruct planar fields without running controllers or overwriting historical results."""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'prototype_3d'))
from check_domain_collars import collar_bound
import franka3d as F
from dock_cover import body_field, tight_level_2d
from dock_shapes import station_body, craft_body
from cspace_cbf_5robots import random_shape
import maze_cover as M

OUT = Path(__file__).resolve().parent
results = {}


def check(name, f, level, side, provenance):
    results[name] = dict(collar_bound(f, level, side * np.sqrt(2) / 2, dim=2),
                         level=level, source=provenance)
    print(name, f'{1000*results[name]["bernstein_lower_above_level_m"]:.3f} mm', flush=True)
    (OUT / 'check_planar_domains.json').write_text(json.dumps(results, indent=2), encoding='utf-8')


def main():
    for name, polygon, h, margin, side in [('dock_station', station_body(), .0125, .15, .01),
                                         ('dock_craft', craft_body(), .01, .06, .01)]:
        f = body_field(polygon, h, margin)
        level = tight_level_2d(f, F.hessian_majorant(f), polygon)
        check(name, f, level, side, 'reconstructed with original polygon/fit settings')
    rng = np.random.default_rng(3)
    shapes = [random_shape(rng) for _ in range(5)]
    for i in (0, 3):
        f = body_field(shapes[i], .04, .1, training_spacing=.005)
        check(f'overview_{i}', f, tight_level_2d(f, F.hessian_majorant(f), shapes[i]), .04,
              'reconstructed Fig. 2 fields')
    stats_gt = np.load(ROOT / 'cache/shapes5.npz')['GT']
    for scene, polygons in [('five_stats', stats_gt), ('five_swap', [2*g for g in shapes])]:
        for i, polygon in enumerate(polygons):
            f = body_field(polygon, .01, .1)
            check(f'{scene}_{i}', f, tight_level_2d(f, F.hessian_majorant(f), polygon), .005,
                  'reconstructed historical field; no trajectory changed')

    M.M.FRAME = M.M.FRAME.copy()
    M.M.FRAME[:, 1] += M.M.PASS_Y[1] - M.M.FRAME[:, 1].mean()
    walls, path = M.M.designed_corridor_walls()
    waypoints = M.M.few_waypoints(path, walls, cut=M.M.WAYPOINT_CUT)
    lo = np.minimum(M.M.FRAME[0], waypoints[:, :2].min(0)) - .5
    hi = np.maximum(M.M.FRAME[1], waypoints[:, :2].max(0)) + .5
    print('Reconstructing symmetric-maze field for collar audit...', flush=True)
    f = M.wall_field(walls, lo, hi, .0125)
    old = json.loads((ROOT / 'results/reference_maze_summed_subdivided_symmetric_ports/results.json').read_text())
    old_level = old['field']['l_wall_mm'] / 1000
    check('symmetric_maze_walls', f, old_level, .005, 'reconstructed field with saved historical level')
    np.savez_compressed(OUT / 'maze_reconstructed_field.npz', d=np.array(F.spline_state(f), dtype=object))
    # Check the saved historical enclosure with direct on-curve lower bounds and shrinking Bezier hulls.
    _, hc = M.cell_bounds(f)
    new_level, chord = M.tight_level_curves(f, M.majorant_from(f, hc), walls)
    (OUT / 'maze_level_recheck.json').write_text(json.dumps(dict(
        historical_level=old_level, direct_cubic_level=new_level, chord_error=chord,
        historical_level_at_least_recomputed=bool(old_level >= new_level)), indent=2))
    print('Maze old/new certified levels:', old_level, new_level, flush=True)
    assert all(data['positive'] for data in results.values()), results


if __name__ == '__main__':
    main()
