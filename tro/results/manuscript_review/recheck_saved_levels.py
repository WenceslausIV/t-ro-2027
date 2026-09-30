"""Independently check old implicit/curve levels with repaired lower-bound routines.

Uses tighter tolerance than the historical experiments; does not replace their fields or levels.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.path.join(str(ROOT), 'planar'))  # planar modules
sys.path.insert(0, str(ROOT / 'prototype_3d'))
sys.path.insert(0, str(ROOT / 'prototype_3d/interactive'))
import franka3d as F
from proto3d import sdf_boxes
import server as S
import maze_cover as M

OUT = Path(__file__).resolve().parent


def main():
    results = {}
    old = np.load(F.CACHE, allow_pickle=True)['d'].item()
    bodies = [(name, o, lambda p, boxes=F.OBSTACLES[name]: sdf_boxes(p, boxes))
              for name, o in old['obst'].items()]
    star = np.load(ROOT/'prototype_3d/interactive/cache_objects.npz', allow_pickle=True)['d'].item()['star tube (rounded)']
    bodies.append(('rounded_star', star, S.shape_spec('star tube (rounded)')['sdf']))
    for name, body, distance in bodies:
        f = F.spline_from(body['f'])
        level, info = F.tight_level_implicit(f, F.hessian_majorant(f), distance,
            body['blo'], body['bhi'], eps=5e-5, return_info=True)
        info.update(historical_level=body['l'], old_strict_enclosure_reconfirmed=bool(level <= body['l']))
        results[name] = info
        print(name, info, flush=True)
        (OUT/'recheck_saved_levels.json').write_text(json.dumps(results,indent=2))

    M.M.FRAME = M.M.FRAME.copy()
    M.M.FRAME[:, 1] += M.M.PASS_Y[1] - M.M.FRAME[:, 1].mean()
    walls, _ = M.M.designed_corridor_walls()
    f = F.spline_from(np.load(OUT/'maze_reconstructed_field.npz',allow_pickle=True)['d'].item())
    _, hc = M.cell_bounds(f)
    level, chord = M.tight_level_curves(f, M.majorant_from(f, hc), walls, eps=1e-4)
    saved = json.loads((ROOT/'results/reference_maze_summed_subdivided_symmetric_ports/results.json').read_text())
    original = saved['field']['l_wall_mm']/1000
    results['symmetric_maze'] = dict(historical_level=original, recomputed_level=level,
        old_strict_enclosure_reconfirmed=bool(level <= original), epsilon=1e-4, chord_error=chord)
    print('symmetric_maze',results['symmetric_maze'],flush=True)
    (OUT/'recheck_saved_levels.json').write_text(json.dumps(results,indent=2))
    assert all(v['old_strict_enclosure_reconfirmed'] for v in results.values()), results


if __name__ == '__main__':
    main()
