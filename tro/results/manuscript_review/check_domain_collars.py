"""Bernstein lower bounds on a 2r-wide domain collar for saved 3D fields.

A cover box with enclosing-ball radius <= r whose ball leaves a rectangular
domain has all its in-domain points within a 2r collar. Coordinate-clamped
extension puts its remaining points on the domain boundary. A strictly positive
collar bound above the certified level makes these omitted points harmless near
the zero barrier. This is a real-arithmetic Bernstein argument evaluated in
double precision, not an outward-rounded computation.
"""
import os
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[name] = '1'
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.path.join(str(ROOT), 'planar'))  # planar modules
import franka3d as F


def collar_bound(field, level, radius, dim=3):
    width = int(np.ceil(2 * radius / field.h))
    bound = np.inf
    for axis in range(dim):
        # Form only boundary strips, avoiding a full coefficient array for the maze.
        for start, stop in ((0, min(width, field.K[axis])),
                            (max(0, field.K[axis] - width), field.K[axis])):
            lo, count = field.lo.copy(), field.K.copy()
            lo[axis] += start * field.h
            count[axis] = stop - start
            index = [slice(None)] * 3
            index[axis] = slice(start, stop + 3)
            strip = F.spline_from(dict(lo=lo, h=field.h, K=count, W=field.W[tuple(index)]))
            bound = min(bound, float(strip.cell_coeffs().min()) - level)
    return dict(radius_m=radius, collar_width_m=2 * radius,
                whole_cell_collar_width_m=width * field.h,
                bernstein_lower_above_level_m=bound, positive=bool(bound > 0))


def main():
    old = np.load(F.CACHE, allow_pickle=True)['d'].item()
    out = {}
    for side in (.006, .012):
        r = side * np.sqrt(3) / 2
        for i, link in enumerate(old['links'], 1):
            key = f'link_{i}_cover_{side*1000:g}mm'
            out[key] = collar_bound(F.spline_from(link['f']), link['l'], r)
        for name, body in old['obst'].items():
            out[f'{name}_cover_{side*1000:g}mm'] = collar_bound(F.spline_from(body['f']), body['l'], r)
    fine = np.load(ROOT / 'results/fine_native_6mm_trial/cache_franka_6mm.npz', allow_pickle=True)['d'].item()
    for i, link in enumerate(fine['links'], 1):
        out[f'fine_link_{i}'] = collar_bound(F.spline_from(link['f']), link['l'], .006 * np.sqrt(3) / 2)
    objects = np.load(ROOT / 'prototype_3d/interactive/cache_objects.npz', allow_pickle=True)['d'].item()
    star = objects['star tube (rounded)']
    out['rounded_star'] = collar_bound(F.spline_from(star['f']), star['l'], .006 * np.sqrt(3) / 2)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    for name, data in out.items():
        print(name, f'{data["bernstein_lower_above_level_m"]*1000:.3f} mm', 'PASS' if data['positive'] else 'UNRESOLVED')
    assert all(data['positive'] for data in out.values())


if __name__ == '__main__':
    main()
