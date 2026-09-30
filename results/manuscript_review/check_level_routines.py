"""Check strict level routines against analytic boundary maxima and saved link fields."""
import os
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[name] = '1'
import json
from pathlib import Path
import sys
import numpy as np
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.path.join(str(ROOT), 'planar'))  # planar modules
import franka3d as F
from dock_cover import tight_level_2d
from maze_cover import tight_level_curves


class Quadratic:
    h = .006

    def eval(self, points, order=0):
        v = np.sum(points ** 2, axis=1)
        return v if not order else (v, 2 * points)


class HessianBound:
    def eval(self, points):
        return np.full(len(points), 2.)


class NegativeQuadratic(Quadratic):
    def eval(self, points, order=0):
        v = -np.sum(points ** 2, axis=1)
        return v if not order else (v, -2 * points)


def main():
    f, bound, eps = Quadratic(), HessianBound(), 1e-5
    polygon = np.array([[-.01, -.01], [.01, -.01], [.01, .01], [-.01, .01]])
    vertices = np.c_[polygon, np.zeros(4)]
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    radius = .01
    results = {}
    results['polygon'] = tight_level_2d(f, bound, polygon, eps)
    results['mesh'] = F.tight_level_mesh(f, bound, vertices, faces, eps)
    results['implicit_sphere'] = F.tight_level_implicit(
        f, bound, lambda p: np.linalg.norm(p, axis=1) - radius,
        np.full(3, -radius), np.full(3, radius), eps=eps)
    for name, true_max in [('polygon', .0002), ('mesh', .0002), ('implicit_sphere', .0001)]:
        assert true_max < results[name] <= true_max + eps + 2e-11, (name, results[name])
    for invalid in (0., -1e-4, float('nan')):
        try:
            tight_level_2d(f, bound, polygon, invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid tolerance accepted')

    # Two abutting boxes have an interior zero seam in min(sdf_1,sdf_2).
    # Its field maximum exceeds the true boundary maximum, so a depth limit
    # must enclose the unresolved seam rather than discard it or loop forever.
    def box_sdf(points, center):
        q = np.abs(points - center) - np.array([.005, .01, .01])
        return np.linalg.norm(np.maximum(q, 0), axis=1) + np.minimum(q.max(axis=1), 0)

    union = lambda p: np.minimum(box_sdf(p, [-.005, 0, 0]), box_sdf(p, [.005, 0, 0]))
    level, info = F.tight_level_implicit(NegativeQuadratic(), bound, union,
        np.full(3, -.01), np.full(3, .01), eps=eps, h0=.004,
        max_depth=2, return_info=True)
    assert info['reason'] == 'depth_limit' and not info['epsilon_tight'], info
    assert info['lower_bound'] <= -.0001 + 1e-12 and level > -.0001, info
    results['abutting_union_conservative_cap'] = info
    level, info = F.tight_level_implicit(f, bound,
        lambda p: np.linalg.norm(p, axis=1) - radius,
        np.full(3, -radius), np.full(3, radius), eps=1e-12,
        max_boxes=1, return_info=True)
    assert info['reason'] == 'box_limit' and level > .0001, info
    results['sphere_conservative_box_cap'] = info

    class HeightField:
        h = .006
        def eval(self, points, order=0):
            v = points[:, 1]
            return v if not order else (v, np.tile([0., 1., 0.], (len(points), 1)))
    class ZeroHessian:
        def eval(self, points):
            return np.zeros(len(points))
    # Exact quadratic arch in cubic Bezier form; maximum y is 0.01 at t=0.5.
    arch = np.array([[[-.01, 0], [-.01/3, .04/3], [.01/3, .04/3], [.01, 0]]])
    wall = SimpleNamespace(bez=arch, chord_error=.002)
    level, chord = tight_level_curves(HeightField(), ZeroHessian(), [wall], eps)
    assert .01 < level <= .01 + eps + 1e-12 and chord == .002, level
    results['cubic_arch'] = level
    saved = np.load(F.CACHE, allow_pickle=True)['d'].item()
    results['historical_link_level_rechecks'] = []
    for i, link in enumerate(saved['links'], 1):
        field = F.spline_from(link['f'])
        new = F.tight_level_mesh(field, F.hessian_majorant(field), link['V'], link['F'])
        delta = new - float(link['l'])
        results['historical_link_level_rechecks'].append(dict(link=i, saved_level=link['l'],
                                                             recomputed_level=new, difference=delta))
        print(f'Link {i}: saved={link["l"]:.12g}, recomputed={new:.12g}, difference={delta:.3g}', flush=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print('Analytic strict-boundary checks passed; saved-field rechecks recorded.', flush=True)


if __name__ == '__main__':
    main()
