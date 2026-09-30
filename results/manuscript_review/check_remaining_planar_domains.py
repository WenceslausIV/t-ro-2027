"""Check compact level sets for the rocket and square-domain Fig. 1 fit."""
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
from dock_cover import body_field, tight_level_2d, refine
from proto3d import Spline3
from cspace_cbf_5robots import random_shape
from rocket_geometry import OUTLINE, local
from scipy.spatial import cKDTree
from matplotlib.path import Path as PolygonPath


def main():
    results = {}
    polygon = local(OUTLINE)
    field = body_field(polygon, .01, .06)
    level = tight_level_2d(field, F.hessian_majorant(field), polygon)
    results['rocket'] = dict(collar_bound(field, level, .005 * np.sqrt(2) / 2, dim=2), level=level)
    rng = np.random.default_rng(3)
    polygon = [random_shape(rng) for _ in range(5)][3]
    h = .04
    midpoint = (polygon.min(0) + polygon.max(0)) / 2
    domain_side = np.ceil((np.ptp(polygon, axis=0).max() + .2) / h) * h
    lo = midpoint - domain_side / 2
    tree, shape = cKDTree(refine(polygon)), PolygonPath(polygon)
    def training_sdf(p):
        d = tree.query(p[:, :2])[0]
        return np.where(shape.contains_points(p[:, :2]), -d, d)
    field = Spline3(np.r_[lo, -1.5*h], np.r_[lo + domain_side - 1e-10, 1.5*h], h).fit(training_sdf, .005)
    level = tight_level_2d(field, F.hessian_majorant(field), polygon)
    results['fig1_square'] = dict(collar_bound(field, level, h * np.sqrt(2) / 2, dim=2), level=level)
    Path(__file__).with_suffix('.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))
    assert all(r['positive'] for r in results.values())


if __name__ == '__main__':
    main()
