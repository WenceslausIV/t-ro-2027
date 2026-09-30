"""Choose a Fig. 1 grid using a 1-cm sampled contour-gap diagnostic."""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from patch_resolution_check import exact_distance
from cspace_cbf_5robots import random_shape
from dock_cover import body_field, tight_level_2d
from franka3d import hessian_majorant


def field_for_figure(G, h):
    return body_field(G, h, .1, training_spacing=.005)


def main():
    out = Path(__file__).resolve().parents[1]/'results'/'overview_1cm'
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    shapes = [random_shape(rng) for _ in range(5)]
    records = []
    for index in (3, 0):
        G = shapes[index]
        for h in (.08, .04, .02, .01):
            f = field_for_figure(G, h)
            level = tight_level_2d(f, hessian_majorant(f), G)
            axes = [np.arange(f.lo[a], f.lo[a]+f.K[a]*h, .001) for a in (0, 1)]
            X, Y = np.meshgrid(*axes)
            P = np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]
            value = f.eval(P).reshape(X.shape)
            fig, ax = plt.subplots()
            contour = ax.contour(X, Y, value, levels=[level])
            curve = np.vstack([s for s in contour.allsegs[0] if len(s)])
            plt.close(fig)
            gaps = exact_distance(curve, G)
            accepted = bool(gaps.max() <= .01)
            item = dict(shape_index=index, h=h, patches=int(np.prod(f.K[:2])),
                        counts=f.K[:2].tolist(), level_mm=level*1000,
                        sampled_gap_max_mm=float(gaps.max()*1000),
                        sampled_gap_min_mm=float(gaps.min()*1000), accepted=accepted)
            records.append(item)
            print(json.dumps(item), flush=True)
            if accepted:
                break
        else:
            raise RuntimeError('No tested resolution met the diagnostic target')
    (out/'selection.json').write_text(json.dumps(dict(
        target='Sampled maximum gap from certified level contour to true polygon <= 10 mm',
        warning='Contour gaps are diagnostics, not a continuous 1-cm clearance certificate.',
        training_spacing_m=.005, contour_grid_spacing_m=.001, records=records), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
