"""Fig. 1 B_i: 10x10 versus 119x90 SDF patches, fixed training data/domain.

Run: python prototype_3d/patch_resolution_check.py
Outputs are standalone diagnostics; no paper or simulation results are replaced.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import sys
import json
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cspace_cbf_5robots import random_shape
from dock_cover import refine, seg_distance, tight_level_2d, body_field
from proto3d import Spline3, _b
from sdf_cbf_utils import Q_B2B


class PhysicalField:
    """Spline3 on normalized knot coordinates, with physical chain rules."""
    def __init__(self, base, lo, step):
        self.base, self.lo, self.step = base, lo, step

    def eval(self, points, order=0):
        ans = self.base.eval((np.asarray(points) - self.lo) / self.step, order)
        if order == 0:
            return ans
        if order == 1:
            return ans[0], ans[1] / self.step
        return ans[0], ans[1] / self.step, ans[2] / self.step[None, :, None] / self.step[None, None, :]

    def planar_hessian_bound(self):
        # Exact restriction to z=0, then Bernstein bounds in physical x,y.
        z = -self.lo[2] / self.step[2]
        iz = int(np.floor(z))
        weights = np.einsum('ijk,k->ij', self.base.W[:, :, iz:iz+4], _b(z-iz))
        ix = np.arange(self.base.K[0])[:, None] + np.arange(4)
        iy = np.arange(self.base.K[1])[:, None] + np.arange(4)
        local = weights[ix[:, None, :, None], iy[None, :, None, :]]
        coeff = np.einsum('ip,jq,abpq->abij', Q_B2B, Q_B2B, local, optimize=True)
        xx = np.max(np.abs(6*np.diff(coeff, n=2, axis=2))) / self.step[0]**2
        yy = np.max(np.abs(6*np.diff(coeff, n=2, axis=3))) / self.step[1]**2
        xy = np.max(np.abs(9*np.diff(np.diff(coeff, axis=2), axis=3))) / np.prod(self.step[:2])
        return float(np.sqrt(xx*xx + yy*yy + 2*xy*xy))


class ConstantBound:
    def __init__(self, value):
        self.value = value

    def eval(self, points):
        return np.full(len(points), self.value)


def fit(G, lo, hi, counts):
    start = time.perf_counter()
    tree, path = cKDTree(refine(G)), PolygonPath(G)
    # Exactly the original 5-mm sample spacing, including its thin z layer.
    axes = [np.arange(lo[a], hi[a]+1e-9, .005) for a in range(3)]
    xyz = np.meshgrid(*axes, indexing='ij')
    points = np.stack([v.ravel() for v in xyz], axis=1)
    distance = tree.query(points[:, :2])[0]
    targets = np.where(path.contains_points(points[:, :2]), -distance, distance).reshape(xyz[0].shape)
    data_s = time.perf_counter() - start
    solve_start = time.perf_counter()
    counts = np.r_[counts, 3]
    step = (hi-lo) / counts
    base = Spline3(np.zeros(3), counts, 1.)
    maps = []
    for a in range(3):
        design = base._design(a, (axes[a]-lo[a])/step[a])
        maps.append(np.linalg.solve(design.T @ design + 1e-4*np.eye(base.n[a]), design.T))
    base.W = np.einsum('ai,bj,ck,ijk->abc', *maps, targets, optimize=True)
    solve_s = time.perf_counter() - solve_start
    return PhysicalField(base, lo, step), dict(data_s=data_s, coefficient_s=solve_s, total_s=data_s+solve_s)


def exact_distance(points, G):
    values = []
    for offset in range(0, len(points), 512):
        values.append(seg_distance(points[offset:offset+512], G, np.roll(G, -1, axis=0)))
    d = np.concatenate(values)
    return np.where(PolygonPath(G).contains_points(points), -d, d)


def main():
    out = Path(__file__).resolve().parents[1] / 'results' / 'patch_resolution_100'
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    G = [random_shape(rng) for _ in range(5)][3]
    lo2, requested_hi = G.min(axis=0)-.1, G.max(axis=0)+.1
    original_counts = np.ceil((requested_hi-lo2)/.01).astype(int)
    lo, hi = np.r_[lo2, -.015], np.r_[lo2+original_counts*.01, .015]
    results, fields = [], []
    for counts in (original_counts, np.array([10, 10])):
        timings = []
        for _ in range(3):
            field, timing = fit(G, lo, hi, counts)
            timings.append(timing)
        start = time.perf_counter()
        bound = field.planar_hessian_bound()
        level = tight_level_2d(field, ConstantBound(bound), G)
        cert_s = time.perf_counter()-start
        result = dict(counts=counts.tolist(), patches=int(np.prod(counts)),
                      patch_mm=(field.step[:2]*1000).tolist(),
                      effective_planar_controls=int(np.prod(counts+3)),
                      actual_stored_controls=int(field.base.W.size), timing_runs=timings,
                      median_total_s=float(np.median([t['total_s'] for t in timings])),
                      median_coefficient_s=float(np.median([t['coefficient_s'] for t in timings])),
                      certified_level_mm=level*1000, certification_s=cert_s,
                      global_planar_hessian_bound=bound)
        fields.append(field)
        results.append(result)
        print(json.dumps(result), flush=True)
    # Check that the normalized fine fit reproduces the production implementation.
    native = body_field(G, .01, .1)
    check_rng = np.random.default_rng(419)
    probe = check_rng.uniform(lo[:2], hi[:2], size=(1000, 2))
    p3 = np.c_[probe, np.zeros(len(probe))]
    reproduction = float(np.max(np.abs(native.eval(p3)-fields[0].eval(p3))))
    assert reproduction < 1e-10, reproduction
    xs = np.linspace(lo[0], hi[0], 481)
    ys = np.linspace(lo[1], hi[1], 361)
    X, Y = np.meshgrid(xs, ys)
    points = np.c_[X.ravel(), Y.ravel()]
    print('Evaluating exact polygon distance on independent diagnostic grid...', flush=True)
    truth = exact_distance(points, G)
    near = np.abs(truth) <= .03
    p3 = np.c_[points, np.zeros(len(points))]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)
    for ax, field, result in zip(axes, fields, results):
        value = field.eval(p3)
        err = value-truth
        result.update(grid_rmse_mm=float(np.sqrt(np.mean(err**2))*1000),
                      grid_max_abs_error_mm=float(np.max(np.abs(err))*1000),
                      near_surface_rmse_mm=float(np.sqrt(np.mean(err[near]**2))*1000),
                      near_surface_max_abs_error_mm=float(np.max(np.abs(err[near]))*1000))
        ax.fill(G[:, 0], G[:, 1], color='.85', label='True polygon')
        ax.plot(*np.vstack([G, G[:1]]).T, color='.25', lw=1)
        ax.contour(X, Y, value.reshape(X.shape), levels=[0], colors=['#d95f02'], linestyles='--', linewidths=1.5)
        contour = ax.contour(X, Y, value.reshape(X.shape), levels=[result['certified_level_mm']/1000], colors=['#0072b2'], linewidths=1.8)
        curve = np.vstack([seg for seg in contour.allsegs[0] if len(seg)])
        gap = exact_distance(curve, G)
        result['certified_contour_sampled_gap_mm'] = dict(min=float(gap.min()*1000), median=float(np.median(gap)*1000), max=float(gap.max()*1000))
        result['contour_components'] = len(contour.allsegs[0])
        # Diagnostics only: contour interpolation and grid errors are not certificates.
        for a in range(2):
            knots = np.linspace(lo[a], hi[a], result['counts'][a]+1)
            for knot in knots:
                (ax.axvline if a == 0 else ax.axhline)(knot, color='.65', alpha=.25, lw=.45)
        ax.plot([], [], '--', color='#d95f02', label='Fitted zero contour')
        ax.plot([], [], color='#0072b2', label='Certified level contour')
        ax.set_title(f"{result['patches']:,} patches ({result['counts'][0]} x {result['counts'][1]})\nCertified level = {result['certified_level_mm']:.3f} mm")
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]')
        ax.legend(loc='upper right', fontsize=8)
    axes[0].set_ylabel('y [m]')
    fig.suptitle('Same Fig. 1 body, domain, training data, cubic degree and regularization')
    fig.tight_layout()
    fig.savefig(out/'comparison.png', dpi=180)
    fig.savefig(out/'comparison.pdf')
    report = dict(shape='five_robot_setup(seed=3, scale=1.0), GT[3]',
                  domain_lo=lo[:2].tolist(), domain_hi=hi[:2].tolist(),
                  training_spacing_m=.005, level_tolerance_m=.0001,
                  diagnostic_grid=[len(xs), len(ys)], near_surface_band_m=.03,
                  fine_reproduction_max_error_m=reproduction,
                  notes=['Timing: imports excluded, one BLAS thread, background activity not excluded.',
                         'Boundary levels certified by edge branch-and-bound with Bernstein Hessian bounds in double precision.',
                         'Distance errors and contour gaps are sampled diagnostics, not continuous certificates.',
                         'No controller simulation performed.'], results=results)
    (out/'results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
