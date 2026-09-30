"""
Curved maze of Sec. VIII-E (reference_maze_escape.py) with the SURFACE-COVER construction instead of the
tabulated configuration-space field: one workspace field for the two wall regions and one for the rocket,
tightest levels certified on the ground truths (Proposition 1; the walls' ground truth is their cubic
Bezier spans, handled through the certified chord bound), a Bernstein cover of the rocket's certified
curve, and corner barriers against the wall field with the C^2 curvature majorant and certified cluster
pruning (prototype_3d/dock_cover.make_cover_fn). Same passage, waypoints, nominal controller, gains,
step, and interval audit as the tabulated run.

    python planar/maze_cover.py [--redraw] [--side MM] [--wall-cell MM] [--robot-cell MM]

Results: results/reference_maze_cover/{results.json, trajectory.npz, field.npz}; with the default
rocket_maze_cover.png there (copied to tro/figs/rocket_maze.png by hand if chosen for the paper).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
import franka3d as F                                                    # noqa: E402
from proto3d import Spline3                                             # noqa: E402
from dock_cover import cover_2d, make_cover_fn, make_summed_fn, refine, body_field, ACT # noqa: E402
import summed as SM                                                     # noqa: E402
import reference_maze_escape as M                                       # noqa: E402
import maze_reference_geometry as G                                     # noqa: E402
import rocket_geometry as Rocket                                        # noqa: E402
from maze_escape import simulate                                        # noqa: E402

OUT = ROOT / 'results' / 'reference_maze_cover'


# ---------------------------------------------------------------------------------------------
# wall field (large: cellwise bounds are computed in slabs)
# ---------------------------------------------------------------------------------------------
def wall_field(walls, lo, hi, h):
    """2D cubic B-spline (thin constant z layer, as dock_cover.body_field) fitted to the signed distance
    of the union of the wall regions (negative inside), sampled on a grid of spacing h/2."""
    from skimage.draw import polygon as raster
    f = Spline3(np.r_[lo, -1.5 * h], np.r_[hi, 1.5 * h], h)
    axes = [np.arange(f.lo[a], f.lo[a] + f.K[a] * h + 1e-9, h / 2) for a in range(3)]
    X, Y = np.meshgrid(axes[0], axes[1], indexing='ij')
    P = np.c_[X.ravel(), Y.ravel()]
    tree = cKDTree(np.vstack([refine(w.dense, .001) for w in walls]))
    d = tree.query(P, workers=-1)[0]
    inside = np.zeros(X.shape, bool)                     # grid points inside a wall (sign of the targets)
    for w in walls:
        rr, cc = raster((w.dense[:, 0] - axes[0][0]) / (h / 2), (w.dense[:, 1] - axes[1][0]) / (h / 2), X.shape)
        inside[rr, cc] = True
    T2 = np.where(inside, -d.reshape(X.shape), d.reshape(X.shape))
    maps = []
    for a in range(3):
        D = f._design(a, axes[a])
        maps.append(np.linalg.solve(D.T @ D + 1e-4 * np.eye(f.n[a]), D.T))
    # the targets are constant in z, so the z map reduces to its row sums
    f.W = np.einsum('ai,bj,ij->ab', maps[0], maps[1], T2, optimize=True)[:, :, None] * maps[2].sum(axis=1)[None, None]
    return f


def slab(f, a0, a1):
    s = Spline3.__new__(Spline3)
    s.lo, s.h, s.W = f.lo, f.h, f.W[a0:a1 + 3]
    s.K = np.array([a1 - a0, f.K[1], f.K[2]])
    s.n = s.K + 3
    return s


def cell_bounds(f, width=48):
    """Per-cell certified bounds of |grad phi| and of the Hessian norm, x-slab by x-slab."""
    Gs, Hs = [], []
    for a0 in range(0, f.K[0], width):
        s = slab(f, a0, min(a0 + width, f.K[0]))
        C = s.cell_coeffs()
        Gs.append(np.sqrt(sum((3 * np.abs(np.diff(C, axis=3 + a)).max(axis=(3, 4, 5)) / f.h) ** 2
                              for a in range(3))))
        Hs.append(s.hessian_bound(cellwise=True))
    return np.concatenate(Gs), np.concatenate(Hs)


def third_bounds(f, width=48):
    """summed.third_cellwise, x-slab by x-slab (the full cell-coefficient array does not fit in memory)."""
    out = []
    for a0 in range(0, f.K[0], width):
        sl = slab(f, a0, min(a0 + width, f.K[0]))
        out.append(SM.third_cellwise(sl, sl.cell_coeffs()))
    return np.concatenate(out)


def majorant_from(f, Hcell):
    """franka3d.hessian_majorant with a precomputed cellwise Hessian bound."""
    W = np.pad(Hcell, 4, mode='edge')
    for a in range(3):
        W = np.lib.stride_tricks.sliding_window_view(W, 6, axis=a).max(axis=-1)
    g = Spline3.__new__(Spline3)
    g.lo, g.h, g.K, g.n, g.W = f.lo, f.h, f.K, f.n, W
    return g


def tight_level_curves(f, Mf, walls, eps=3e-4):
    """Certify the actual cubic boundary, using its Bezier control hulls and on-curve evaluations.

    Unlike chord-midpoint lower bounds, all values entering LB are attained on the true curve.
    Subdivision shrinks the hull radii to zero; the plotting polygon's chord error is returned
    for the independent distance audit only. Historical cached experiments are not rewritten.
    """
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('eps must be finite and positive')
    z = lambda P: np.c_[P, np.zeros(len(P))]
    spans = np.concatenate([w.bez for w in walls])

    def split(T):
        a = (T[:, :-1] + T[:, 1:]) / 2
        b = (a[:, :-1] + a[:, 1:]) / 2
        c = b.mean(axis=1)
        return np.concatenate((np.stack((T[:, 0], a[:, 0], b[:, 0], c), axis=1),
                               np.stack((c, b[:, 1], a[:, 2], T[:, 3]), axis=1)))

    LB = float(f.eval(z(spans[:, [0, 3]].reshape(-1, 2))).max())
    while len(spans):
        center = np.einsum('i,kid->kd', np.array([1, 3, 3, 1]) / 8, spans)
        radius = np.linalg.norm(spans - center[:, None], axis=2).max(axis=1)
        supported = radius <= f.h if hasattr(f, 'h') else np.ones(len(spans), bool)
        pending = spans[~supported]
        T, center, radius = spans[supported], center[supported], radius[supported]
        if len(T):
            v, g = f.eval(z(center), order=1)
            LB = max(LB, float(v.max()))
            UB = v + np.linalg.norm(g[:, :2], axis=1) * radius + .5 * Mf.eval(z(center)) * radius ** 2
            pending = np.concatenate((pending, T[UB >= LB + eps]))
        spans = split(pending) if len(pending) else pending
    return LB + eps, max(w.chord_error for w in walls)


def tight_level_polygon(f, Mf, P, eps=1e-4):
    """dock_cover.tight_level_2d (polygon ground truth)."""
    from dock_cover import tight_level_2d
    return tight_level_2d(f, Mf, P, eps)


# ---------------------------------------------------------------------------------------------
def level_outline(f, level, res=.004):
    """Polygons of {phi = level} (drawing only)."""
    from skimage import measure
    xs = np.arange(f.lo[0], f.lo[0] + f.K[0] * f.h, res)
    ys = np.arange(f.lo[1], f.lo[1] + f.K[1] * f.h, res)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    P = np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]
    V = np.concatenate([f.eval(P[k:k + 200000]) for k in range(0, len(P), 200000)]).reshape(X.shape)
    out = []
    for c in measure.find_contours(V, level):
        out.append(np.c_[np.interp(c[:, 0], np.arange(len(xs)), xs), np.interp(c[:, 1], np.arange(len(ys)), ys)])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--redraw', action='store_true')
    ap.add_argument('--side', type=float, default=None)       # default: native robot SDF patches
    ap.add_argument('--wall-cell', type=float, default=12.5)   # [mm]
    ap.add_argument('--robot-cell', type=float, default=10.)   # [mm]
    ap.add_argument('--corner', action='store_true')           # previous corner barriers (default: summed-field)
    ap.add_argument('--symmetric-ports', action='store_true',
                    help='center the frame on the middle pass, with equal entrance/exit straight lengths')
    args = ap.parse_args()
    global OUT
    if not args.corner:
        OUT = ROOT / 'results' / ('reference_maze_summed_native' if args.side is None else 'reference_maze_summed_subdivided')
    if args.symmetric_ports:
        # Change the actual wall geometry, then refit, simulate, and audit it.
        # Keep all historical asymmetric-scene results in their original directories.
        M.FRAME = M.FRAME.copy()
        M.FRAME[:, 1] += M.PASS_Y[1] - M.FRAME[:, 1].mean()
        OUT = OUT.with_name(OUT.name + '_symmetric_ports')
    OUT.mkdir(parents=True, exist_ok=True)
    robot, gt, robot_meta = Rocket.setup()
    walls, path = M.designed_corridor_walls()
    waypoints = M.few_waypoints(path, walls, cut=M.WAYPOINT_CUT)
    print(f'{len(walls)} wall regions, {len(waypoints)} waypoints', flush=True)
    lo = np.minimum(M.FRAME[0], waypoints[:, :2].min(0)) - .5
    hi = np.maximum(M.FRAME[1], waypoints[:, :2].max(0)) + .5

    t0 = time.perf_counter()
    fW = wall_field(walls, lo, hi, args.wall_cell / 1e3)
    t_fit = time.perf_counter() - t0
    Gc, Hc = cell_bounds(fW)
    MW = majorant_from(fW, Hc)
    t0 = time.perf_counter()
    lW, chord = tight_level_curves(fW, MW, walls)
    t_lW = time.perf_counter() - t0
    fR = body_field(gt, args.robot_cell / 1e3, .06)
    MR = F.hessian_majorant(fR)
    lR = tight_level_polygon(fR, MR, gt)
    PC, sz = cover_2d(fR, lR, None if args.side is None else args.side / 1e3)
    r = sz * np.sqrt(2) / 2
    _, _, crad = F.clusters(PC, r)
    k = int(np.ceil(float(crad.max()) / fW.h)) + 1
    bounds = (maximum_filter(Gc, size=2 * k + 1), maximum_filter(Hc, size=2 * (k + 4) + 1))
    fn = make_cover_fn(fW, MW, lW, PC, sz, bounds=bounds)
    zero = np.zeros(3)

    def rows_fn(x):
        R6, hv = fn(zero, x)
        return R6[:, 3:6], hv                          # walls are static: only the rocket's inputs

    summed_fn = None
    if not args.corner:
        fW._M3 = SM.field_bound(fW, sz * np.sqrt(2) / 2, cell=third_bounds(fW))
        fW._M3_reach = sz * np.sqrt(2) / 2
        body = SM.prepare_body(fR, lR, PC, sz, 2)
        sfn = make_summed_fn(fW, MW, lW, body, bounds=bounds, umax=(0., 0., 0., .8, .8, 1.4))

        def summed_fn(x):
            res = sfn(zero, x)
            if res is None:
                return None
            A6, Tk, C6, nb, hl, E6 = res
            return A6[:, 3:6], Tk, C6, nb, hl, E6[:, 3:6]   # walls are static: only the rocket's inputs

    meta = dict(wall_cell_mm=args.wall_cell, robot_cell_mm=args.robot_cell, wall_cells=[int(k) for k in fW.K[:2]],
                l_wall_mm=1e3 * lW, l_robot_mm=1e3 * lR, wall_chord_bound_mm=1e3 * chord, n_boxes=len(PC),
                side_mm=1e3 * sz, fit_s=t_fit, level_wall_s=t_lW, act_m=ACT)
    meta.update(frame_m=M.FRAME.tolist(), symmetric_ports=args.symmetric_ports,
                straight_port_lengths_m=[float(M.FRAME[1, 1] - (M.PASS_Y[0] + M.TURN_R)),
                                        float((M.PASS_Y[2] - M.TURN_R) - M.FRAME[0, 1])])
    print(meta, flush=True)
    if args.redraw:
        arrays = dict(np.load(OUT / 'trajectory.npz'))
    else:
        arrays, result = simulate(None, None, gt, walls, waypoints=waypoints, domain=M.DOMAIN,
                                  gap_function=G.physical_gaps, output_dir=OUT, cruise=True, steps=18000,
                                  switch_radius=.40, dyn='si', rows_fn=None if summed_fn else rows_fn,
                                  summed_fn=summed_fn)
        result.update(construction='summed-field barriers' if summed_fn else 'surface cover (corners)',
                      field=meta, robot=robot_meta, wall_count=len(walls),
                      gap_definition='Polygon distance minus the certified chord bound of the wall cubics.')
        (OUT / 'results.json').write_text(json.dumps(result, indent=2, default=float))
        np.savez_compressed(OUT / 'trajectory.npz', **arrays)
    # figure in the style of the tabulated one; red / blue outlines: certified level curves of the two fields
    walls_lv = [type('W', (), dict(dense=c))() for c in level_outline(fW, lW, .01) if len(c) > 50]
    rob_lv = max(level_outline(fR, lR, .001), key=len)

    class Outline:                                      # rocket's certified curve, moved with the pose
        rho = robot.rho

        @staticmethod
        def world(pose):
            import cspace_sdf_cbf_compare as C
            return rob_lv @ C.rot(pose[2]).T + pose[:2]

    M.draw(Outline, gt, walls, arrays, fitted_walls=walls_lv, out_dir=OUT, fig_path=OUT / 'rocket_maze_cover.png')


if __name__ == '__main__':
    main()
