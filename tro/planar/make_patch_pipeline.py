"""Introductory figure from an actual fitted SDF and its native Bernstein cover.

Run: python planar/make_patch_pipeline.py
The dense display grid is for rendering only; cover selection uses spline coefficients.
"""
import os
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[name] = '1'

import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch
from matplotlib.path import Path as PolygonPath
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
from cspace_cbf_5robots import random_shape
from dock_cover import tight_level_2d, cover_2d, seg_distance, refine
from proto3d import Spline3
from scipy.spatial import cKDTree
import franka3d as F


def main():
    rng = np.random.default_rng(3)
    shape = [random_shape(rng) for _ in range(5)][3]
    h = .04
    midpoint = (shape.min(0) + shape.max(0)) / 2
    domain_side = np.ceil((np.ptp(shape, axis=0).max() + .2) / h) * h
    square_lo = midpoint - domain_side / 2
    # Fit on a genuinely square domain, not a square crop of a rectangular field.
    tree, polygon = cKDTree(refine(shape)), PolygonPath(shape)
    def training_sdf(points):
        distance = tree.query(points[:, :2])[0]
        return np.where(polygon.contains_points(points[:, :2]), -distance, distance)
    field = Spline3(np.r_[square_lo, -1.5*h],
                    np.r_[square_lo + domain_side - 1e-10, 1.5*h], h).fit(training_sdf, .005)
    assert field.K[0] == field.K[1]
    level = tight_level_2d(field, F.hessian_majorant(field), shape)
    centers, side = cover_2d(field, level)
    assert np.isclose(side, field.h)
    lo, hi = field.lo[:2], (field.lo + field.K * field.h)[:2]
    metadata_dir = ROOT / 'results' / 'patch_pipeline'
    metadata_dir.mkdir(parents=True, exist_ok=True)
    xs, ys = [np.linspace(a, b, 501) for a, b in zip(lo, hi)]
    X, Y = np.meshgrid(xs, ys)
    values = field.eval(np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]).reshape(X.shape)
    # Exact point-to-segment distances to the ground-truth polygon, not distances
    # to its vertices or to the spline. The cache only speeds up figure rendering.
    cache = metadata_dir / 'exact_polygon_sdf_square.npz'
    if cache.exists():
        with np.load(cache) as data:
            assert np.array_equal(data['shape'], shape)
            assert np.array_equal(data['xs'], xs) and np.array_equal(data['ys'], ys)
            truth = data['truth']
    else:
        points = np.c_[X.ravel(), Y.ravel()]
        distance = np.concatenate([seg_distance(chunk, shape, np.roll(shape, -1, axis=0))
                                   for chunk in np.array_split(points, 128)])
        truth = np.where(PolygonPath(shape).contains_points(points), -distance, distance).reshape(X.shape)
        np.savez_compressed(cache, shape=shape, xs=xs, ys=ys, truth=truth)
    blue, green, teal, black = '#126eaa', '#4b9275', '#25867d', '#252b31'
    plt.rcParams.update({'font.family': 'serif', 'font.size': 8,
                         'mathtext.fontset': 'dejavuserif', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.16, 2.65), facecolor='white')
    positions = [.005, .251, .497]
    axes = [fig.add_axes([left, .32, .223, .54]) for left in positions]
    for ax in axes:
        ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), aspect='equal')
        ax.set_axis_off()
    grid = [[(x, lo[1]), (x, hi[1])] for x in lo[0] + np.arange(field.K[0]+1)*side]
    grid += [[(lo[0], y), (hi[0], y)] for y in lo[1] + np.arange(field.K[1]+1)*side]
    for ax in axes:
        ax.add_patch(Rectangle(lo, *(hi-lo), facecolor='#f2f2f2', edgecolor='none', zorder=-3))
    for ax in axes[1:]:
        ax.add_collection(LineCollection(grid, colors='#c5c8cc', linewidths=.25, zorder=-1))
    level_spacing = .04
    first = int(np.floor(min(truth.min(), values.min()) / level_spacing))
    last = int(np.ceil(max(truth.max(), values.max()) / level_spacing))
    levels = level_spacing * np.array([k for k in range(first, last+1) if k != 0])
    for ax, data in zip(axes[:2], [truth, values]):
        ax.contour(X, Y, data, levels=levels, colors='#747a80',
                   linewidths=.45, linestyles='solid')
    axes[0].add_patch(Polygon(shape, fill=False, edgecolor=black, lw=1))
    axes[0].text(*shape.mean(0), r'$\mathcal{G}_i$', ha='center', va='center', fontsize=9)
    axes[1].add_patch(Polygon(shape, fill=False, edgecolor=black, lw=.55, linestyle='--'))
    axes[1].contour(X, Y, values, levels=[0], colors=[teal], linewidths=.75, linestyles='--')
    axes[1].contour(X, Y, values, levels=[level], colors=[blue], linewidths=1.15)

    patches = [Rectangle(c-side/2, side, side) for c in centers]
    axes[2].add_collection(PatchCollection(patches, facecolor='#d4e5dc', edgecolor='none', zorder=1))
    axes[2].contour(X, Y, values, levels=[level], colors=[blue], linewidths=1.15)
    axes[2].add_collection(PatchCollection(patches, facecolor='none', edgecolor=green,
                                          linewidths=.5, zorder=1.5))
    axes[2].text(*shape.mean(0), r'$\mathcal{B}_i$', ha='center', va='center', fontsize=10, color=blue, zorder=6)
    for x, title in zip(positions, ['(a) Ground-truth SDF', '(b) Patch-wise SDF', '(c) Reuse surface patches']):
        fig.text(x+.1115, .94, title, ha='center', va='center')
    for x, text in zip(positions, [r'Exact SDF of $\mathcal{G}_i$', r'Fit $\phi_i$; certify level $l_i$',
                                  r'$\mathcal{S}_i=\{\phi_i=l_i\}$']):
        fig.text(x+.1115, .28, text, ha='center', fontsize=7.5,
                 color=blue if x == positions[2] else '#343c44')
    for x in (.238, .484):
        fig.add_artist(FancyArrowPatch((x-.008, .54), (x+.012, .54),
                                      transform=fig.transFigure, arrowstyle='-|>',
                                      mutation_scale=9, color='#56616c', lw=.8))

    # Schematic statement of Proposition coef, not samples of the surface.
    fig.text(.86, .94, '(d) Finite control constraints', ha='center')
    fig.text(.86, .78, 'Bernstein + Taylor bounds', ha='center', fontsize=7)
    fig.text(.86, .65, r'$L_{kj}(\mathbf{u},\mathbf{a})\geq 0$', ha='center', fontsize=10)
    fig.text(.86, .52, r'$j=1,\ldots,4$ per patch', ha='center', fontsize=7.5)
    fig.text(.86, .38, r'$\Downarrow$', ha='center', fontsize=12, color=blue)
    fig.text(.86, .25, r'$\dot h(\mathbf{x})+\gamma h(\mathbf{x})\geq0$', ha='center', fontsize=8)
    fig.text(.86, .135, r'$\forall\,\mathbf{x}\in\mathcal{S}_i\cap\mathcal{Q}_k$', ha='center', fontsize=8, color=blue)
    fig.add_artist(FancyArrowPatch((.722, .54), (.747, .54), transform=fig.transFigure,
                                  arrowstyle='-|>', mutation_scale=9, color='#56616c', lw=.8))
    # A boundary detail makes millimetric fitting/enclosure differences visible
    # without changing the fitted field or exaggerating the separation of curves.
    surface_values = field.eval(np.c_[shape, np.zeros(len(shape))])
    center = shape[np.argmin(surface_values)]
    radius = .024
    zoom = fig.add_axes([.253, .035, .105, .205])
    zx, zy = [np.linspace(a-radius, a+radius, 301) for a in center]
    ZX, ZY = np.meshgrid(zx, zy)
    zv = field.eval(np.c_[ZX.ravel(), ZY.ravel(), np.zeros(ZX.size)]).reshape(ZX.shape)
    zoom.add_patch(Polygon(shape, facecolor='#e7e9eb', edgecolor=black, lw=.9))
    zoom.contour(ZX, ZY, zv, levels=[0], colors=[teal], linewidths=1, linestyles='--')
    zoom.contour(ZX, ZY, zv, levels=[level], colors=[blue], linewidths=1.2)
    zoom.set(xlim=(zx[0], zx[-1]), ylim=(zy[0], zy[-1]), aspect='equal', xticks=[], yticks=[])
    for spine in zoom.spines.values():
        spine.set_color('#888888')
        spine.set_linewidth(.5)
    axes[1].add_patch(Rectangle(center-radius, 2*radius, 2*radius, fill=False,
                                 edgecolor=black, linewidth=.65, zorder=6))
    fig.text(.3055, .248, '(b) boundary detail', ha='center', fontsize=6)
    fig.text(.116, .16, 'Uniformly spaced level sets', ha='center', fontsize=7)
    fig.text(.116, .09, r'$0,\ \pm n,\ \pm 2n,\ \pm 3n,\ldots$', ha='center', fontsize=7)
    handles = [Line2D([], [], color=black, lw=.9, label=r'Ground truth $\partial\mathcal{G}_i$'),
               Line2D([], [], color=teal, lw=1, ls='--', label=r'Fitted zero level $\phi_i=0$'),
               Line2D([], [], color=blue, lw=1.2, label=r'Certified surface $\partial\mathcal{B}_i=\mathcal{S}_i$')]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(.372, .13),
               fontsize=6.5, frameon=False, handlelength=2, labelspacing=.5)

    out = ROOT / 'paper' / 'figs'
    out.mkdir(exist_ok=True)
    # Replace atomically so an open IDE image preview does not block writing.
    rendered = metadata_dir / 'patch_pipeline_render.png'
    fig.savefig(rendered, dpi=350)
    os.replace(rendered, out / 'patch_pipeline.png')
    plt.close(fig)
    meta = {'shape_seed': 3, 'shape_index': 3, 'sdf_cell_m': field.h,
            'training_spacing_m': .005, 'certified_level_m': level,
            'native_cover_count': len(centers), 'cover_side_m': side,
            'render_grid': [501, 501], 'subdivision': False,
            'contour_spacing_m': level_spacing, 'contour_levels_m': levels.tolist(),
            'domain_lo_m': lo.tolist(), 'domain_hi_m': hi.tolist(),
            'zoom_center_m': center.tolist(),
            'zoom_width_m': 2*radius, 'ground_truth': 'Exact signed point-to-segment polygon distance',
            'note': 'Display grid only renders the field; cover selection uses Bernstein coefficients. '
                    'Panel d is the sufficient vertex certificate, with epigraph constraints understood.'}
    (metadata_dir / 'metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    main()
