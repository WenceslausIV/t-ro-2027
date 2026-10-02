"""Figures for results/star_random_goals (for presenting the comparison).

    python prototype_3d/star_random_goals_fig.py spheres OUT.png   # the three sphere decompositions at the start pose
    python prototype_3d/star_random_goals_fig.py trial K OUT.png   # one goal: our arm vs the best sphere setting
"""
import json
import os
import sys

import numpy as np

import star_random_goals as G
import continuous_baselines as CBL
import franka3d as F
from franka_gif import decimate, shade
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection         # noqa: E402

ROBOT, TUBE, RSPH, TSPH = (.66, .68, .71), (.12, .47, .71), (.55, .55, .58), (.20, .55, .85)
VIEW, LIM = (18, -35), dict(x=(-.40, .50), y=(-.40, .40), z=(.25, .90))
SETTINGS = [('64 / link, tube 256', 64, 256, .010), ('64 / link, tube 1024', 64, 1024, .005),
            ('128 / link, tube 1024', 128, 1024, .005)]


def sphere_mesh(C, r, n=10):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, np.pi, n // 2 + 1))
    P = np.stack([np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)], -1)          # (m, n, 3)
    quads = np.stack([P[:-1, :-1], P[:-1, 1:], P[1:, 1:], P[1:, :-1]], 2).reshape(-1, 4, 3)
    return (C[:, None, None] + r[:, None, None, None] * quads[None]).reshape(-1, 4, 3)


def setup_axes(ax, title):
    ax.set_xlim(*LIM['x']); ax.set_ylim(*LIM['y']); ax.set_zlim(*LIM['z'])
    ax.set_box_aspect((1., .9, .85)); ax.view_init(*VIEW); ax.set_axis_off()
    ax.set_title(title, fontsize=10, y=.98)


def arm_tris(q, links, meshes):
    T, _, _ = F.fk(q)
    return np.concatenate([(V @ T[L['frame']][:3, :3].T + T[L['frame']][:3, 3])[Fc] for L, (V, Fc) in zip(links, meshes)])


def fig_spheres(out):
    q0, links, B, RB, pB = G.scene()
    Vt, Ft = G.SV.exact_surface(B)
    tube = (Vt @ RB.T + pB)[Ft]
    meshes = [decimate(L['V'], L['F'], .01) for L in links]
    fig = plt.figure(figsize=(16, 5.2))
    ax = fig.add_subplot(1, 4, 1, projection='3d')
    arm = arm_tris(q0, links, meshes)
    ax.add_collection3d(Poly3DCollection(np.concatenate([tube, arm]), edgecolor='none', facecolors=np.concatenate(
        [np.c_[shade(tube, TUBE), np.full(len(tube), .8)], np.c_[shade(arm, ROBOT), np.ones(len(arm))]])))
    setup_axes(ax, 'ground truth (start pose)\nours certifies these surfaces')
    T, _, _ = F.fk(q0)
    for j, (name, kl, kt, vox) in enumerate(SETTINGS):
        prim = CBL.primitives(links, 'spheres_kmeans', .01, kl)
        Ct, rt = G.tube_spheres(B, vox, kt)
        Cw = Ct @ RB.T + pB
        X = np.concatenate([P @ T[L['frame']][:3, :3].T + T[L['frame']][:3, 3] for L, (P, _) in zip(links, prim)])
        rl = np.concatenate([r for _, r in prim])
        h0 = float((np.linalg.norm(X[:, None] - Cw[None], axis=2) - rl[:, None] - rt[None]).min())
        rs, ts = sphere_mesh(X, rl), sphere_mesh(Cw, rt)
        ax = fig.add_subplot(1, 4, j + 2, projection='3d')
        ax.add_collection3d(Poly3DCollection(ts, facecolor=TSPH + (.18,), edgecolor='none'))
        ax.add_collection3d(Poly3DCollection(rs, facecolor=RSPH + (.35,), edgecolor='none'))
        state = 'start INSIDE its safe set' if h0 >= 0 else f'start OUTSIDE its safe set (overlap {-1e3 * h0:.0f} mm)'
        setup_axes(ax, f'spheres {name}\n{len(rl)} robot + {len(rt)} tube spheres, '
                       f'radii {1e3 * rt.min():.0f}-{1e3 * rt.max():.0f} mm (tube)\n{state}')
    fig.subplots_adjust(left=0, right=1, bottom=0, top=.86, wspace=0)
    fig.savefig(out, dpi=170)
    print(out)


def fig_trial(k, out, ours_folder='ours', sph_folder='spheres_L64_T1024', ncol=4):
    q0, links, B, RB, pB = G.scene()
    Vt, Ft = G.SV.exact_surface(B)
    tube = (Vt @ RB.T + pB)[Ft]
    meshes = [decimate(L['V'], L['F'], .01) for L in links]
    goals = json.loads((G.OUT / 'goals.json').read_text())
    qg = np.array(goals['goals'][k])
    fig = plt.figure(figsize=(3.2 * ncol, 6.4))
    for row, (label, folder) in enumerate([('ours', ours_folder), ('spheres', sph_folder)]):
        z = np.load(G.OUT / folder / f'trial_{k:02d}.npz')
        m = json.loads((G.OUT / folder / f'trial_{k:02d}.json').read_text())['metrics']
        q = z['q']
        steps = [int(x) for x in np.linspace(0, len(q) - 1, ncol)]
        for c, s in enumerate(steps):
            ax = fig.add_subplot(2, ncol, row * ncol + c + 1, projection='3d')
            arm = arm_tris(q[s], links, meshes)
            goal = arm_tris(qg, links, meshes)
            ax.add_collection3d(Poly3DCollection(np.concatenate([tube, arm, goal]), edgecolor='none', facecolors=np.concatenate(
                [np.c_[shade(tube, TUBE), np.full(len(tube), .75)], np.c_[shade(arm, ROBOT), np.ones(len(arm))],
                 np.c_[shade(goal, (.2, .7, .3)), np.full(len(goal), .12)]])))
            res = 'reached' if m['reached_s'] is not None else 'stopped'
            setup_axes(ax, f"{label}: t = {s * F.DT:.1f} s" + (f"  ({res})" if c == ncol - 1 else ''))
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == '__main__':
    if sys.argv[1] == 'spheres':
        fig_spheres(sys.argv[2])
    else:
        fig_trial(int(sys.argv[2]), sys.argv[3], *sys.argv[4:])
