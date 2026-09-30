"""
GIF of one Franka trial: left without the safety filter, right with it (franka3d.py barriers).
Real shapes filled (links light gray, obstacles dark gray; a link intersecting an obstacle turns red),
certified SDF level sets {phi = l} translucent (links orange, obstacles blue).

    python prototype_3d/franka_gif.py [trial] [--png]   -> prototype_3d/franka_trial<k>.gif (.png preview)
"""
import json
import os
import sys

import numpy as np

import franka3d as F
from skimage import measure

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter        # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402

EVERY, HOLD = 6, 15
LIGHT = np.array([.4, -.5, .77])


def decimate(V, Fc, cell=.015):
    """Vertex clustering: merge vertices on a grid, drop degenerate triangles (for drawing only)."""
    key = np.floor(V / cell).astype(int)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    Vn = np.zeros((inv.max() + 1, 3)); np.add.at(Vn, inv, V)
    Vn /= np.bincount(inv)[:, None]
    Fn = inv[Fc]
    Fn = Fn[(Fn[:, 0] != Fn[:, 1]) & (Fn[:, 1] != Fn[:, 2]) & (Fn[:, 0] != Fn[:, 2])]
    _, keep = np.unique(np.sort(Fn, axis=1), axis=0, return_index=True)     # drop duplicates, keep winding
    return Vn, Fn[np.sort(keep)]


def shade(tris, rgb):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    k = .45 + .55 * np.abs(n @ LIGHT)
    return np.clip(k[:, None] * np.array(rgb), 0, 1)


def level_mesh(f, l, res=.005):
    """Triangles of the level set {phi = l} of a spline SDF (marching cubes, drawing only)."""
    ax = [np.arange(f.lo[a], f.lo[a] + f.K[a] * f.h, res) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    V = f.eval(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)
    v, fc, _, _ = measure.marching_cubes(V, l, spacing=(res,) * 3)
    return decimate(v + f.lo, fc, .012)


def box_mesh(c, b):
    s = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float) * b + c
    quads = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]
    return [s[qd] for qd in quads]


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 0
    import summed as SM
    import summed3d as S3
    SM.REFINE_DEPTH = 2  # Reproduce the historical subdivided experiment.
    links, obst, _ = F.build(cover_side=.006)
    S3.prep_all(links, obst.values())
    tr = json.load(open(os.path.join(F.HERE, 'franka_results.json')))['trials'][k]
    q0, qg = np.array(tr['q0']), np.array(tr['qg'])
    runs = [F.simulate(q0, qg, links, obst, m, record=True) for m in ('nominal', 'ours')]
    meshes = [decimate(L['V'], L['F']) for L in links]
    obs_faces = [f for B in F.OBSTACLES.values() for c, b in B for f in box_mesh(np.array(c), np.array(b))]
    OBS_RGB = shade(np.array([f[:3] for f in obs_faces]), (.35,) * 3)
    lvl_links = [level_mesh(F.spline_from(L['f']), L['l']) for L in links]
    lvl_obs = []
    for o in obst.values():
        V, Fc = level_mesh(o['f'], o['l'])
        lvl_obs.append(V[Fc])
    lvl_obs = np.concatenate(lvl_obs)
    n = max(len(r['q']) for r in runs)

    plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})
    fig = plt.figure(figsize=(11, 5.5))
    axes = [fig.add_axes([.5 * i - .04, 0., .58, .9], projection='3d') for i in range(2)]
    titles = ['no safety filter', 'ours (all link-obstacle pairs)']

    def draw(fr):
        s = min(fr * EVERY, n - 1)
        for ax, r, ttl in zip(axes, runs, titles):
            ax.cla()
            q = r['q'][min(s, len(r['q']) - 1)]
            T, _, _ = F.fk(q)
            gaps = F.real_gap(q, links, obst, per_link=True)
            tris = [np.array(obs_faces)[:, :3], np.array(obs_faces)[:, [0, 2, 3]], lvl_obs]
            cols = [np.c_[OBS_RGB, np.ones(len(OBS_RGB))]] * 2 + [np.tile((.12, .47, .71, .18), (len(lvl_obs), 1))]
            for L, (V, Fc), g, (Vl, Fl) in zip(links, meshes, gaps, lvl_links):
                R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
                tri = (V @ R.T + p)[Fc]
                tris.append(tri)
                cols.append(np.c_[shade(tri, (1, .1, .1) if g < 0 else (.75,) * 3), np.ones(len(tri))])
                tris.append((Vl @ R.T + p)[Fl])
                cols.append(np.tile((1., .5, .05, .1), (len(Fl), 1)))
            ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols),
                                                 edgecolor='none'))            # one collection: global depth sort
            ax.set_xlim(-.45, .75); ax.set_ylim(-.8, .8); ax.set_zlim(0, 1.)
            ax.set_box_aspect((1.2, 1.6, 1.)); ax.view_init(24, -58); ax.set_axis_off()
            gmin = min(r['gap'][:min(s, len(r['gap']) - 1) + 1])
            ax.set_title(f"{ttl}\nt = {s * F.DT:4.2f} s   min. real gap so far "
                         + ("COLLISION" if gmin < 0 else f"{1e3 * gmin:5.1f} mm"))
        return []

    if '--png' in sys.argv:                          # preview: the closest approach of the filtered run
        draw(int(np.argmin(runs[1]['gap'])) // EVERY)
        out = os.path.join(F.HERE, f'franka_trial{k}.png')
        fig.savefig(out, dpi=110, bbox_inches='tight')
        print('>>>', out)
        return
    out = os.path.join(F.HERE, f'franka_trial{k}.gif')
    FuncAnimation(fig, draw, frames=n // EVERY + HOLD, blit=False).save(out, writer=PillowWriter(fps=15), dpi=85)
    print('>>>', out)


if __name__ == '__main__':
    main()
