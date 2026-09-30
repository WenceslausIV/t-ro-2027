"""
GIF of one dual-arm trial (dual_results.json): left without the filter, right with it.
Arm 1 light gray, arm 2 darker gray; certified link surfaces of arm 1 translucent orange.

    python prototype_3d/dual_gif.py [trial]      -> prototype_3d/dual_trial<k>.gif
"""
import json
import os
import sys

import numpy as np

import dual3d as D
import trimesh
import franka3d as F
from franka_gif import decimate, level_mesh, shade

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter        # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402

EVERY, HOLD = 6, 15


def main():
    trials = json.load(open(os.path.join(F.HERE, 'dual_results.json')))['trials']
    k = int(sys.argv[1]) if len(sys.argv) > 1 else min(range(len(trials)),
                                                        key=lambda i: trials[i]['certified_min_gap_mm'])
    import summed as SM
    import summed3d as S3
    SM.REFINE_DEPTH = 2  # Reproduce the historical subdivided experiment.
    links, _ = D.build(cover_side=.006)
    S3.prep_all(links)
    q0, qg = np.array(trials[k]['q0']), np.array(trials[k]['qg'])
    runs = [D.simulate(q0, qg, links, m) for m in ('nominal', 'ours')]
    meshes = [decimate(L['V'], L['F']) for L in links]
    lvl = [level_mesh(L['fs'], L['l']) for L in links]
    n = max(len(r['q']) for r in runs)
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    Vb, Fb = decimate(np.asarray(base.vertices), np.asarray(base.faces), .012)
    bases = [Vb[Fb], (Vb @ D.BASE2[:3, :3].T + D.BASE2[:3, 3])[Fb]]
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})
    fig = plt.figure(figsize=(11, 5.5))
    axes = [fig.add_axes([.5 * i - .04, 0., .58, .9], projection='3d') for i in range(2)]
    titles = ['no safety filter', 'ours (64 link pairs)']

    def draw(fr):
        s = min(fr * EVERY, n - 1)
        for ax, r, ttl in zip(axes, runs, titles):
            ax.cla()
            q = r['q'][min(s, len(r['q']) - 1)]
            (T1, _, _), (T2, _, _) = D.fk2(q)
            gap = D.approx_gap(q, links)
            tris = list(bases)
            cols = [np.c_[shade(b, (.7,) * 3), np.ones(len(b))] for b in bases]
            for T, rgb in ((T1, (.8,) * 3), (T2, (.55,) * 3)):
                for L, (V, Fc), (Vl, Fl) in zip(links, meshes, lvl):
                    R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
                    tri = (V @ R.T + p)[Fc]
                    tris.append(tri)
                    cols.append(np.c_[shade(tri, (1, .1, .1) if gap < .002 else rgb), np.ones(len(tri))])
                    if T is T1:
                        tris.append((Vl @ R.T + p)[Fl])
                        cols.append(np.tile((1., .5, .05, .1), (len(Fl), 1)))
            ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols),
                                                 edgecolor='none'))
            ax.set_xlim(-.15, 1.05); ax.set_ylim(-.6, .6); ax.set_zlim(0, 1.)
            ax.set_box_aspect((1.2, 1.2, 1.)); ax.view_init(22, -70); ax.set_axis_off()
            gmin = min(r['gap'][:min(s, len(r['gap']) - 1) + 1])
            txt = 'CONTACT' if gmin < .002 else ('> 30 mm' if not np.isfinite(gmin) else f'{1e3 * gmin:4.1f} mm')
            ax.set_title(f"{ttl}\nt = {s * F.DT:4.2f} s   min. distance so far {txt}")
        return []

    out = os.path.join(F.HERE, f'dual_trial{k}.gif')
    FuncAnimation(fig, draw, frames=n // EVERY + HOLD, blit=False).save(out, writer=PillowWriter(fps=15), dpi=85)
    print('>>>', out)


if __name__ == '__main__':
    main()
