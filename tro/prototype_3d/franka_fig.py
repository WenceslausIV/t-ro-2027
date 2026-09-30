"""
Paper figure for the Franka experiment (paper/figs/franka_paper.png):
(a) poses of one filtered trial with the certified surfaces {phi = l} (links orange, obstacles blue),
(b) certified lower bound of the mesh-obstacle distance over time, with and without the filter.

    python prototype_3d/franka_fig.py [trial]
"""
import json
import os
import sys

import numpy as np

import franka3d as F
import trimesh
from franka_gif import box_mesh, decimate, level_mesh, shade

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402

OUT = os.path.join(os.path.dirname(F.HERE), 'paper', 'figs', 'franka_paper.png')


def main():
    res = json.load(open(os.path.join(F.HERE, 'franka_results.json')))
    trials = res['trials']
    if len(sys.argv) > 1:
        k = int(sys.argv[1])
    else:                                            # a trial that reaches its goal after touching h = 0
        k = min((t for t in trials if t['reached'] is not None), key=lambda t: t['min_h_mm'])['trial']
    tr = trials[k]
    import summed as SM
    import summed3d as S3
    SM.REFINE_DEPTH = 2  # Reproduce the historical subdivided experiment.
    links, obst, _ = F.build(cover_side=.006)
    S3.prep_all(links, obst.values())
    q0, qg = np.array(tr['q0']), np.array(tr['qg'])
    runs = {m: F.simulate(q0, qg, links, obst, m, record=True) for m in ('nominal', 'ours')}
    L = runs['ours']
    meshes = [decimate(Lk['V'], Lk['F'], .012) for Lk in links]
    lvl_links = [level_mesh(F.spline_from(Lk['f']), Lk['l']) for Lk in links]
    obs_faces = np.array([f for B in F.OBSTACLES.values() for c, b in B for f in box_mesh(np.array(c), np.array(b))])
    obs_tri = np.concatenate([obs_faces[:, :3], obs_faces[:, [0, 2, 3]]])
    lvl_obs = np.concatenate([(lambda VF: VF[0][VF[1]])(level_mesh(o['f'], o['l'])) for o in obst.values()])

    kc = int(np.argmin(L['gap']))                    # closest approach
    n = len(L['q'])
    idx = sorted(set([0, kc, n - 1] + [int(x) for x in np.linspace(0, n - 1, 5)]))

    plt.rcParams.update({'font.size': 8, 'font.family': 'serif'})
    fig = plt.figure(figsize=(3.5, 4.1))
    ax = fig.add_axes([-.08, .3, 1.16, .74], projection='3d')
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    Vb, Fb = decimate(np.asarray(base.vertices), np.asarray(base.faces), .012)
    tris = [obs_tri, lvl_obs, Vb[Fb]]
    cols = [np.c_[shade(obs_tri, (.12, .47, .71)), np.ones(len(obs_tri))],
            np.tile((.12, .47, .71, .15), (len(lvl_obs), 1)),
            np.c_[shade(Vb[Fb], (.72,) * 3), np.ones(len(Fb))]]
    for s in idx:
        T, _, _ = F.fk(L['q'][s])
        fade = 1. if s == kc else .35
        for Lk, (V, Fc), (Vl, Fl) in zip(links, meshes, lvl_links):
            R, p = T[Lk['frame']][:3, :3], T[Lk['frame']][:3, 3]
            tri = (V @ R.T + p)[Fc]
            tris.append(tri)
            cols.append(np.c_[shade(tri, (.72,) * 3), np.full(len(tri), fade)])
            if s == kc:
                tris.append((Vl @ R.T + p)[Fl])
                cols.append(np.tile((1., .5, .05, .25), (len(Fl), 1)))
    ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols), edgecolor='none'))
    ax.set_xlim(-.6, .8); ax.set_ylim(-.8, .8); ax.set_zlim(0, 1.1)
    ax.set_box_aspect((1.4, 1.6, 1.1)); ax.view_init(24, -58)
    ax.set_axis_off()
    fig.text(.02, .95, '(a)')

    ax2 = fig.add_axes([.17, .1, .79, .2])
    for m, c, lab in (('nominal', '0.5', 'no filter'), ('ours', 'C0', 'ours')):
        g = np.minimum(1e3 * np.array(runs[m]['gap']), 50.)     # real_gap reports > 5 cm as inf
        ax2.plot(np.arange(len(g)) * F.DT, g, color=c, lw=1.1, label=lab)
    ax2.axhline(0, color='k', lw=.6)
    ax2.set_xlabel('time [s]'); ax2.set_ylabel('distance [mm]')
    ax2.set_ylim(-45, 55); ax2.legend(loc='lower right', frameon=False, ncol=2)
    fig.text(.02, .3, '(b)')
    fig.savefig(OUT, dpi=300)
    print('>>>', OUT, 'trial', k, 'closest step', kc)


if __name__ == '__main__':
    main()
