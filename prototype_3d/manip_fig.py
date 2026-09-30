"""
Paper figure for the manipulator experiments (tro/figs/manip_paper.png), three columns:
(a) Franka among non-convex obstacles (franka3d, one trial), (b) inserting the hand into the star tube
(dolphin3d), (c) two Pandas in a shared workspace (dual3d, one trial). Top: poses (earlier poses faded,
closest approach opaque); bottom: certified lower bound of the mesh distance over time, with and without
the filter.

    python prototype_3d/manip_fig.py [franka_trial] [dual_trial]
"""
import json
import os
import sys

import numpy as np
import trimesh

import dolphin3d as Dp
import dual3d as D
import franka3d as F
from franka_gif import box_mesh, decimate, shade

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402

OUT = os.path.join(os.path.dirname(F.HERE), 'tro', 'figs', 'manip_paper.png')
GRAY, BLUE, PURPLE = (.72,) * 3, (.12, .47, .71), (.55, .45, .75)


def surface(fun, lo, hi, res=.004):
    from skimage import measure
    ax = [np.arange(lo[a] - 2 * res, hi[a] + 2 * res + 1e-9, res) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    V = fun(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)
    v, fc, _, _ = measure.marching_cubes(V, 0., spacing=(res,) * 3)
    return decimate(v + np.array([a[0] for a in ax]), fc, .006)


def arm_tris(meshes, links, T, rgb, alpha, base=None):
    tris, cols = [], []
    for L, (V, Fc) in zip(links, meshes):
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        tri = (V @ R.T + p)[Fc]
        tris.append(tri); cols.append(np.c_[shade(tri, rgb), np.full(len(tri), alpha)])
    if base is not None:
        tris.append(base); cols.append(np.c_[shade(base, rgb), np.full(len(base), alpha)])
    return tris, cols


def panel(ax, tris, cols, lim, view, title):
    ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols), edgecolor='none'))
    ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.set_zlim(*lim[2])
    ax.set_box_aspect([l[1] - l[0] for l in lim]); ax.view_init(*view); ax.set_axis_off()
    ax.set_title(title, y=.97)


def main():
    kf = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    links, obst, _ = F.build()
    meshes = [decimate(L['V'], L['F'], .01) for L in links]
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    Vb, Fb = decimate(np.asarray(base.vertices), np.asarray(base.faces), .012)
    B0 = Vb[Fb]
    plt.rcParams.update({'font.size': 8, 'font.family': 'serif'})
    fig = plt.figure(figsize=(7.16, 3.6))
    axes3 = [fig.add_axes([i / 3 - .02, .3, 1 / 3 + .04, .68], projection='3d') for i in range(3)]
    axes2 = [fig.add_axes([i / 3 + .055, .09, 1 / 3 - .08, .19]) for i in range(3)]

    # (a) obstacles
    tr = json.load(open(os.path.join(F.HERE, 'franka_results.json')))['trials'][kf]
    q0, qg = np.array(tr['q0']), np.array(tr['qg'])
    runs = {m: F.simulate(q0, qg, links, obst, m, record=True) for m in ('nominal', 'ours')}
    L = runs['ours']; kc = int(np.argmin(L['gap'])); n = len(L['q'])
    faces = np.array([f for Bx in F.OBSTACLES.values() for c, b in Bx for f in box_mesh(np.array(c), np.array(b))])
    ob = np.concatenate([faces[:, :3], faces[:, [0, 2, 3]]])
    tris, cols = [ob, B0], [np.c_[shade(ob, BLUE), np.ones(len(ob))], np.c_[shade(B0, GRAY), np.ones(len(B0))]]
    for s in sorted(set([0, n - 1] + [int(x) for x in np.linspace(0, n - 1, 5)])):
        t, c = arm_tris(meshes, links, F.fk(L['q'][s])[0], GRAY, .3); tris += t; cols += c
    t, c = arm_tris(meshes, links, F.fk(L['q'][kc])[0], GRAY, 1.); tris += t; cols += c
    panel(axes3[0], tris, cols, ((-.5, .8), (-.8, .8), (0, 1.05)), (24, -58), '(a) non-convex obstacles')
    plots = [('a', runs)]

    # (b) star tube
    hl, ho = Dp.build_hoop()
    q0 = Dp.start_config()
    runs_b = {m: Dp.simulate(q0, hl, ho, m) for m in ('nominal', 'ours')}
    Lb = runs_b['ours']; kb = len(Lb['q']) - 1
    lo, hi = Dp.HOOP_C - [Dp.HALF_T, Dp.R_OUT, Dp.R_OUT], Dp.HOOP_C + [Dp.HALF_T, Dp.R_OUT, Dp.R_OUT]
    Vt, Ft = surface(Dp.hoop_sdf, lo, hi)
    tube = Vt[Ft]
    tris, cols = [tube, B0], [np.c_[shade(tube, BLUE), np.full(len(tube), .55)], np.c_[shade(B0, GRAY), np.ones(len(B0))]]
    for s in np.linspace(0, kb, 4).astype(int)[:-1]:
        t, c = arm_tris(meshes, links, F.fk(Lb['q'][s])[0], GRAY, .3); tris += t; cols += c
    t, c = arm_tris(meshes, links, F.fk(Lb['q'][kb])[0], GRAY, 1.); tris += t; cols += c
    panel(axes3[1], tris, cols, ((-.2, 1.), (-.6, .6), (0, 1.0)), (18, -40), '(b) star tube insertion')
    plots.append(('b', runs_b))

    # (c) dual arm
    dl, _ = D.build()
    dres = json.load(open(os.path.join(F.HERE, 'dual_results.json')))['trials']
    kd = int(sys.argv[2]) if len(sys.argv) > 2 else min(range(len(dres)), key=lambda i: dres[i]['certified_min_gap_mm'])
    q0, qg = np.array(dres[kd]['q0']), np.array(dres[kd]['qg'])
    runs_c = {m: D.simulate(q0, qg, dl, m, record=True) for m in ('nominal', 'ours')}
    Lc = runs_c['ours']; kc = int(np.argmin(Lc['gap']))
    B2 = (Vb @ D.BASE2[:3, :3].T + D.BASE2[:3, 3])[Fb]
    tris, cols = [B0, B2], [np.c_[shade(B0, GRAY), np.ones(len(B0))], np.c_[shade(B2, PURPLE), np.ones(len(B2))]]
    for s, a in ((0, .25), (kc, 1.)):
        (T1, _, _), (T2, _, _) = D.fk2(Lc['q'][s])
        t, c = arm_tris(meshes, links, T1, GRAY, a); tris += t; cols += c
        t, c = arm_tris(meshes, links, T2, PURPLE, a); tris += t; cols += c
    panel(axes3[2], tris, cols, ((-.2, 1.1), (-.65, .65), (0, 1.05)), (20, -70), '(c) two arms')
    for ax, (tag, rr) in zip(axes2, plots + [('c', runs_c)]):
        for m, col, lab in (('nominal', '0.55', 'no filter'), ('ours', 'C0', 'ours')):
            g = 1e3 * np.array(rr[m]['gap'])
            if tag == 'c':                                   # sampled distance, 0 at contact, capped at 30 mm
                g = np.where(np.isfinite(g), g, 30.)
            g = np.minimum(np.where(np.isfinite(g), g, 50.), 50.)
            ax.plot(np.arange(len(g)) * F.DT, g, color=col, lw=1., label=lab)
        ax.axhline(0, color='k', lw=.5)
        ax.set_ylim(-45, 55); ax.set_xlabel('time [s]')
        if tag == 'a':
            ax.set_ylabel('distance [mm]'); ax.legend(loc='lower right', frameon=False, fontsize=7, ncol=2)
    fig.savefig(OUT, dpi=300)
    print('>>>', OUT, 'franka trial', kf, 'dual trial', kd)


if __name__ == '__main__':
    main()
