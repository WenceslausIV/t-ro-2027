"""
Paper figure for the dual-arm experiment (paper/figs/dual_paper.png): two Pandas move to their goals and cross;
top row without the safety filter (red stars mark sampled near-contact), bottom row with it; columns are time instants.

    python prototype_3d/dual_fig.py [trial] [n_columns]
"""
import json
import os
import sys

import numpy as np
import trimesh

import dual3d as D
import franka3d as F
from franka_gif import decimate, shade

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection             # noqa: E402
from mpl_toolkits.mplot3d import proj3d
from matplotlib.lines import Line2D

OUT = os.path.join(os.path.dirname(F.HERE), 'paper', 'figs', 'dual_paper.png')
ARM1, ARM2 = (.66, .68, .71), (.12, .47, .71)


def contact_points(q, links):
    """Nearest sampled surface-pair midpoints for link pairs within 2 mm."""
    (T1, _, _), (T2, _, _) = D.fk2(q)
    points = []
    for ia, A in enumerate(links):
        for ib, B in enumerate(links):
            pA, pB = T1[A['frame']][:3, 3], T2[B['frame']][:3, 3]
            if np.linalg.norm(pA - pB) > A['rad'] + B['rad']:
                continue
            Rab = T2[B['frame']][:3, :3].T @ T1[A['frame']][:3, :3]
            tab = T2[B['frame']][:3, :3].T @ (pA - pB)
            P = A['smp'] @ Rab.T + tab
            P = P[np.all((P > B['bb'][0]) & (P < B['bb'][1]), axis=1)]
            if len(P):
                distances, indices = B['tree'].query(P, distance_upper_bound=.002)
                k = int(np.argmin(distances))
                if distances[k] < .002:
                    midpoint = .5 * (P[k] + B['smp'][indices[k]])
                    world = T2[B['frame']][:3, :3] @ midpoint + pB
                    if all(np.linalg.norm(world - other) > .05 for other in points):
                        points.append(world)
    return points


def pick_trial(links, trials):
    """Among trials in which both arms reach their goals with the filter, the one whose unfiltered run stays in
    contact the longest."""
    best, score = None, -1
    for k, t in enumerate(trials):
        if t['reached'] is None:
            continue
        nom = D.simulate(np.array(t['q0']), np.array(t['qg']), links, 'nominal')
        c = int(np.sum(np.array(nom['gap']) < .002))
        if c > score:
            best, score = k, c
    return best


def main():
    import summed3d as S3
    import summed as SM
    SM.REFINE_DEPTH = 2  # Reproduce the historical subdivided experiment.
    links, _ = D.build(cover_side=.006)
    S3.prep_all(links)                                  # summed-field barriers ('ours' in dual3d.simulate)
    trials = json.load(open(os.path.join(F.HERE, 'dual_results.json')))['trials']
    k = int(sys.argv[1]) if len(sys.argv) > 1 else pick_trial(links, trials)
    ncol = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    q0, qg = np.array(trials[k]['q0']), np.array(trials[k]['qg'])
    runs = {m: D.simulate(q0, qg, links, m) for m in ('nominal', 'ours')}
    T_end = max(len(r['q']) for r in runs.values()) - 1
    nom_gap = np.array(runs['nominal']['gap'])
    k_first = int(np.argmax(nom_gap < .002)) if (nom_gap < .002).any() else T_end // 2
    steps = [int(x) for x in np.linspace(0, T_end, ncol)]
    steps[int(np.argmin([abs(s - k_first) for s in steps[1:-1]])) + 1] = k_first   # show first contact
    meshes = [decimate(L['V'], L['F'], .01) for L in links]
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    Vb, Fb = decimate(np.asarray(base.vertices), np.asarray(base.faces), .012)
    bases = [Vb[Fb], (Vb @ D.BASE2[:3, :3].T + D.BASE2[:3, 3])[Fb]]

    plt.rcParams.update({'font.size': 8, 'font.family': 'serif'})
    Y0, AH, B1 = 1.12, 1.5, -.25           # row spacing, axes height, bottom of the second row [in]
    FH = B1 + Y0 + AH
    fig = plt.figure(figsize=(7.16, FH))
    w = (1. - .035) / ncol
    for row, m in enumerate(('nominal', 'ours')):
        r = runs[m]
        for c, s in enumerate(steps):
            ax = fig.add_axes([.035 + c * w - .02, (B1 + Y0 * (1 - row)) / FH, w + .04, AH / FH], projection='3d')
            q = r['q'][min(s, len(r['q']) - 1)]
            (T1, _, _), (T2, _, _) = D.fk2(q)
            contacts = contact_points(q, links) if m == 'nominal' else []
            tris = list(bases)
            cols = [np.c_[shade(b, rgb), np.ones(len(b))] for b, rgb in zip(bases, (ARM1, ARM2))]
            for T, rgb in ((T1, ARM1), (T2, ARM2)):
                for i, (L, (V, Fc)) in enumerate(zip(links, meshes)):
                    R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
                    tri = (V @ R.T + p)[Fc]
                    tris.append(tri)
                    cols.append(np.c_[shade(tri, rgb), np.ones(len(tri))])
            coll = Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols), edgecolor='none')
            coll.set_clip_on(False)                         # arms reaching past the axes box stay visible
            ax.add_collection3d(coll)
            ax.set_xlim(-.25, 1.0); ax.set_ylim(-.45, .45); ax.set_zlim(0, .95)
            ax.set_box_aspect((1.25, .9, .95)); ax.view_init(16, -78); ax.set_axis_off()
            for point in contacts:
                x, y, _ = proj3d.proj_transform(*point, ax.get_proj())
                ax.add_artist(Line2D([x], [y], marker='*', color='crimson', markersize=9,
                                     linestyle='none', transform=ax.transData, zorder=100,
                                     clip_on=False))
            ax.patch.set_alpha(0.)                          # neighbouring axes overlap: do not hide their arms
            if row == 0:
                ax.set_title(f't = {s * F.DT:.1f} s', y=.93, fontsize=8)
        fig.text(.012, (B1 + Y0 * (1 - row) + AH / 2) / FH, 'no filter' if m == 'nominal' else 'ours', rotation=90, va='center',
                 ha='center', fontsize=8)
    render_dir = os.path.join(os.path.dirname(F.HERE), 'results', 'figure_style')
    os.makedirs(render_dir, exist_ok=True)
    temp = os.path.join(render_dir, 'dual_paper.png')
    fig.savefig(temp, dpi=300)
    os.replace(temp, OUT)
    print('>>>', OUT, 'trial', k, 'first contact step', k_first, 'steps', steps,
          'reached (ours)', runs['ours']['reached'])


if __name__ == '__main__':
    main()
