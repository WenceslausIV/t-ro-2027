"""
Side-by-side GIFs: tabulated configuration-space barrier (left) vs surface-cover barriers (right).

    python prototype_3d/cover_vs_tab_gif.py dock    # docking into the tapered receptacle (Sec. VIII-C)
    python prototype_3d/cover_vs_tab_gif.py swap    # five-robot swap of Fig. 5 (needs five_cover.py swap)

Output: ../results/gif/{dock,swap}_tab_vs_cover.gif and a final-frame PNG.
"""
import itertools
import json
import os
import sys

import numpy as np

import franka3d as F
from dock_cover import body_field, cover_2d, make_cover_fn, tight_level_2d

ROOT = os.path.dirname(F.HERE)
sys.path.insert(0, ROOT)
import cspace_sdf_cbf_compare as C                                      # noqa: E402
import cspace_experiments as E                                          # noqa: E402

import matplotlib                                                       # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                         # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter            # noqa: E402
from matplotlib.patches import Polygon                                  # noqa: E402

OUT = os.path.join(ROOT, 'results', 'gif')
DT = .01
TAB, COV = '#1f77b4', '#d62728'
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})


def world(P, x):
    return P @ C.rot(x[2]).T + x[:2]


def tab_active(traj, fields):
    """Number of active tabulated barriers (pairs inside their activation region) per step."""
    n = np.zeros(len(traj), int)
    for (i, j), (f, D, l) in fields.items():
        t = np.einsum('kab,kb->ka', np.stack([C.rot(-th) for th in traj[:, j, 2]]), traj[:, i, :2] - traj[:, j, :2])
        n += np.all(np.abs(t) < D - C.ACT_DELTA, axis=1)
    return n


# ---------------------------------------------------------------------------------------------
def gif_dock(side=.005, steps=900, every=8, hold=25):
    from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH
    GA, GB, A, B, fields, _ = E.dock2_setup()
    f_tab, D, l_tab = fields['bspline']
    kw = dict(steps=steps, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
    Lt = E.simulate_team([GA, GB], [A, B], {(0, 1): fields['bspline']}, E.DOCK2_START, E.DOCK2_GOAL, 'cspace', **kw)

    fS, fC = body_field(GA, .0125, .15), body_field(GB, .01, .06)
    MS, MC = F.hessian_majorant(fS), F.hessian_majorant(fC)
    lS, lC = tight_level_2d(fS, MS, GA), tight_level_2d(fC, MC, GB)
    PC, sz = cover_2d(fC, lC, side)
    fn, counts = make_cover_fn(fS, MS, lS, PC, sz), []

    def counted(xi, xj):
        R, hv = fn(xi, xj)
        counts.append(len(hv))
        return R, hv
    Lc = E.simulate_team([GA, GB], [A, B], {(0, 1): counted}, E.DOCK2_START, E.DOCK2_GOAL, 'cover', **kw)

    runs = []
    for L, nb, col, name in ((Lt, None, TAB, f'tabulated C-space field (1 barrier), $l^\\star$ = {1e3 * l_tab:.1f} mm'),
                             (Lc, np.array(counts), COV, f'surface cover ({len(PC)} boxes of {1e3 * sz:.0f} mm), '
                              f'levels {1e3 * lS:.1f}/{1e3 * lC:.1f} mm')):
        traj = np.array(L['traj'])[:, 1]
        tip = traj[:, :2] + NOSE_TIP_X * np.stack([np.cos(traj[:, 2]), np.sin(traj[:, 2])], axis=1)
        gap = np.array([C.true_distance(GA, world(GB, x)) for x in traj])
        if nb is None:
            nb = tab_active(np.array(L['traj']), {(0, 1): fields['bspline']})
        runs.append(dict(traj=traj, depth=MOUTH_X - tip[:, 0], gap=gap, nb=nb, col=col, name=name,
                         t_med=1e3 * np.median(L['t_ctrl'])))
        print(f"{name}: depth {runs[-1]['depth'][-1]:.3f} gap {1e3 * gap[-1]:.1f} mm min gap {1e3 * gap.min():.1f} mm "
              f"t_med {runs[-1]['t_med']:.2f} ms bars {nb.max()}", flush=True)

    fig = plt.figure(figsize=(12, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.0, 1.0], hspace=.25, wspace=.12)
    arts = []
    for c, r in enumerate(runs):
        ax = fig.add_subplot(gs[0, c])
        az = ax.inset_axes([.50, .02, .49, .45])
        items = {}
        for a, lw in ((ax, 1.0), (az, .8)):
            a.add_patch(Polygon(GA, closed=True, color='.6', lw=0))
            if c == 0:
                a.add_patch(Polygon(A.dense, closed=True, fill=False, ec='k', lw=lw))
            items.setdefault('gt', []).append(a.add_patch(Polygon(GB, closed=True, color='.6', alpha=.8, lw=0)))
            if c == 0:
                items.setdefault('fit', []).append(a.add_patch(Polygon(B.dense, closed=True, fill=False, ec=r['col'], lw=lw)))
            else:
                items.setdefault('box', []).append(a.plot([], [], '.', ms=1.2 if a is az else .6, color=r['col'])[0])
            a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
        ax.set_xlim(-1.0, 4.3); ax.set_ylim(-1.1, 1.25)
        az.set_xlim(MOUTH_X - SEAT_DEPTH - .08, MOUTH_X + .12); az.set_ylim(-.22, .22)
        for s in az.spines.values():
            s.set_color(r['col'])
        ax.set_title(r['name'], fontsize=9, color=r['col'])
        items['tr'] = ax.plot([], [], color=r['col'], lw=.8)[0]
        items['txt'] = ax.text(.02, .97, '', transform=ax.transAxes, va='top', fontsize=9, family='monospace')
        arts.append(items)

    t = np.arange(steps) * DT
    axd = fig.add_subplot(gs[1, 0])
    axn = fig.add_subplot(gs[1, 1])
    lines = []
    for r in runs:
        n = len(r['depth'])
        axd.plot(t[:n], r['depth'], color=r['col'], lw=.6, alpha=.25)
        axn.plot(t[:n], np.maximum(r['nb'][:n], .8), color=r['col'], lw=.6, alpha=.25)
        lines.append((axd.plot([], [], color=r['col'], lw=1.6)[0], axn.plot([], [], color=r['col'], lw=1.6)[0]))
    axd.axhline(SEAT_DEPTH, color='k', lw=.5, ls='--'); axd.axhline(0, color='k', lw=.5, ls=':')
    axd.text(t[-1], SEAT_DEPTH + .05, 'full seating', ha='right', fontsize=8)
    axd.set_xlim(0, t[-1]); axd.set_ylim(-2.5, SEAT_DEPTH + .3)
    axd.set_xlabel('$t$ [s]'); axd.set_ylabel('nose depth [m]')
    axn.set_yscale('log'); axn.set_xlim(0, t[-1]); axn.set_ylim(.8, 3e3)
    axn.set_xlabel('$t$ [s]'); axn.set_ylabel('active barriers')
    for a in (axd, axn):
        a.spines[['top', 'right']].set_visible(False)
    n_frames = steps // every + hold

    def update(fr):
        for r, it, (ld, ln) in zip(runs, arts, lines):
            j = min(fr * every, len(r['traj']) - 1)
            x = r['traj'][j]
            for p in it['gt']:
                p.set_xy(world(GB, x))
            for p in it.get('fit', []):
                p.set_xy(B.world(x))
            for p in it.get('box', []):
                P = world(PC, x)
                p.set_data(P[:, 0], P[:, 1])
            it['tr'].set_data(r['traj'][:j + 1, 0], r['traj'][:j + 1, 1])
            it['txt'].set_text(f"nose depth {r['depth'][j]:+.3f} / {SEAT_DEPTH:.3f} m\n"
                               f"true gap   {1e3 * r['gap'][j]:6.1f} mm\n"
                               f"barriers   {r['nb'][j]:6d}\n"
                               f"time/step  {r['t_med']:6.2f} ms (median)")
            ld.set_data(t[:j + 1], r['depth'][:j + 1])
            ln.set_data(t[:j + 1], np.maximum(r['nb'][:j + 1], .8))
        return []

    save(fig, update, n_frames, 'dock_tab_vs_cover')


# ---------------------------------------------------------------------------------------------
def gif_swap(side=.005, every=10, hold=25):
    from cspace_cbf_5robots import five_robot_setup
    cfg = json.load(open(os.path.join(ROOT, 'results', 'five_jitter_si.json')))
    GT, shapes, fields, starts, goals, _ = five_robot_setup(
        scale=cfg['scale'], radius=cfg['radius'], jitter_seed=cfg['chosen']['jitter_seed'],
        jitter=tuple(cfg['jitter']), goal_shift=cfg.get('goal_shift'), n_robots=cfg.get('n_robots', 5))
    res = json.load(open(os.path.join(ROOT, 'results', f'five_cover_swap_{1e3 * side:g}.json')))
    N = len(GT)
    cols = plt.cm.tab10(np.arange(N))
    runs = []
    for m, name in (('cspace', 'tabulated C-space fields (1 barrier per pair)'),
                    ('cover', f"surface covers ({res['cover_setup']['side_mm']:g}-mm boxes)")):
        d = np.load(os.path.join(ROOT, 'results', f'five_cover_swap_{m}_{1e3 * side:g}.npz'))
        traj = d['traj']
        nb = tab_active(traj, fields) if m == 'cspace' else d['nb']
        PCs = [d[f'PC_{i}'] for i in range(N)] if m == 'cover' else None
        runs.append(dict(traj=traj, d=d['d'], nb=nb, PC=PCs, name=name, r=res[m],
                         col=TAB if m == 'cspace' else COV))
        print(m, res[m], flush=True)
    steps = max(len(r['traj']) for r in runs)
    lim = cfg['radius'] + 1.5

    fig = plt.figure(figsize=(11, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.6, 1.0], hspace=.28, wspace=.1)
    arts = []
    for c, r in enumerate(runs):
        ax = fig.add_subplot(gs[0, c])
        it = dict(gt=[], out=[], tr=[])
        for i in range(N):
            ax.plot(*goals[i, :2], 'x', color=cols[i], ms=6, mew=1.5)
            it['gt'].append(ax.add_patch(Polygon(GT[i], closed=True, color=cols[i], alpha=.35, lw=0)))
            if r['PC'] is None:
                it['out'].append(ax.add_patch(Polygon(shapes[i].dense, closed=True, fill=False, ec=cols[i], lw=1.0)))
            else:
                it['out'].append(ax.plot([], [], '.', ms=.8, color=cols[i])[0])
            it['tr'].append(ax.plot([], [], color=cols[i], lw=.8)[0])
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(r['name'], fontsize=10, color=r['col'])
        it['txt'] = ax.text(.02, .98, '', transform=ax.transAxes, va='top', fontsize=8.5, family='monospace')
        arts.append(it)
    t = np.arange(steps) * DT
    axd, axn = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    lines = []
    for r in runs:
        n = len(r['d'])
        axd.plot(t[:n], 1e3 * np.minimum(r['d'], .3), color=r['col'], lw=.6, alpha=.25)
        axn.plot(t[:n], np.maximum(r['nb'][:n], .8), color=r['col'], lw=.6, alpha=.25)
        lines.append((axd.plot([], [], color=r['col'], lw=1.6)[0], axn.plot([], [], color=r['col'], lw=1.6)[0]))
    axd.axhline(0, color='k', lw=.5)
    axd.set_xlim(0, t[-1]); axd.set_ylim(-5, 300)
    axd.set_xlabel('$t$ [s]'); axd.set_ylabel('min. true distance [mm]')
    axn.set_yscale('log'); axn.set_xlim(0, t[-1]); axn.set_ylim(.8, 3e4)
    axn.set_xlabel('$t$ [s]'); axn.set_ylabel('active barriers')
    for a in (axd, axn):
        a.spines[['top', 'right']].set_visible(False)
    n_frames = (steps + every - 1) // every + hold

    def update(fr):
        for r, it, (ld, ln) in zip(runs, arts, lines):
            j = min(fr * every, len(r['traj']) - 1)
            X = r['traj'][j]
            for i in range(N):
                it['gt'][i].set_xy(world(GT[i], X[i]))
                if r['PC'] is None:
                    it['out'][i].set_xy(shapes[i].world(X[i]))
                else:
                    P = world(r['PC'][i], X[i])
                    it['out'][i].set_data(P[:, 0], P[:, 1])
                it['tr'][i].set_data(r['traj'][:j + 1, i, 0], r['traj'][:j + 1, i, 1])
            it['txt'].set_text(f"t = {j * DT:5.2f} s\nmin true dist {1e3 * min(r['d'][:j + 1].min(), 9.99):6.1f} mm\n"
                               f"barriers {r['nb'][j]:6d}\ntime/step {r['r']['t_med_ms']:5.2f} ms (median)")
            ld.set_data(t[:j + 1], 1e3 * np.minimum(r['d'][:j + 1], .3))
            ln.set_data(t[:j + 1], np.maximum(r['nb'][:j + 1], .8))
        return []

    save(fig, update, n_frames, 'swap_tab_vs_cover')


def save(fig, update, n_frames, name):
    os.makedirs(OUT, exist_ok=True)
    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    anim.save(os.path.join(OUT, name + '.gif'), writer=PillowWriter(fps=20), dpi=80)
    update(n_frames - 1)
    fig.savefig(os.path.join(OUT, name + '_final.png'), dpi=110, bbox_inches='tight')
    print('>>>', os.path.join(OUT, name + '.gif'), flush=True)


if __name__ == '__main__':
    {'dock': gif_dock, 'swap': gif_swap}[sys.argv[1]]()
