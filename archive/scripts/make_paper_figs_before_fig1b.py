"""
Publication figures for tro/main.tex (replace the GIF screenshots).

    python make_paper_figs.py        -> tro/figs/dock_paper.png, tro/figs/five_paper.png

Uses cached certified fields and the shared five_robot_setup for the five random shapes.
The statistical experiments use the original cache/shapes5.* cache.
"""
import os
import sys

import numpy as np

import cspace_sdf_cbf_compare as C
from cspace_experiments import FIGS, TWO_PI, ab_field, simulate_team

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
from matplotlib.patches import Circle, Polygon                           # noqa: E402

plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
                     'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
BLUE, ORANGE = '#1f77b4', '#ff7f0e'


def bezier_ab_field(A, B):
    """Piecewise Bezier field of the same pair (same training data as the B-spline field)."""
    D = round(A.rho + B.rho + 0.3, 2)
    th = np.arange(96) * TWO_PI / 96
    g, S = C.cspace_sdf_slices(A, B, th, D)
    gxy = np.linspace(-D, D, 100)
    f = C.Field3D('bezier', D, 16, 16)
    f.fit(C.sample_slices(g, S, gxy, gxy).transpose(1, 2, 0), gxy, gxy, th)
    l_star, ok, _, _ = C.Certifier(f, A, B).certify(verbose=False)
    assert ok
    return f, D, l_star


def draw_pose(ax, shape, gt, x, color, alpha=1.0, fill=True, zorder=None):
    kw = {} if zorder is None else dict(zorder=zorder)
    if fill:
        ax.add_patch(Polygon(gt @ C.rot(x[2]).T + x[:2], closed=True, color='0.6', alpha=alpha, lw=0, **kw))
    ax.add_patch(Polygon(shape.world(x), closed=True, fill=False, ec=color, lw=1.0, alpha=alpha, **kw))


def fading_path(ax, xy, color, lw=2.4, a0=0., a1=1., zorder=2, ds=.01):
    """Path whose opacity grows in proportion to the arclength travelled: fully transparent at the start,
    fully opaque at the end. The path is resampled every ds and drawn as non-overlapping butt-capped
    segments, so the transparency does not accumulate where segments would overlap."""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import to_rgba
    s = np.r_[0, np.linalg.norm(np.diff(xy, axis=0), axis=1).cumsum()]
    if s[-1] < 1e-9:
        return
    q = np.arange(0, s[-1] + ds, ds).clip(max=s[-1])
    p = np.column_stack([np.interp(q, s, xy[:, k]) for k in range(2)])
    seg = np.stack([p[:-1], p[1:]], axis=1)
    rgba = np.tile(to_rgba(color), (len(seg), 1))
    rgba[:, 3] = a0 + (a1 - a0) * (q[:-1] + q[1:]) / (2 * s[-1])
    ax.add_collection(LineCollection(seg, colors=rgba, linewidths=lw, capstyle='butt', joinstyle='miter',
                                     zorder=zorder))


def fig_dock():
    A, B, fb, D, lb = ab_field()
    fz, Dz, lz = bezier_ab_field(A, B)
    GT2 = [A.dense, B.dense]
    starts = np.array([[0.0, 0.0, 0.0], [2.6, 0.25, -0.6]])
    goals = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, np.pi / 3]])
    runs = {}
    for name, (f, DD, l) in (('B-spline ($C^2$)', (fb, D, lb)), ('piecewise Bézier ($C^0$)', (fz, Dz, lz))):
        runs[name] = simulate_team(GT2, [A, B], {(0, 1): (f, DD, l)}, starts, goals, 'cspace', steps=700,
                                   v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        L = runs[name]
        print(f"dock {name}: final x {L['final'][1, 0]:.4f} min h {np.nanmin(L['h_t']):+.4f} tv {L['tv']:.2f}")

    fig = plt.figure(figsize=(7.16, 3.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.8], wspace=0.2, hspace=0.5)
    t = np.arange(700) * 0.01
    for k, (name, L) in enumerate(runs.items()):
        ax = fig.add_subplot(gs[0, k])
        traj = np.array(L['traj'])
        for j in range(0, len(traj), 70):
            draw_pose(ax, B, B.dense, traj[j, 1], ORANGE if k else BLUE, alpha=0.25, fill=False)
        draw_pose(ax, A, A.dense, traj[-1, 0], 'k')
        draw_pose(ax, B, B.dense, traj[-1, 1], ORANGE if k else BLUE)
        ax.plot(traj[:, 1, 0], traj[:, 1, 1], color=ORANGE if k else BLUE, lw=0.7)
        ax.set_xlim(-0.7, 3.2)
        ax.set_ylim(-0.9, 0.9)
        ax.set_aspect('equal')
        ax.set_title(f"({'ab'[k]}) {name}")
        ax.set_xlabel('$x$ [m]')
        if k == 0:
            ax.set_ylabel('$y$ [m]')
    axh = fig.add_subplot(gs[1, 0])
    axw = fig.add_subplot(gs[1, 1])
    for k, (name, L) in enumerate(runs.items()):
        col = ORANGE if k else BLUE
        n = len(L['h_t'])
        axh.plot(t[:n], 1e3 * np.array(L['h_t']), color=col, lw=0.9, label=name)
        axw.plot(t[:n], np.array(L['u_t'])[:, 5], color=col, lw=0.6, label=name)
    axh.axhline(0, color='k', lw=0.5)
    axh.set_ylim(-10, 60)
    axh.set_xlabel('$t$ [s]')
    axh.set_ylabel('$h$ [mm]')
    axh.set_title('(c) barrier value')
    axh.legend(loc='upper right', frameon=False, fontsize=6)
    axh.set_xlim(0, 7)
    axw.set_xlim(0, 7)
    axw.set_xlabel('$t$ [s]')
    axw.set_ylabel(r'$\omega_B$ [rad/s]')
    axw.set_title('(d) angular input of B')
    fig.savefig(os.path.join(FIGS, 'dock_paper.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('>>> tro/figs/dock_paper.png')


def fig_dock2():
    """Docking into a tapered receptacle: ours vs baselines."""
    from cspace_experiments import DOCK2_GOAL, DOCK2_START, dock2_setup
    from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH
    GA, GB, A, B, fields, meta = dock2_setup()
    # surface cover (prototype_3d/dock_cover.py): fields of the ground-truth polygons, tightest levels,
    # 5-mm boxes on the craft's certified curve
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'prototype_3d'))
    import franka3d as F3
    from dock_cover import body_field, tight_level_2d, cover_2d, make_summed_fn
    import summed as SM
    fS, fC = body_field(GA, .0125, .15), body_field(GB, .01, .06)
    MS, MC = F3.hessian_majorant(fS), F3.hessian_majorant(fC)
    lS, lC = tight_level_2d(fS, MS, GA), tight_level_2d(fC, MC, GB)
    PC, sz = cover_2d(fC, lC, .005)
    # surface cover with summed-field barriers and Bernstein coefficient constraints (prototype_3d/summed.py)
    fields = dict(fields, cover=make_summed_fn(fS, MS, lS, SM.prepare_body(fC, lC, PC, sz, 2)))

    def level_curve(f, l, res=.002):
        xs_ = np.arange(f.lo[0], f.lo[0] + f.K[0] * f.h, res)
        ys_ = np.arange(f.lo[1], f.lo[1] + f.K[1] * f.h, res)
        XX, YY = np.meshgrid(xs_, ys_, indexing='ij')
        V = f.eval(np.c_[XX.ravel(), YY.ravel(), np.zeros(XX.size)]).reshape(XX.shape)
        cs = plt.figure().add_subplot().contour(xs_, ys_, V.T, levels=[l])
        seg = max(cs.allsegs[0], key=len)
        plt.close(cs.axes.figure)
        return seg

    class Level:                                   # certified level curve, moved like a fitted boundary
        def __init__(self, P):
            self.dense = P

        def world(self, x):
            return self.dense @ C.rot(x[2]).T + x[:2]

    certS, certC = Level(level_curve(fS, lS)), Level(level_curve(fC, lC, .001))
    runs = {}
    # categorical palette in fixed order (blue, orange, aqua, red); line styles as secondary encoding
    STYLE = {'ours, surface cover': dict(lw=1.4, ls='-'), 'ours, compiled field': dict(lw=1.1, ls='--'),
             'closest point': dict(lw=1.1, ls='-.'), 'circle': dict(lw=1.2, ls=':')}
    for name, method, kind, col in (('ours, surface cover', 'summed', 'cover', '#2a78d6'),
                                    ('ours, compiled field', 'cspace', 'bspline', '#eb6834'),
                                    ('closest point', 'closest', 'bspline', '#1baf7a'),
                                    ('circle', 'circle', 'bspline', '#e34948')):
        L = simulate_team([GA, GB], [A, B], {(0, 1): fields[kind]}, DOCK2_START, DOCK2_GOAL, method,
                          steps=900, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
        traj = np.array(L['traj'])
        tip = traj[:, 1, :2] + np.stack([np.cos(traj[:, 1, 2]), np.sin(traj[:, 1, 2])], axis=1) * NOSE_TIP_X
        runs[name] = (L, traj, MOUTH_X - tip[:, 0], col)
        print(f"dock2 {name}: depth {MOUTH_X - tip[-1, 0]:.3f} m")

    fig = plt.figure(figsize=(7.16, 3.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.85], wspace=0.2, hspace=0.5)
    for k, name in enumerate(('ours, surface cover', 'circle')):
        L, traj, _, col = runs[name]
        SA, SB = (certS, certC) if name == 'ours, surface cover' else (A, B)   # certified boundaries
        ax = fig.add_subplot(gs[0, k])
        for j in range(0, len(traj), 60):
            ax.add_patch(Polygon(SB.world(traj[j, 1]), closed=True, fill=False, ec=col, lw=0.6, alpha=0.3))
        ax.add_patch(Polygon(GA, closed=True, color='0.6', lw=0))
        ax.add_patch(Polygon(SA.dense, closed=True, fill=False, ec='k', lw=0.9))
        draw_pose(ax, SB, GB, traj[-1, 1], col)
        ax.plot(traj[:, 1, 0], traj[:, 1, 1], color=col, lw=0.7)
        if name == 'circle':                                   # minimum enclosing circles of the barrier
            for S, xf in ((A, traj[-1, 0]), (B, traj[-1, 1])):
                c, r = C.min_enclosing_circle(S.dense)
                ax.add_patch(Circle(xf[:2] + C.rot(xf[2]) @ c, r, fill=False, ec=col, lw=0.9, ls='--'))
        ax.set_xlim(-1.65, 4.3)
        ax.set_ylim(-1.32, 1.32)
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [m]')
        if k == 0:
            ax.set_ylabel('$y$ [m]')
        ax.set_title(f"({'ab'[k]}) {name}")
    t = np.arange(900) * 0.01
    axh = fig.add_subplot(gs[1, 0])
    axd = fig.add_subplot(gs[1, 1])
    for name, (L, traj, depth, col) in runs.items():
        n = len(depth)
        if name != 'circle':
            axh.plot(t[:n], 1e3 * np.array(L['h_t']), color=col, label=name, **STYLE[name])
        axd.plot(t[:n], depth, color=col, label=name, **STYLE[name])
    for a_ in (axh, axd):
        a_.spines[['top', 'right']].set_visible(False)
    axh.axhline(0, color='0.6', lw=0.6, zorder=0)
    axh.set_ylim(-25, 80)
    axh.set_xlim(0, 9)
    axh.set_xlabel('$t$ [s]')
    axh.set_ylabel('$h$ [mm]')
    axh.set_title('(c) barrier value')
    axh.legend(loc='upper right', frameon=False, fontsize=6)
    axd.axhline(0, color='0.6', lw=0.6, zorder=0)
    axd.axhline(SEAT_DEPTH, color='0.6', lw=0.6, zorder=0)
    axd.text(8.9, SEAT_DEPTH + 0.07, 'full seating', ha='right', fontsize=6)
    axd.text(8.9, 0.06, 'mouth', ha='right', fontsize=6)
    axd.set_xlim(0, 9)
    axd.set_ylim(-2.5, SEAT_DEPTH + 0.35)
    axd.set_xlabel('$t$ [s]')
    axd.set_ylabel('nose depth [m]')
    axd.set_title('(d) depth of the nose tip past the mouth')
    axd.legend(loc='lower right', frameon=False, fontsize=6, ncol=2, columnspacing=1.0, handlelength=1.5)
    fig.savefig(os.path.join(FIGS, 'dock_paper.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('>>> tro/figs/dock_paper.png')


def fig_five(construction='cover'):
    import json
    from cspace_cbf_5robots import five_robot_setup
    # configuration chosen by five_jitter_search.py (2x shapes, perturbed antipodal swap)
    # paper version: single integrators, perturbed antipodal swap (five_jitter_si.json)
    cfg = json.load(open(os.path.join(os.path.dirname(__file__), 'results', 'five_jitter_si.json')))
    GT, shapes, fields, starts, goals, _ = five_robot_setup(
        scale=cfg['scale'], radius=cfg['radius'], jitter_seed=cfg['chosen']['jitter_seed'],
        jitter=tuple(cfg['jitter']), goal_shift=cfg.get('goal_shift'), n_robots=cfg.get('n_robots', 5))
    steps, lim = cfg['steps'], cfg['radius'] + 1.5
    if cfg.get('starts') is not None:                  # explicit poses (four_touch_search.py)
        starts, goals = np.array(cfg['starts']), np.array(cfg['goals'])
    N = len(shapes)
    # Reuse the exact GIF run when its geometry, poses, and control horizon match.
    saved_path = os.path.join(os.path.dirname(__file__), 'results', 'five_random_trajectory.npz')
    L = None
    if os.path.exists(saved_path):
        with np.load(saved_path, allow_pickle=False) as saved:
            matches = (len(saved["x"]) == steps and float(saved["dt"]) == 0.01
                       and np.array_equal(saved['x'][0], starts)
                       and np.array_equal(saved['goals'], goals)
                       and all(np.array_equal(saved[f'ctrl_{i}'], shapes[i].ctrl)
                               and np.array_equal(saved[f'GT_{i}'], GT[i]) for i in range(N)))
            if matches:
                L = dict(traj=saved['x'], h_t=saved['h'], d_t=saved['d'], min_gt=saved['d'].min())
    # surface-cover run of the same swap (prototype_3d/five_cover.py swap): trajectory, min h, min distance;
    # colored outlines are then the certified level curves of the bodies' fields
    cover_path = os.path.join(os.path.dirname(__file__), 'results', 'five_cover_swap_summed_5.npz')
    if construction == 'cover':
        with np.load(cover_path) as z:
            assert np.allclose(z['traj'][0], starts) and np.allclose(z['goals'], goals)
            L = dict(traj=z['traj'], h_t=z['h'], d_t=z['d'], min_gt=z['d'].min())
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'prototype_3d'))
        import franka3d as F3
        from dock_cover import body_field, tight_level_2d

        class Level:
            def __init__(self, G):
                f = body_field(G, .01, .10)
                l_ = tight_level_2d(f, F3.hessian_majorant(f), G)
                xs_ = np.arange(f.lo[0], f.lo[0] + f.K[0] * f.h, .002)
                ys_ = np.arange(f.lo[1], f.lo[1] + f.K[1] * f.h, .002)
                XX, YY = np.meshgrid(xs_, ys_, indexing='ij')
                V = f.eval(np.c_[XX.ravel(), YY.ravel(), np.zeros(XX.size)]).reshape(XX.shape)
                cs = plt.figure().add_subplot().contour(xs_, ys_, V.T, levels=[l_])
                self.dense = max(cs.allsegs[0], key=len)
                plt.close(cs.axes.figure)

            def world(self, x):
                return self.dense @ C.rot(x[2]).T + x[:2]

        shapes = [Level(G) for G in GT]
    if L is None:
        L = simulate_team(GT, shapes, fields, starts, goals, "cspace", steps=steps, record=True,
                          dyn=cfg.get("dyn", "si"), v_max=cfg.get("v_max", 1.0), w_max=cfg.get("w_max", 2.0))
    print(f"five ({construction}): min GT {L['min_gt']:+.4f} min h {np.nanmin(L['h_t']):+.5f}")
    traj = np.array(L['traj'])
    cols = plt.cm.tab10(np.arange(N))
    fig = plt.figure(figsize=(3.5, 4.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.62)
    ax = fig.add_subplot(gs[0])
    for i in range(N):                              # layers: faded start < fading paths < solid final poses
        draw_pose(ax, shapes[i], GT[i], traj[0, i], cols[i], alpha=0.35, zorder=1)
        fading_path(ax, traj[:, i, :2], cols[i], lw=1.8, zorder=2)
        draw_pose(ax, shapes[i], GT[i], traj[-1, i], cols[i], zorder=3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$y$ [m]')
    ax.set_title('(a) start (faded) and final poses')
    axd = fig.add_subplot(gs[1])
    t = np.arange(len(L['d_t'])) * 0.01
    axd.axhline(0, color='0.75', lw=0.6, zorder=0)
    axd.plot(t, 1e3 * np.minimum(np.array(L['d_t']), 1.0), color='#2b2b2b', lw=1.2, label='min. true distance')
    axd.plot(t, 1e3 * np.array(L['h_t']), color='#119c99', lw=1.2, ls='--', label='min. active $h$')
    axd.spines[['top', 'right']].set_visible(False)
    axd.set_ylim(-5, 80)
    axd.set_xlim(0, t[-1])
    axd.set_xlabel('$t$ [s]')
    axd.set_ylabel('[mm]')
    axd.set_title('(b) safety margins', pad=13)
    axd.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, fontsize=6, handlelength=2.2, borderaxespad=0.2, columnspacing=1.0)
    fig.savefig(os.path.join(FIGS, 'five_paper.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('>>> tro/figs/five_paper.png')


def fig_overview():
    """(a) first contact of the two bodies, (b) the surface cover of B_i's certified curve near B_j,
    (c) the tabulated construction: stacked C-space slices with the certified level set."""
    import torch
    from cspace_experiments import closest_pair
    fig = plt.figure(figsize=(7.16, 2.5))
    # (a) two random non-convex bodies (as in the multi-robot experiments): slide B_i (rotated by 60 deg)
    # towards B_j along -x until first contact
    from cspace_cbf_5robots import five_robot_setup, pair_field
    rGT, rshapes, _, _, _, _ = five_robot_setup(scale=1.0)
    Ai, Bj = rshapes[3], rshapes[0]
    # (c) uses the certified configuration-space field of this same pair (B_i moving relative to B_j), cached
    A, B = Ai, Bj
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'overview_pair_3_0.npz')
    if os.path.exists(cache):
        z = np.load(cache)
        D, l_star = float(z['D']), float(z['l'])
        f = C.Field3D('bspline', D, int(z['Kxy']), int(z['Kth']))
        f.set_W(z['W'])
    else:
        f, D, l_star, ok, info = pair_field(A, B)
        assert ok, info
        np.savez(cache, W=f.W, D=D, l=l_star, Kxy=f.bx.K, Kth=f.bt.K)
    print(f'overview pair: D {D}, l* {l_star:+.4f}')
    # three equal panels; titles on one baseline. Colors: gray = reference / ground truth (thin, on top),
    # orange = certified level set and the contact point (complementary to the blue body)
    REF, CERT = '0.3', '#eb6834'
    W3, TY = 1. / 3, .93
    ax = fig.add_axes([.01, .05, W3 - .02, .82])
    th, ty = np.pi / 3, 0.5          # offset chosen so that the first contact is a single point
    tx = next(x for x in np.arange(2.2, 0.0, -0.001)
              if C.true_distance(Ai.world(np.array([x, ty, th])), Bj.dense) < 1e-3)
    xa = np.array([tx, ty, th])
    _, pa, pb, _ = closest_pair(torch.as_tensor(Ai.world(xa), device=C.DEV), torch.as_tensor(Bj.dense, device=C.DEV))
    ax.add_patch(Polygon(Bj.dense, closed=True, fc='0.85', ec='0.35', lw=1.0))
    ax.add_patch(Polygon(Ai.world(xa), closed=True, fc='#c6dbef', ec=BLUE, lw=1.0))
    ax.plot(*pb, marker='*', color='crimson', ms=9, zorder=5)
    lo = np.minimum(Bj.dense.min(0), Ai.world(xa).min(0)) - .1
    hi = np.maximum(Bj.dense.max(0), Ai.world(xa).max(0)) + .1
    # label above or below both bodies, arrow pointing at the star through free space only: the whole arrow,
    # tip included, keeps >= 4 cm from B_i and B_j (the bodies nearly touch around the star, so the tip
    # stops where the gap between them has opened); the label box does not overlap either body
    from matplotlib.path import Path as MPath
    bodies = [Bj.dense, Ai.world(xa)]
    CLR, LW_, LH_ = .04, .72, .09                       # clearance (incl. line widths); label width, height [m]

    def seg_clear(p, q):
        s = np.linspace(0, 1, 300)[:, None] * (q - p) + p
        if any(MPath(P).contains_points(s).any() for P in bodies):
            return -1.
        return min(np.min(np.linalg.norm(s[:, None] - P[None], axis=2)) for P in bodies)

    def box_free(x, y0, y1):
        g = np.stack(np.meshgrid(np.linspace(x - LW_ / 2, x + LW_ / 2, 25), np.linspace(y0, y1, 5)), -1).reshape(-1, 2)
        return not any(MPath(P).contains_points(g).any() for P in bodies) and \
            min(np.min(np.linalg.norm(g[:, None] - P[None], axis=2)) for P in bodies) > CLR

    best = None
    for side, y in (('below', lo[1] + .02), ('above', hi[1] - .02)):
        for x in np.linspace(lo[0] + LW_ / 2, hi[0] - LW_ / 2, 50):
            tail = np.array([x, y])
            if not box_free(x, *((y - LH_, y) if side == 'below' else (y, y + LH_))):
                continue
            u = (tail - pb) / np.linalg.norm(tail - pb)
            for stop in np.arange(.04, .26, .01):          # shortest stop at which the arrow is clear
                if seg_clear(pb + stop * u, tail) >= CLR:
                    if best is None or stop < best[0]:
                        best = (stop, side, tail, pb + stop * u)
                    break
    assert best is not None
    stop, side, tail, tip = best
    print(f'overview (a): label {side}, arrow stops {100 * stop:.0f} cm from the contact point')
    ax.annotate(r'$b(\tau)=\mathbf{t}+\mathbf{R}(\theta)a(\sigma)$', xy=tip, xytext=tail, ha='center',
                va='top' if side == 'below' else 'bottom', fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6, relpos=(0.5, 1.0 if side == 'below' else 0.0),
                                shrinkA=0, shrinkB=0))
    ax.text(*(Bj.dense.mean(0) - [0.1, 0.05]), '$\\mathcal{B}_j$', fontsize=8)
    ax.text(*(Ai.world(xa).mean(0) - [0.05, 0.05]), '$\\mathcal{B}_i$', fontsize=8, color=BLUE)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(*((lo[1] - .3, hi[1]) if side == 'below' else (lo[1], hi[1] + .2)))
    ax.set_aspect('equal')
    ax.set_axis_off()
    fig.text(W3 / 2, TY, '(a) first contact', ha='center')
    # (b) surface cover (Sec. IV): B_i moved 3 cm back from the contact of (a), zoom on the contact region.
    #     One field per body fitted to its ground-truth polygon, tightest certified levels (Prop. 1), and the
    #     true 5-mm Bernstein cover of B_i's certified curve S_i; dots: boxes whose coefficient constraints are active.
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'prototype_3d'))
    import franka3d as F3
    from dock_cover import body_field, tight_level_2d, cover_2d, CORNERS2
    fi, fj = body_field(rGT[3], .01, .1), body_field(rGT[0], .01, .1)
    Mi, Mj = F3.hessian_majorant(fi), F3.hessian_majorant(fj)
    li, lj = tight_level_2d(fi, Mi, rGT[3]), tight_level_2d(fj, Mj, rGT[0])
    PC, sz = cover_2d(fi, li, .005)
    print(f'overview (b): levels {1e3 * li:.2f} / {1e3 * lj:.2f} mm, {len(PC)} boxes of {1e3 * sz:g} mm')
    # first contact of the ground truths along the same approach, then 3 cm back
    Gi_w = lambda x: rGT[3] @ C.rot(th).T + np.array([x, ty])
    xg = next(x for x in np.arange(tx + .1, tx - .2, -.0005) if C.true_distance(Gi_w(x), rGT[0]) < 5e-4)
    pc = rGT[0][np.argmin(np.linalg.norm(rGT[0][:, None] - Gi_w(xg)[None], axis=2).min(axis=1))]
    xb = np.array([xg + .03, ty, th])
    Rb = C.rot(xb[2])
    cen, HW = pc + np.array([.015, 0.]), .12                    # zoom window (center, half width)

    def level_curve(f, l, lo_, hi_, res=.001):
        xs_, ys_ = np.arange(lo_[0], hi_[0], res), np.arange(lo_[1], hi_[1], res)
        XX, YY = np.meshgrid(xs_, ys_, indexing='ij')
        V = f.eval(np.c_[XX.ravel(), YY.ravel(), np.zeros(XX.size)]).reshape(XX.shape)
        cs = plt.figure().add_subplot().contour(xs_, ys_, V.T, levels=[l])
        segs = cs.allsegs[0]
        plt.close(cs.axes.figure)
        return segs

    ax = fig.add_axes([W3 + .01, .05, W3 - .02, .82])
    ax.add_patch(Polygon(rGT[0], closed=True, fc='0.85', ec='none'))
    ax.add_patch(Polygon(rGT[3] @ Rb.T + xb[:2], closed=True, fc='#c6dbef', ec='none'))
    for seg in level_curve(fj, lj, cen - HW - .02, cen + HW + .02):
        ax.plot(seg[:, 0], seg[:, 1], color=CERT, lw=1.2)
    for seg in level_curve(fi, li, fi.lo[:2], (fi.lo + fi.K * fi.h)[:2]):
        w_ = seg @ Rb.T + xb[:2]
        ax.plot(w_[:, 0], w_[:, 1], color=BLUE, lw=0.8)
    Cw = PC @ Rb.T + xb[:2]
    vis = np.flatnonzero(np.all(np.abs(Cw - cen) < HW + .01, axis=1))
    r_ = sz * np.sqrt(2) / 2
    for k in vis:
        corners = (PC[k] + CORNERS2[[0, 2, 3, 1]] * sz) @ Rb.T + xb[:2]
        ax.add_patch(Polygon(corners, closed=True, fill=False, ec=BLUE, lw=0.35))
    # dots: centers of the boxes whose Bernstein coefficient constraints are active (lower bound of the lifted
    # barrier below eta = 5 cm; summed-field barriers, prototype_3d/summed.py)
    import summed as SM
    body_i = SM.prepare_body(fi, li, PC, sz, 2)
    res = SM.rows(body_i, vis, Rb, Cw[vis], fj, lj, SM.field_bound(fj, r_), np.zeros((len(vis), 1, 2)),
                  np.zeros((1, 2, 2)), np.zeros((len(vis), 3)), .05, np.zeros(1), np.zeros(3), prune_rows=False,
                  depth=0)
    act = Cw[vis[res[5]]]
    ax.plot(act[:, 0], act[:, 1], 'o', color='0.15', ms=0.9, mew=0)
    ax.set_xlim(cen[0] - HW, cen[0] + HW)
    ax.set_ylim(cen[1] - HW, cen[1] + HW)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color('0.6'); sp.set_linewidth(0.6)
    ax.text(cen[0] - HW + .01, cen[1] - HW + .01, r'$\mathcal{B}_j$', fontsize=8)
    ax.text(cen[0] + HW - .025, cen[1] + HW - .03, r'$\mathcal{B}_i$', fontsize=8, color=BLUE)
    fig.text(1.5 * W3, TY, r'(b) cover of $\mathcal{S}_i$ (zoom, 5-mm boxes)', ha='center')
    # (c) compiled construction: stacked C-space slices
    ax = fig.add_axes([2 * W3 - .025, .12, W3 - .01, .8], projection='3d')   # room for the theta label
    ax.computed_zorder = False                  # thin reference drawn over the wider certified curve
    xs = np.linspace(-D, D, 200)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    for thk in np.linspace(0, TWO_PI, 8, endpoint=False):
        g, S = C.cspace_sdf_slices(A, B, [thk], D, res=0.01)
        Z = f.eval(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, thk)])).reshape(X.shape)
        for field, xx, lev, col, lw in ((S[0], g, 0.0, REF, 0.6), (Z, xs, l_star, CERT, 1.5)):
            cs = plt.figure().add_subplot().contour(xx, xx, field.T if field is S[0] else field.T, levels=[lev])
            for seg in cs.allsegs[0]:
                ax.plot(seg[:, 0], seg[:, 1], np.full(len(seg), np.degrees(thk)), color=col, lw=lw,
                        zorder=2 if col == REF else 1)
            plt.close(cs.axes.figure)
    ax.set_xlabel('$t_x$', labelpad=-8)
    ax.set_ylabel('$t_y$', labelpad=-8)
    ax.set_zlabel(r'$\theta$ [deg]', labelpad=-6)
    ax.set_zticks([0, 100, 200, 300])
    ax.tick_params(pad=-3)
    ax.view_init(elev=22, azim=-60)
    fig.text(2.5 * W3, TY, r'(c) compiled field over $SE(2)$', ha='center')
    fig.savefig(os.path.join(FIGS, 'overview.png'), dpi=300)
    plt.close(fig)
    print('>>> tro/figs/overview.png')


if __name__ == '__main__':
    import sys
    todo = sys.argv[1:] or ['dock', 'five', 'overview']
    for w in todo:
        {'dock': fig_dock, 'dock2': fig_dock2, 'five': fig_five, 'overview': fig_overview}[w]()
