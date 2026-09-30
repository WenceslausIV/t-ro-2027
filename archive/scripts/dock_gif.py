"""
Docking GIF with our B-spline barrier (single certified level l*): the vehicle approaches, rotates,
and seats its nose in the tapered receptacle.

    python dock_gif.py      -> dock_craft.gif, dock_craft_final.png
"""
import numpy as np

import cspace_sdf_cbf_compare as C
from cspace_experiments import DOCK2_GOAL, DOCK2_START, dock2_setup, simulate_team
from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter             # noqa: E402
from matplotlib.patches import Polygon                                   # noqa: E402

STEPS, DT, EVERY, HOLD = 900, 0.01, 8, 25
COL = '#1f77b4'


def main():
    GA, GB, A, B, fields, meta = dock2_setup()
    f, D, l_star = fields['bspline']
    L = simulate_team([GA, GB], [A, B], {(0, 1): fields['bspline']}, DOCK2_START, DOCK2_GOAL, 'cspace',
                      steps=STEPS, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
    traj = np.array(L['traj'])[:, 1]
    tip = traj[:, :2] + NOSE_TIP_X * np.stack([np.cos(traj[:, 2]), np.sin(traj[:, 2])], axis=1)
    depth = MOUTH_X - tip[:, 0]
    gap = np.array([C.true_distance(GA, GB @ C.rot(x[2]).T + x[:2]) for x in traj])
    print(f"depth {depth[-1]:.3f} m (full seating {SEAT_DEPTH:.3f} m), final gap {1e3 * gap[-1]:.1f} mm, "
          f"min gap {1e3 * gap.min():.1f} mm, min h {1e3 * np.nanmin(L['h_t']):.2f} mm")

    plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
    fig = plt.figure(figsize=(9, 6.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.3, 1.0], hspace=0.22)
    ax = fig.add_subplot(gs[0])
    ax.add_patch(Polygon(GA, closed=True, color='0.6', lw=0))
    ax.add_patch(Polygon(A.dense, closed=True, fill=False, ec='k', lw=0.9))
    pg = ax.add_patch(Polygon(GB, closed=True, color='0.6', alpha=0.8, lw=0))
    pf = ax.add_patch(Polygon(B.dense, closed=True, fill=False, ec=COL, lw=1.2))
    tr, = ax.plot([], [], color=COL, lw=0.8)
    ax.set_xlim(-1.0, 4.3)
    ax.set_ylim(-1.1, 1.25)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'B-spline configuration-space barrier, certified level $l^\\star$ = {1e3 * l_star:.1f} mm')
    az = ax.inset_axes([0.58, 0.02, 0.41, 0.38])                      # zoom on the receptacle
    az.add_patch(Polygon(GA, closed=True, color='0.6', lw=0))
    az.add_patch(Polygon(A.dense, closed=True, fill=False, ec='k', lw=0.7))
    zg = az.add_patch(Polygon(GB, closed=True, color='0.6', alpha=0.8, lw=0))
    zf = az.add_patch(Polygon(B.dense, closed=True, fill=False, ec=COL, lw=1.0))
    az.set_xlim(MOUTH_X - SEAT_DEPTH - 0.08, MOUTH_X + 0.12)
    az.set_ylim(-0.22, 0.22)
    az.set_aspect('equal')
    az.set_xticks([])
    az.set_yticks([])
    for s in az.spines.values():
        s.set_color(COL)
    txt = ax.text(0.02, 0.97, '', transform=ax.transAxes, va='top', fontsize=9)

    axd = fig.add_subplot(gs[1])
    t = np.arange(STEPS) * DT
    axd.plot(t[:len(depth)], depth, color=COL, lw=0.6, alpha=0.25)
    ln, = axd.plot([], [], color=COL, lw=1.6)
    axd.axhline(0, color='k', lw=0.5, ls=':')
    axd.axhline(SEAT_DEPTH, color='k', lw=0.5, ls='--')
    axd.text(t[-1], SEAT_DEPTH + 0.05, 'full seating', ha='right', fontsize=8)
    axd.text(t[-1], 0.05, 'mouth', ha='right', fontsize=8)
    axd.set_xlim(0, t[-1])
    axd.set_ylim(-2.5, SEAT_DEPTH + 0.3)
    axd.set_xlabel('$t$ [s]')
    axd.set_ylabel('nose depth [m]')

    n_frames = STEPS // EVERY + HOLD

    def update(fr):
        j = min(fr * EVERY, len(traj) - 1)
        x = traj[j]
        for p in (pg, zg):
            p.set_xy(GB @ C.rot(x[2]).T + x[:2])
        for p in (pf, zf):
            p.set_xy(B.world(x))
        tr.set_data(traj[:j + 1, 0], traj[:j + 1, 1])
        txt.set_text(f"nose depth {depth[j]:+.3f} m / {SEAT_DEPTH:.3f} m\ntrue gap {1e3 * gap[j]:5.1f} mm")
        ln.set_data(t[:j + 1], depth[:j + 1])
        return []

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    anim.save('dock_craft.gif', writer=PillowWriter(fps=20), dpi=80)
    update(n_frames - 1)
    fig.savefig('dock_craft_final.png', dpi=110, bbox_inches='tight')
    print('>>> dock_craft.gif, dock_craft_final.png')


if __name__ == '__main__':
    main()
