"""
GIF of the four-robot unicycle position exchange chosen for review (four_touch_search.poses(perm, seed)).
Red segments join robot pairs whose barrier is within 5 mm of zero (safety filter acting).

    python four_touch_gif.py [perm0 perm1 perm2 perm3 seed]   -> four_robots_uni.gif
"""
import itertools
import sys

import numpy as np

import cspace_sdf_cbf_compare as C
from cspace_cbf_5robots import five_robot_setup
from cspace_experiments import simulate_team
from four_touch_search import DYN, N, RADIUS, SCALE, STEPS, V_MAX, W_MAX, poses

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter             # noqa: E402
from matplotlib.patches import Polygon                                   # noqa: E402

EVERY, HOLD, NEAR = 5, 20, .005


def pair_h(fields, x):
    out = {}
    for i, j in itertools.combinations(range(len(x)), 2):
        f, D, l = fields[(i, j)]
        t = C.rot(-x[j, 2]) @ (x[i, :2] - x[j, :2])
        if np.any(np.abs(t) >= D - C.ACT_DELTA):
            continue
        out[(i, j)] = float(f.eval(np.r_[t, (x[i, 2] - x[j, 2]) % C.TWO_PI])[0] - l)
    return out


def main():
    args = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else [1, 3, 0, 2, 9]
    perm, seed = tuple(args[:4]), args[4]
    GT, shapes, fields, _, _, _ = five_robot_setup(scale=SCALE, radius=RADIUS, n_robots=N)
    starts, goals = poses(perm, seed)
    L = simulate_team(GT, shapes, fields, starts, goals, 'cspace', dyn=DYN, steps=STEPS, v_max=V_MAX,
                      w_max=W_MAX, record=True)
    traj = np.array(L['traj'])
    print(f"perm {perm} seed {seed}: reached {L['reached']}, min true dist {1e3 * L['min_gt']:.1f} mm")
    cols = plt.cm.tab10(np.arange(N))
    lim = RADIUS + 1.6
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    gt_p, fit_p, trails = [], [], []
    for i in range(N):
        ax.plot(goals[i, 0], goals[i, 1], 'x', color=cols[i], ms=9, mew=2)
        gt_p.append(ax.add_patch(Polygon(GT[i], closed=True, color='0.6', lw=0)))
        fit_p.append(ax.add_patch(Polygon(shapes[i].dense, closed=True, fill=False, ec=cols[i], lw=1.3)))
        trails.append(ax.plot([], [], color=cols[i], lw=.8)[0])
    links = [ax.plot([], [], color='red', lw=2)[0] for _ in range(6)]
    txt = ax.text(.02, .97, '', transform=ax.transAxes, va='top')
    n_frames = len(traj) // EVERY + HOLD

    def update(fr):
        k = min(fr * EVERY, len(traj) - 1)
        x = traj[k]
        for i in range(N):
            gt_p[i].set_xy(GT[i] @ C.rot(x[i, 2]).T + x[i, :2])
            fit_p[i].set_xy(shapes[i].world(x[i]))
            trails[i].set_data(traj[:k + 1, i, 0], traj[:k + 1, i, 1])
        hs = pair_h(fields, x)
        act = [(i, j) for (i, j), h in hs.items() if h < NEAR]
        for m, ln in enumerate(links):
            if m < len(act):
                i, j = act[m]
                ln.set_data([x[i, 0], x[j, 0]], [x[i, 1], x[j, 1]])
            else:
                ln.set_data([], [])
        hmin = min(hs.values()) if hs else np.nan
        txt.set_text(f"t = {k * .01:5.2f} s   min h = {1e3 * hmin:6.1f} mm\n"
                     f"red: pairs with h < {1e3 * NEAR:.0f} mm (filter acting)")
        return []

    FuncAnimation(fig, update, frames=n_frames, blit=False).save('four_robots_uni.gif',
                                                                 writer=PillowWriter(fps=20), dpi=80)
    print('>>> four_robots_uni.gif')


if __name__ == '__main__':
    main()
