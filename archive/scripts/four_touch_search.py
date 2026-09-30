"""
Fig. 6 setup search: four unicycle robots (2x random shapes) exchange positions on a circle.
Goals are permutations of the start positions (the all-antipodal one, which jams at the center, is
excluded), with small start/goal noise. Accept the first run in which every robot reaches its goal
without collision AND every robot comes within TOUCH of the barrier boundary with at least one other
robot (so that the safety filter visibly acts on all four).

    python four_touch_search.py   -> results/five_jitter.json (used by make_paper_figs.py five)
"""
import itertools
import time

import numpy as np

from cspace_cbf_5robots import five_robot_setup
from cspace_experiments import dump, simulate_team

N, SCALE, RADIUS, STEPS = 4, 2.0, 4.0, 3000
DYN, V_MAX, W_MAX = 'uni', 1.0, np.pi / 2
NOISE_ANG, NOISE_RAD = 0.25, 0.6          # start/goal noise: angle [rad], radius [m]
TOUCH = 0.005                             # "touching": barrier value below 5 mm (red links in the GIF)


def poses(perm, seed):
    rng = np.random.default_rng(seed)
    ang = np.pi / 2 + 2 * np.pi * np.arange(N) / N
    a_s = ang + rng.uniform(-NOISE_ANG, NOISE_ANG, N)
    r_s = RADIUS + rng.uniform(-NOISE_RAD, NOISE_RAD, N)
    a_g = ang[list(perm)] + rng.uniform(-NOISE_ANG, NOISE_ANG, N)
    r_g = RADIUS + rng.uniform(-NOISE_RAD, NOISE_RAD, N)
    ps = np.stack([r_s * np.cos(a_s), r_s * np.sin(a_s)], axis=1)
    pg = np.stack([r_g * np.cos(a_g), r_g * np.sin(a_g)], axis=1)
    hd = np.arctan2(*(pg - ps).T[::-1])                  # every robot starts facing its goal
    starts, goals = np.column_stack([ps, hd]), np.column_stack([pg, hd])
    return starts, goals


def main():
    GT, shapes, fields, _, _, meta = five_robot_setup(scale=SCALE, radius=RADIUS, n_robots=N)
    perms = [p for p in itertools.permutations(range(N))
             if all(p[i] != i for i in range(N)) and p != (2, 3, 0, 1)]
    tried, chosen = [], None
    for seed in range(0, 40):
        for perm in perms:
            starts, goals = poses(perm, seed)
            t0 = time.perf_counter()
            L = simulate_team(GT, shapes, fields, starts, goals, 'cspace', dyn=DYN, steps=STEPS,
                              v_max=V_MAX, w_max=W_MAX, stall_stop=200)
            P = np.minimum(L['pair_min_h'], L['pair_min_h'].T)
            touch = [float(np.min(np.delete(P[i], i))) for i in range(N)]
            ok = L['reached'] is not None and L['min_gt'] > 0 and max(touch) <= TOUCH
            rec = dict(perm=list(perm), seed=seed, reached=L['reached'], min_gt=L['min_gt'],
                       min_h=L['min_h'], touch=touch, ok=ok)
            tried.append(rec)
            print(f"perm {perm} seed {seed}: reached {L['reached']}, min gt {1e3 * L['min_gt']:.1f} mm, "
                  f"closest-h per robot {[round(1e3 * t, 1) for t in touch]} mm "
                  f"({time.perf_counter() - t0:.0f} s){'  <== chosen' if ok else ''}", flush=True)
            if ok:
                chosen = dict(rec, starts=starts.tolist(), goals=goals.tolist())
                break
        if chosen:
            break
    if chosen is None:                                  # best near-success: all reached, closest to touching
        ok = [t for t in tried if t['reached'] is not None and t['min_gt'] > 0]
        if ok:
            best = min(ok, key=lambda t: max(t['touch']))
            st, go = poses(tuple(best['perm']), best['seed'])
            print('best near-success:', best, flush=True)
            dump('four_touch_best', dict(best, starts=st.tolist(), goals=go.tolist()))
    dump('four_touch' if chosen is None else 'five_jitter', dict(scale=SCALE, radius=RADIUS, steps=STEPS, dyn=DYN, v_max=V_MAX, w_max=W_MAX,
                             n_robots=N, jitter=[0.0, 0.0], touch=TOUCH, tried=tried, chosen=chosen,
                             starts=chosen['starts'] if chosen else None,
                             goals=chosen['goals'] if chosen else None,
                             offsets=meta['offsets'], pairs=meta['pairs']))


if __name__ == '__main__':
    main()
