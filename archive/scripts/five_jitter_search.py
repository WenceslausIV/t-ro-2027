"""
Fig. 6 setup search: five random shapes at twice the size, antipodal swap on a circle with slightly
perturbed start/goal angles; picks the first perturbation in which every robot reaches its goal
without collision (the exact antipodal swap deadlocks at the center).

    python five_jitter_search.py   -> results/five_jitter.json
"""
import json
import time

import numpy as np

from cspace_cbf_5robots import five_robot_setup
from cspace_experiments import dump, simulate_team

SCALE, RADIUS, JITTER, STEPS = 2.0, 4.0, (0.30, 0.30, 0.6), 3000   # start angle, goal angle [rad], start/goal radius [m]
DYN, V_MAX, W_MAX = 'uni', 1.0, np.pi / 2      # unicycle, as in the unicycle experiments of the paper
GOAL_SHIFT = None     # None: (perturbed) antipodal goals; 2: pentagram crossings
N_ROBOTS = 4


def main(max_tries=60):
    GT, shapes, fields, starts, goals, meta = five_robot_setup(scale=SCALE, radius=RADIUS, n_robots=N_ROBOTS)
    print('pairs:', {k: round(v['l'] * 1e3, 1) for k, v in meta['pairs'].items()}, flush=True)
    tried = []
    for js in range(max_tries):
        _, _, _, starts, goals, _ = five_robot_setup(scale=SCALE, radius=RADIUS, jitter_seed=js, jitter=JITTER,
                                                     goal_shift=GOAL_SHIFT, n_robots=N_ROBOTS)
        t0 = time.perf_counter()
        L = simulate_team(GT, shapes, fields, starts, goals, 'cspace', dyn=DYN, steps=STEPS,
                          v_max=V_MAX, w_max=W_MAX)
        ok = L['reached'] is not None and L['min_gt'] > 0
        tried.append(dict(jitter_seed=js, reached=L['reached'], min_gt=L['min_gt'], min_h=L['min_h'],
                          slack=L['slack_steps']))
        print(f"jitter seed {js}: reached {L['reached']}, min gt {1e3 * L['min_gt']:.1f} mm, "
              f"min h {1e3 * L['min_h']:.2f} mm ({time.perf_counter() - t0:.0f} s)", flush=True)
        if ok:
            break
    dump('five_jitter', dict(scale=SCALE, radius=RADIUS, jitter=JITTER, steps=STEPS, dyn=DYN,
                             v_max=V_MAX, w_max=W_MAX, goal_shift=GOAL_SHIFT, n_robots=N_ROBOTS, tried=tried,
                             chosen=tried[-1] if ok else None, offsets=meta['offsets'], pairs=meta['pairs']))


if __name__ == '__main__':
    main()
