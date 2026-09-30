"""Scan start poses for the star-tube insertion: keep those where the unfiltered motion collides and ours
reaches the goal inside the tube (for choosing the paper scenario)."""
import itertools
import numpy as np
import dolphin3d as Dp

links, obst = Dp.build_hoop()
for y, z in itertools.product((-.14, -.09, .09, .14), (.42, .47, .63, .68)):
    Dp.START_EE = np.array([.38, y, z])
    Dp.Q_START = np.array([0., -.35, 0., -2.55, 0., 2.2, .785])
    q0 = Dp.start_config()
    if np.linalg.norm(Dp.ee_pose(q0)[0] - Dp.START_EE) > .01:
        continue
    n = Dp.simulate(q0, links, obst, 'nominal', steps=500)
    if min(n['gap']) >= 0:
        print(y, z, 'nominal safe %.1f' % (1e3 * min(n['gap'])), flush=True)
        continue
    o = Dp.simulate(q0, links, obst, 'ours', steps=700)
    print(y, z, 'nominal %.1f' % (1e3 * min(n['gap'])), 'ours gap %.1f reached %s err %.0f' % (
        1e3 * min(o['gap']), o['reached'], 1e3 * np.linalg.norm(np.array(o['ee'])[-1] - Dp.GOAL_EE)), flush=True)
