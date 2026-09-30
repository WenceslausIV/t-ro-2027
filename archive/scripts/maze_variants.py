"""
Try nominal-controller variants on the unicycle rocket maze (same walls, same certified field) and
report which reach the exit and how often the barrier acts. No figures are written.

    python maze_variants.py
"""
import json

import numpy as np

import maze_escape as M
import maze_reference_geometry as G
import reference_maze_escape as R
import rocket_geometry as Rocket


def route_points(path, step=.1):
    arc = np.r_[0, np.linalg.norm(np.diff(path, axis=0), axis=1).cumsum()]
    sub = path[np.searchsorted(arc, np.r_[np.arange(0, arc[-1], step), arc[-1]])[:-1].tolist() + [len(path) - 1]]
    hd = np.arctan2(*np.diff(sub, axis=0).T[::-1])
    return np.column_stack([sub, np.r_[hd[0], hd]])


def main():
    robot, gt, meta = Rocket.setup()
    walls = G.make_walls()
    path = R.plan(robot, walls, meta['width_m'] / 2)
    walls, center = R.smooth_corridor_walls(walls, path)
    fitted = [G.fit_wall(w)[0] for w in walls]
    f, fmeta = R.field(robot, fitted, False)
    wp = route_points(center)
    out = []
    for look in (.6, .8, 1.2):                         # wall following, no recovery (never backs up)
        name = f'center, look-ahead {look}, wall following, no recovery'
        try:
            _, res = M.simulate(f, fmeta['level'], gt, walls, waypoints=wp, domain=R.DOMAIN,
                                gap_function=G.physical_gaps, output_dir=R.OUT, cruise=True, steps=12000,
                                switch_radius=.40, dyn='uni', lookahead=look, recovery=False, wall_follow=.10)
            ok = True
        except RuntimeError:
            res = json.loads((R.OUT / 'last_result.json').read_text())
            ok = False
        out.append(dict(name=name, look=look, ok=ok, **{k: res.get(k) for k in (
            'reached_s', 'barrier_intervention_steps', 'steps', 'min_physical_gap_m', 'min_h_m', 'path_length_m')}))
        print(out[-1], flush=True)
    (R.OUT / 'variants.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
