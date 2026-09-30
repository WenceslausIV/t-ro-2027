"""
Unicycle rocket maze with a wider passage (separate output folder; the paper's single-integrator
maze results are not touched). No recovery maneuver: the rocket never backs up.

    python maze_uni_wide.py   -> results/maze_uni_wide/variants.json
"""
import json

import maze_escape as M
import maze_reference_geometry as G
import reference_maze_escape as R
import rocket_geometry as Rocket
from maze_variants import route_points

R.OUT = R.ROOT / 'results' / 'maze_uni_wide'
WIDTH_SCALE, WIDTH_MIN, SIGMA_CENTER = 2.0, .60, 1.0


def main():
    R.OUT.mkdir(parents=True, exist_ok=True)
    robot, gt, meta = Rocket.setup()
    walls = G.make_walls()
    path = R.plan(robot, walls, meta['width_m'] / 2)
    walls, center = R.smooth_corridor_walls(walls, path, sigma_center=SIGMA_CENTER,
                                            width_scale=WIDTH_SCALE, width_min=WIDTH_MIN)
    fitted = [G.fit_wall(w)[0] for w in walls]
    f, fmeta = R.field(robot, fitted, False)
    wp = route_points(center)
    out = []
    for look, wf in ((.8, 0.), (1.0, .15), (1.5, .15)):
        try:
            arrays, res = M.simulate(f, fmeta['level'], gt, walls, waypoints=wp, domain=R.DOMAIN,
                                     gap_function=G.physical_gaps, output_dir=R.OUT, cruise=True, steps=12000,
                                     switch_radius=.40, dyn='uni', lookahead=look, recovery=False, wall_follow=wf)
            ok = True
        except RuntimeError:
            res = json.loads((R.OUT / 'last_result.json').read_text())
            ok = False
        out.append(dict(look=look, wall_follow=wf, ok=ok, level=fmeta['level'], **{k: res.get(k) for k in (
            'reached_s', 'barrier_intervention_steps', 'steps', 'min_physical_gap_m',
            'verified_interval_gap_lower_m', 'min_h_m', 'path_length_m')}))
        print(out[-1], flush=True)
    (R.OUT / 'variants.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
