"""
Option E: smaller rocket for the unicycle maze. For each rocket size, rebuild the pipeline exactly as
reference_maze_escape.main does (route, smooth corridor, wall fits, certified field) and run the
unicycle pure-pursuit controller without recovery for several look-ahead distances.

    python maze_variants_size.py   -> results/reference_maze_rocket/variants_size.json
"""
import json

import numpy as np

import maze_escape as M
import maze_reference_geometry as G
import reference_maze_escape as R
import rocket_geometry as Rocket
from maze_variants import route_points


def main():
    out = []
    for radius in (.26, .23):
        Rocket.RADIUS = radius
        Rocket.SCALE = radius / np.linalg.norm(Rocket.RAW - Rocket.CENTER, axis=1).max()
        robot, gt, meta = Rocket.setup()
        walls = G.make_walls()
        path = R.plan(robot, walls, meta['width_m'] / 2)
        walls, center = R.smooth_corridor_walls(walls, path)
        fitted = [G.fit_wall(w)[0] for w in walls]
        f, fmeta = R.field(robot, fitted, False)
        wp = route_points(center)
        for look in (.5, .8):
            try:
                _, res = M.simulate(f, fmeta['level'], gt, walls, waypoints=wp, domain=R.DOMAIN,
                                    gap_function=G.physical_gaps, output_dir=R.OUT, cruise=True, steps=12000,
                                    switch_radius=.40, dyn='uni', lookahead=look, recovery=False)
                ok = True
            except RuntimeError:
                res = json.loads((R.OUT / 'last_result.json').read_text())
                ok = False
            out.append(dict(radius=radius, length=meta['length_m'], width=meta['width_m'], level=fmeta['level'],
                            look=look, ok=ok, **{k: res.get(k) for k in (
                                'reached_s', 'barrier_intervention_steps', 'steps', 'min_physical_gap_m',
                                'min_h_m', 'path_length_m')}))
            print(out[-1], flush=True)
    (R.OUT / 'variants_size.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
