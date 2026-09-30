"""
Static-obstacle example from a scene saved in the interactive demo (interactive/setups/*.json): the arm starts
inserted deep in the star tube (the saved configuration) and pulls out to the demo's home pose. Nominal:
damped-least-squares IK (as in the demo); filter: surface-cover barriers of all eight links against the star
tube (as in the demo). Writes startube_setup_results.json, startube_setup_<method>.npz and, with --fig,
tro/figs/star_paper.png (one row of time instants).

    python prototype_3d/startube_setup.py [setup.json] [--fig]
"""
import glob
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'interactive'))
import server as S                                      # noqa: E402
import franka3d as F                                    # noqa: E402
from sdf_cbf_utils import solve_ldp_qp                  # noqa: E402
from goal_ball import ball_triangles                    # noqa: E402
import summed as SM                                     # noqa: E402
import summed3d as S3                                   # noqa: E402

DT, STEPS = .01, int(os.environ.get('STAR_STEPS', 1500))
UMAX = np.full(7, 1.)                                   # |qd_j| <= QD_MAX (rad/s)
STAR_OBJECT = os.environ.get('STAR_OBJECT', 'star tube (rounded)')   # the paper's scene: rounded edges
# robot neutral gray (as the reference arm of the dual-arm figure), tube in its original blue, goal exit-sign green
ROBOT, TUBE = (.66, .68, .71), (.12, .47, .71)
VIEW = (18, -35)                                        # camera (elev, azim), also used for the goal ball's G


def load(path):
    d = json.load(open(path))
    o = d['objects'][0]
    return d, np.array(d['q1']), np.array(o['p']), np.array(o['R'])


def gap(q, links, B, RB, pB):
    """Whole-mesh lower bound using 4-mm triangles and only geometric AABB pruning."""
    T, _, _ = F.fk(q)
    cap = .05
    g = cap
    for L in links:
        Rl, tl = RB.T @ T[L['frame']][:3, :3], RB.T @ (T[L['frame']][:3, 3] - pB)
        P = L['gt'] @ Rl.T + tl
        outside = np.linalg.norm(np.maximum(np.maximum(B['blo'] - P, P - B['bhi']), 0.), axis=1)
        m = outside - L['gt_r'] < cap
        if m.any():
            g = min(g, float((B['sdf'](P[m]) - L['gt_r'][m]).min()))
    return g


def simulate(q0, goal_p, goal_R, links, B, RB, pB, method, waypoints=(), q_goal=None,
             steps=None, record_controls=False):
    """Nominal: IK toward each waypoint in turn (switch within 3 cm), then toward the goal pose."""
    q = q0.copy()
    log = dict(q=[], gap=[], h=[], rows=[], boxes=[], t=[], slack=0, reached=None)
    u, z_prev = None, None
    targets = [np.asarray(w, float) for w in waypoints] + [goal_p]
    it = 0
    for k in range(STEPS if steps is None else steps):
        T1, Z1, O1 = F.fk(q)
        p_ee = T1[8][:3, :3] @ S.TCP + T1[8][:3, 3]
        if it < len(targets) - 1 and np.linalg.norm(targets[it] - p_ee) < .03:
            it += 1
        J = np.vstack([np.cross(Z1, p_ee - O1).T, Z1.T])
        e = np.r_[np.clip(2. * (targets[it] - p_ee), -.5, .5), np.clip(1.5 * S._log(goal_R @ T1[8][:3, :3].T), -1., 1.)]
        Jp = J.T @ np.linalg.inv(J @ J.T + 1e-3 * np.eye(6))
        u_nom = np.clip(Jp @ e + (np.eye(7) - Jp @ J) @ (.5 * (S.Q_HOME - q)), -S.QD_MAX, S.QD_MAX)
        if q_goal is not None:                            # joint-space nominal toward a goal configuration
            u_nom = np.clip(1.5 * (q_goal - q), -S.QD_MAX, S.QD_MAX)
        t0 = time.perf_counter()
        rows, rhs = [], []
        lo = np.clip(-2. * (q - F.Q_MIN), -S.QD_MAX, 0.)
        hi = np.clip(2. * (F.Q_MAX - q), 0., S.QD_MAX)
        if method == 'ours':                              # summed-field barriers, Bernstein coefficient rows
            res = [S3.pair_rows(A, T1[A['frame']][:3, :3], T1[A['frame']][:3, 3], B, RB, pB,
                                [(Z1, O1, A['n_joints'], 1., 0)], F.ACT, UMAX) for A in links]
            Ar, Tr, Cr, nb, hl, Er = S3.stack(res, 7)
            u, s, z_prev = SM.solve(u_nom, Ar, Tr, Cr, Er, np.vstack([np.eye(7), -np.eye(7)]), np.r_[lo, -hi],
                                    z_prev)
            log['t'].append(time.perf_counter() - t0); log['slack'] += s > 0
            log['rows'].append(len(Cr)); log['boxes'].append(nb); log['h'].append(hl)
        else:
            if method == 'corner':                        # previous corner barriers (safe set depends on the boxes)
                Z1xO1 = np.cross(Z1, O1)
                for A in links:
                    res = S.body_rows(A, T1[A['frame']], Z1, Z1xO1, B, RB, pB)
                    if res is not None:
                        rows.append(res[0]); rhs.append(-S.GAMMA * res[1])
            n_bar = int(sum(len(r) for r in rows))
            G = np.vstack(rows + [np.eye(7), -np.eye(7)])
            hv = np.r_[np.concatenate(rhs) if rhs else np.zeros(0), lo, -hi]
            u, s = solve_ldp_qp(u_nom, G, hv, n_bar, u)
            log['t'].append(time.perf_counter() - t0); log['slack'] += s > 0
            log['rows'].append(n_bar); log['boxes'].append(0)
            log['h'].append(float(min(r.min() for r in rhs)) / -S.GAMMA if rhs else np.inf)
        log['gap'].append(gap(q, links, B, RB, pB)); log['q'].append(q.copy())
        if record_controls:
            log.setdefault('u', []).append(u.copy())
            log.setdefault('u_nom', []).append(u_nom.copy())
        if np.linalg.norm(goal_p - p_ee) < .01 and log['reached'] is None:
            log['reached'] = k * DT
        q = q + DT * u
    log['q'].append(q.copy())
    return log


def figure(links, B, RB, pB, qs, out, ncol=4, goal_p=None):
    import trimesh
    from franka_gif import decimate, shade
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    Vt, Ft = S.exact_surface(B)
    tube = (Vt @ RB.T + pB)[Ft]
    meshes = [decimate(L['V'], L['F'], .01) for L in links]
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    Vb, Fb = decimate(np.asarray(base.vertices), np.asarray(base.faces), .012)
    steps = [int(x) for x in np.linspace(0, len(qs) - 1, ncol)]
    plt.rcParams.update({'font.size': 8, 'font.family': 'serif'})
    fig = plt.figure(figsize=(7.16, 1.75))
    w = 1. / ncol
    for c, s in enumerate(steps):
        ax = fig.add_axes([c * w - .025, -.04, w + .05, 1.], projection='3d')
        T, _, _ = F.fk(qs[s])
        tris = [tube, Vb[Fb]]
        cols = [np.c_[shade(tube, TUBE), np.full(len(tube), .75)], np.c_[shade(Vb[Fb], ROBOT), np.ones(len(Fb))]]
        for L, (V, Fc) in zip(links, meshes):
            tri = (V @ T[L['frame']][:3, :3].T + T[L['frame']][:3, 3])[Fc]
            tris.append(tri); cols.append(np.c_[shade(tri, ROBOT), np.ones(len(tri))])
        if goal_p is not None:                  # goal: exit-sign green ball, G wrapped on its surface
            gt, gc = ball_triangles(goal_p, .045, *VIEW)
            tris.append(gt); cols.append(gc)
        ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols), edgecolor='none'))
        ax.set_xlim(-.5, .5); ax.set_ylim(-.45, .45); ax.set_zlim(.05, .95)
        ax.set_box_aspect((1., .9, .9)); ax.view_init(*VIEW); ax.set_axis_off()
        ax.set_title(f't = {s * DT:.1f} s', y=.9, fontsize=8)
    fig.savefig(out, dpi=300)
    print('>>>', out, steps)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    path = args[0] if args else sorted(glob.glob(os.path.join(HERE, 'interactive', 'setups', '*.json')))[-1]
    d, q0, pB, RB = load(path)
    sim = S.Sim()
    links, B = sim.links, sim.objlib[STAR_OBJECT]
    Th, _, _ = F.fk(S.Q_HOME)
    goal_p, goal_R = Th[8][:3, :3] @ S.TCP + Th[8][:3, 3], Th[8][:3, :3]
    rev = '--reverse' in sys.argv
    if rev:          # start threaded through the tube (end of the forward run), leave to the saved hand pose
        q0 = np.load(os.path.join(HERE, 'startube_forward_ours.npz'))['q'][-1]
        goal_p, goal_R = np.array(d['ee']['p']), np.array(d['ee']['R'])
    tag = 'reverse_native' if rev else 'setup_native'
    ax_ = RB[:, 0]                                   # tube axis (the hand is on the +axis side when threaded)
    wps = []
    q_goal = np.array(d['q1']) if rev else None      # reverse: joint-space nominal back to the saved configuration
    out = dict(setup=os.path.basename(path), direction=tag, goal=dict(p=goal_p.tolist()))
    if '--redraw' in sys.argv:                       # figure only, from the saved run
        q = np.load(os.path.join(HERE, f'startube_{tag}_ours.npz'))['q']
        v = np.linalg.norm(np.diff(q, axis=0), axis=1) / DT
        k_stop = int(np.flatnonzero(v > .02)[-1]) + 1
        figure(links, B, RB, pB, q[:k_stop + 1], os.path.join(os.path.dirname(HERE), 'tro', 'figs', 'star_paper.png'),
               goal_p=goal_p)
        return
    for A in links:
        S3.prep_link(A)
    S3.prep_field(B)
    methods = ('nominal', 'ours', 'corner') if '--corner' in sys.argv else ('nominal', 'ours')
    for m in methods:
        L = simulate(q0, goal_p, goal_R, links, B, RB, pB, m, waypoints=wps, q_goal=q_goal)
        r = dict(min_gap_mm=1e3 * min(L['gap']), reached=L['reached'], slack_steps=int(L['slack']),
                 min_h_mm=1e3 * min(L['h']), rows_max=int(max(L['rows'])), boxes_max=int(max(L['boxes'])),
                 t_p95_ms=1e3 * float(np.percentile(L['t'], 95)),
                 t_med_ms=1e3 * float(np.median(L['t'])), t_max_ms=1e3 * float(np.max(L['t'])),
                 unresolved_gap_steps=int(np.sum(np.array(L['gap']) <= 0)),
                 distance_scope='Whole-mesh bounds at evaluated states; nonpositive bounds do not prove collision')
        out[m] = r
        print(m, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        np.savez_compressed(os.path.join(HERE, f'startube_{tag}_{m}.npz'), q=np.array(L['q']), gap=np.array(L['gap']))
        if m == 'ours' and '--fig' in sys.argv:
            q = np.array(L['q'])
            v = np.linalg.norm(np.diff(q, axis=0), axis=1) / DT
            k_stop = int(np.flatnonzero(v > .02)[-1]) + 1              # show the moving interval
            figure(links, B, RB, pB, q[:k_stop + 1], os.path.join(os.path.dirname(HERE), 'tro', 'figs',
                                                                   'star_paper.png'), goal_p=goal_p)
    json.dump(out, open(os.path.join(HERE, f'startube_{tag}_results.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
