"""Shared maze simulation and the archived branched-wall setup.

Run: python maze_escape.py [--rebuild] (current reference-image maze).
The original obstacle-count study in results/static.json remains a separate experiment.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from matplotlib.path import Path as MplPath

import cspace_sdf_cbf_compare as C
import cspace_experiments as E
from sdf_cbf_utils import solve_ldp_qp

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'maze_escape'
FIG = ROOT / 'tro' / 'figs'
DOMAIN = 6.2
KXY, KTH = 160, 48
# Each group is one connected wall. Overlapping capsule beams form genuine
# non-convex crosses/Ts; the eight groups remain mutually disjoint.
WALL_GROUPS = [
    dict(name='bottom', beams=[((-4.1,-4.9),(4.15,-4.9),.22)]),
    dict(name='top', beams=[((-4.2,4.9),(3.85,4.9),.26)]),
    dict(name='left', beams=[((-5.05,-4.2),(-5.05,4.05),.23)]),
    dict(name='right', beams=[((5.,-4.1),(5.,2.2),.30)]),
    dict(name='long_cross', beams=[((-4.,-1.5),(.5,-1.5),.24),
                                  ((-2.,-3.2),(-2.,.3),.25)]),
    dict(name='left_facing_T', beams=[((3.25,-3.5),(3.25,-.15),.20),
                                     ((1.7,-2.7),(3.25,-2.7),.22)]),
    dict(name='downward_T', beams=[((-.2,2.5),(3.7,2.5),.29),
                                 ((1.8,.8),(1.8,2.5),.27)]),
    dict(name='upward_T', beams=[((-4.25,1.5),(-2.7,1.5),.20),
                               ((-3.6,1.5),(-3.6,3.55),.23)]),
]
WALLS = [beam for group in WALL_GROUPS for beam in group['beams']]
GROUP_IDS = [i for i,group in enumerate(WALL_GROUPS) for _ in group['beams']]
WAYPOINTS = np.array([
    [-4.1,-4.,0.], [.9,-4.,0.], [.9,-.4,np.pi/2],
    [-.25,1.05,np.pi], [-2.3,1.05,np.pi], [-1.25,3.7,np.pi/2],
    [5.3,3.7,0.],
])


class Capsule:
    """Exact capsule enclosed by six cubic Bezier spans (circular-arc approximation)."""
    def __init__(self, a, b, radius):
        self.a, self.b, self.radius = np.array(a), np.array(b), radius
        e = (self.b-self.a) / np.linalg.norm(self.b-self.a)
        n = np.array([-e[1], e[0]])
        k = 4*(np.sqrt(2)-1)/3
        def line(p, q):
            return np.array([p, (2*p+q)/3, (p+2*q)/3, q])
        def arc(center, u, v):
            return center + radius*np.array([u, u+k*v, v+k*u, v])
        self.bez = np.array([
            line(self.a-radius*n, self.b-radius*n),
            arc(self.b, -n, e), arc(self.b, e, n),
            line(self.b+radius*n, self.a+radius*n),
            arc(self.a, n, -e), arc(self.a, -e, -n),
        ])
        self.dense = np.einsum('si,kid->ksd', C.bern3(np.arange(80)/80), self.bez).reshape(-1, 2)
        # Exact capsule outline is used for the figure and physical distance audit.
        angles = np.linspace(-np.pi/2, np.pi/2, 100)
        self.outline = np.vstack([
            self.b + radius*(np.cos(angles)[:, None]*e + np.sin(angles)[:, None]*n),
            self.a + radius*(-np.cos(angles)[:, None]*e - np.sin(angles)[:, None]*n),
        ])


def physical_gaps(pose, robot_gt, walls):
    """Exact disjoint polygon/capsule distance; nonpositive means intersection.

    For disjoint planar segments their closest pair includes an endpoint.
    Explicit crossings and containment cover the remaining intersection cases.
    """
    p = robot_gt @ C.rot(pose[2]).T + pose[:2]
    q = np.roll(p, -1, axis=0)
    pq = q-p
    def cross(u, v):
        return u[..., 0]*v[..., 1] - u[..., 1]*v[..., 0]
    gaps = []
    for wall in walls:
        a, b, r = wall.a, wall.b, wall.radius
        ab = b-a
        t = np.clip(((p-a)@ab)/(ab@ab), 0, 1)
        d = np.linalg.norm(p-a-t[:, None]*ab, axis=1).min()
        for endpoint in (a, b):
            s = np.clip(np.sum((endpoint-p)*pq, axis=1)/np.sum(pq*pq, axis=1), 0, 1)
            d = min(d, np.linalg.norm(endpoint-p-s[:, None]*pq, axis=1).min())
        den = cross(pq, ab)
        mask = np.abs(den) > 1e-14
        u = np.divide(cross(a-p, ab), den, out=np.full(len(p), np.inf), where=mask)
        v = np.divide(cross(a-p, pq), den, out=np.full(len(p), np.inf), where=mask)
        crossing = np.any((u >= 0) & (u <= 1) & (v >= 0) & (v <= 1))
        if crossing or MplPath(p).contains_points([a, b]).any():
            d = 0.
        gaps.append(float(d-r))
    return np.array(gaps)


def validate_geometry(walls):
    # Verify the standard cubic quarter-circle encloses the true unit arc:
    # evaluate all extrema of ||p(t)||^2 on [0,1], plus tangent-angle monotonicity.
    k = 4*(np.sqrt(2)-1)/3
    quarter = np.array([[1, 0], [1, k], [k, 1], [0, 1]])
    assert np.all(np.diff(quarter[:,0]) <= 0) and np.all(np.diff(quarter[:,1]) >= 0)
    coeff = C.MB @ quarter
    norm2 = sum(np.polynomial.polynomial.polymul(c, c) for c in coeff.T)
    roots = np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(norm2))
    ts = [0., 1.] + [z.real for z in roots if abs(z.imag) < 1e-8 and 0 < z.real < 1]
    assert np.polynomial.polynomial.polyval(ts, norm2).min() >= 1-1e-12
    assert all(np.isfinite(w.dense).all() for w in walls)
    # Check group connectivity separately from separation BETWEEN groups.
    adjacency = np.eye(len(walls),dtype=bool)
    min_separation = np.inf
    for i, wi in enumerate(walls):
        for j in range(i+1,len(walls)):
            wj = walls[j]
            amin, amax = np.minimum(wi.a, wi.b), np.maximum(wi.a, wi.b)
            bmin, bmax = np.minimum(wj.a, wj.b), np.maximum(wj.a, wj.b)
            separation = np.linalg.norm(np.maximum(np.maximum(amin-bmax, bmin-amax), 0))-wi.radius-wj.radius
            if GROUP_IDS[i] == GROUP_IDS[j]:
                adjacency[i,j] = adjacency[j,i] = separation <= 0
            else:
                assert separation > 0, (wi.a, wj.a, separation)
                min_separation = min(min_separation, separation)
    for gid in range(len(WALL_GROUPS)):
        ids = np.flatnonzero(np.array(GROUP_IDS)==gid)
        visited = {int(ids[0])}
        while True:
            grown = visited | {int(j) for i in visited for j in ids if adjacency[i,j]}
            if grown == visited:
                break
            visited = grown
        assert visited == set(ids), 'Disconnected beams within a declared wall'
    # A long wall-axis crossing must be detected even when all robot vertices are outside.
    square = np.array([[-.1, -.1], [.1, -.1], [.1, .1], [-.1, .1]])
    probe = Capsule((-2., 0.), (2., 0.), .02)
    assert physical_gaps(np.zeros(3), square, [probe])[0] < 0
    assert abs(physical_gaps(np.array([0., 1., 0.]), square, [probe])[0]-.88) < 1e-10
    return float(min_separation)


def make_field(robot, walls, rebuild):
    signature = hashlib.sha256(json.dumps(dict(walls=WALLS, domain=DOMAIN, cells=[KXY,KTH],
        ctrl=robot.ctrl.tolist(), raster=.02, ntheta=96, version=1), sort_keys=True).encode()).hexdigest()
    cache = OUT / 'field.npz'
    meta_path = OUT / 'field.json'
    f = C.Field3D('bspline', DOMAIN, KXY, KTH)
    if cache.exists() and meta_path.exists() and not rebuild:
        meta = json.loads(meta_path.read_text())
        if meta['signature'] == signature and all(meta[k] for k in ('contact_ok','band_ok','regular_ok')):
            with np.load(cache) as data:
                f.set_W(data['W'])
            return f, meta
    start = time.perf_counter()
    th = np.arange(96)*C.TWO_PI/96
    g, slices = E.union_cspace_slices(robot, [w.dense for w in walls], th, DOMAIN+.1, res=.02)
    gxy = np.linspace(-DOMAIN, DOMAIN, 310)
    f.fit(C.sample_slices(g, slices.astype(float), gxy, gxy).transpose(1,2,0), gxy, gxy, th)
    fit_time = time.perf_counter()-start
    print(f'Field fitted in {fit_time:.2f} s; certifying contact...', flush=True)
    cert = C.Certifier(f, robot, E.Union(walls))
    level, ok, regions, cert_time = cert.certify(verbose=True)
    band_ok, band_min = cert.certify_band(level)
    reg_ok, reg_depth, reg_open = cert.certify_regular(level,max_depth=10)
    meta = dict(signature=signature, fit_s=fit_time, contact_s=cert_time, level=float(level),
        contact_ok=bool(ok), band_ok=bool(band_ok), band_min=float(band_min),
        regular_ok=bool(reg_ok), regular_depth=reg_depth, regular_open_boxes=reg_open,
        regions=regions, cells=[KXY,KXY,KTH], coefficients=f.n_coef)
    print(meta, flush=True)
    if not (ok and band_ok and reg_ok):
        np.savez_compressed(OUT/'uncertified_field.npz',W=f.W)
        (OUT/'uncertified_field.json').write_text(json.dumps(meta,indent=2))
        raise RuntimeError('Field certification incomplete; no paper result will be exported.')
    np.savez_compressed(cache, W=f.W)
    meta_path.write_text(json.dumps(meta, indent=2))
    return f, meta


def simulate(f, level, robot_gt, walls, waypoints=None, domain=None, gap_function=None,
             output_dir=None, cruise=False, steps=7500, switch_radius=.20, dyn='si', lookahead=.8, recovery=True, wall_follow=0.,
             rows_fn=None, summed_fn=None):
    """rows_fn(x) -> (rows (k, 3), h (k,)): barrier rows of the surface-cover construction instead of the field f.
    summed_fn(x) -> (A (k, 3), T (k, 3), C (k,), boxes, h_low, E (3, 3)) or None: summed-field barriers
    (prototype_3d/summed.py), Bernstein coefficient rows A u + T a + C >= 0 with a >= |E u| (single integrator only)."""
    waypoints = WAYPOINTS if waypoints is None else waypoints
    domain = DOMAIN if domain is None else domain
    gap_function = physical_gaps if gap_function is None else gap_function
    output_dir = OUT if output_dir is None else output_dir
    dt, vmax, wmax, gamma = .01, .8, 1.4, 5.
    Gbox, hbox = E.box_rows(1, dyn, vmax, wmax)
    min_row_near = np.inf                    # norm of the barrier row whenever h < 5 cm
    x = waypoints[0].copy()
    traj, controls, timings, barriers, residuals, nominal_residuals = [x.copy()], [], [], [], [], []
    n_active, n_boxes, z_sum = [], [], None
    gaps = [gap_function(x, robot_gt, walls)]
    rho = np.linalg.norm(robot_gt, axis=1).max()
    wp = 1
    slack_steps, unresolved, collisions = 0, 0, 0
    verified_lower = np.inf
    def audit(a, b, ga, gb, depth=0):
        nonlocal verified_lower, unresolved, collisions
        motion = np.linalg.norm(b[:2]-a[:2])+rho*abs(b[2]-a[2])
        lower = min(ga.min(), gb.min())-motion/2
        if min(ga.min(), gb.min()) <= 0:
            collisions += 1
        elif lower > 0:
            verified_lower = min(verified_lower, lower)
        elif depth >= 16:
            unresolved += 1
        else:
            mid = (a+b)/2
            gm = gap_function(mid, robot_gt, walls)
            audit(a, mid, ga, gm, depth+1)
            audit(mid, b, gm, gb, depth+1)
    reached = None
    # unicycle: pure pursuit of a carrot `lookahead` ahead of the progress along the waypoint polyline
    poly = waypoints[:, :2]
    seg_len = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.r_[0, seg_len.cumsum()]
    s_prog = 0.
    stall, backup, recoveries, careful = 0, 0, 0, 0

    def along(sv):
        return np.array([np.interp(sv, cum, poly[:, 0]), np.interp(sv, cum, poly[:, 1])])
    for step in range(steps):
        t0 = time.perf_counter()
        if wp < len(waypoints)-1 and np.linalg.norm(x[:2]-waypoints[wp,:2]) < switch_radius:
            wp += 1
        goal = waypoints[wp]
        d = goal[:2]-x[:2]
        v = 1.5*d
        if cruise and wp < len(waypoints)-1:
            v = vmax*d/max(np.linalg.norm(d),1e-12)
        v *= min(1., vmax/max(np.linalg.norm(v), 1e-12))
        if dyn == 'uni':                     # pure pursuit: steer toward the carrot, drive along the heading
            for k in range(len(seg_len)):    # projection onto the polyline near the current progress
                if cum[k+1] < s_prog-.5 or cum[k] > s_prog+2.:
                    continue
                t = np.clip((x[:2]-poly[k])@(poly[k+1]-poly[k])/max(seg_len[k]**2, 1e-12), 0, 1)
                sk = cum[k]+t*seg_len[k]
                if sk > s_prog and np.linalg.norm(x[:2]-along(sk)) < 1.:
                    s_prog = sk
            look = .3 if careful > 0 else lookahead         # shorter look-ahead right after a recovery
            careful = max(careful-1, 0)
            dc = along(min(s_prog+look, cum[-1]))-x[:2]
            if wall_follow > 0:              # near a wall, blend the carrot direction into the wall tangent
                ph, gr = f.eval(np.r_[x[:2], x[2] % C.TWO_PI], grad=True)
                hn, gt_ = float(ph[0]-level), gr[0, :2]
                if hn < wall_follow and np.linalg.norm(gt_) > 1e-9:
                    tan = np.r_[-gt_[1], gt_[0]]/np.linalg.norm(gt_)
                    tan = tan if tan@dc >= 0 else -tan
                    wgt = np.clip(1-max(hn, 0.)/wall_follow, 0, 1)
                    dc = (1-wgt)*dc/max(np.linalg.norm(dc), 1e-12)+wgt*tan
            al = C.wrap(np.arctan2(dc[1], dc[0])-x[2])
            speed = vmax if s_prog+lookahead < cum[-1] else min(vmax, 1.5*np.linalg.norm(d))
            nominal = np.r_[speed*max(np.cos(al), 0.), np.clip(2*al, -wmax, wmax)]
            # recovery from the filter's undesired equilibria: if the filtered input has stayed
            # (almost) zero for 0.5 s, back up slowly for 1 s while turning toward the carrot
            if controls and np.linalg.norm(controls[-1]) < .02:
                stall += 1
            else:
                stall = 0
            if recovery and stall >= 50:
                backup, stall = 100, 0
                recoveries += 1
            if backup > 0:                       # back up while aligning with the route tangent
                backup -= 1
                tg = along(min(s_prog+.05, cum[-1]))-along(max(s_prog-.05, 0.))
                gth = f.eval(np.r_[x[:2], x[2] % C.TWO_PI], grad=True)[1][0, 2]
                turn = (.8*wmax*np.sign(gth) if abs(gth) > 1e-3 else    # rotate away from the wall
                        np.clip(2*C.wrap(np.arctan2(tg[1], tg[0])-x[2]), -wmax, wmax))
                nominal = np.r_[-.25, turn]
                if backup == 0:
                    careful = 200
        else:
            nominal = np.r_[v, np.clip(2*C.wrap(goal[2]-x[2]), -wmax, wmax)]
        if summed_fn is not None:            # summed-field barriers: Bernstein coefficient rows, epigraph t >= |u|
            import summed as SM
            res = summed_fn(x)
            A_, T_, C_, E_ = ((res[0], res[1], res[2], res[5]) if res is not None else
                              (np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0), np.zeros((3, 3))))
            u, slack, z_sum = SM.solve(nominal, A_, T_, C_, E_, Gbox, hbox, z_sum)
            timings.append(time.perf_counter()-t0)
            n_active.append(len(C_))
            n_boxes.append(res[3] if res is not None else 0)
            slack_steps += int(slack > 0)
            barriers.append(res[4] if res is not None else np.inf)
            residuals.append(float((A_@u+T_@np.abs(E_@u)+C_).min()) if len(C_) else np.inf)
            nominal_residuals.append(float((A_@nominal+T_@np.abs(E_@nominal)+C_).min()) if len(C_) else np.inf)
            nxt = x+dt*u
            gn = gap_function(nxt, robot_gt, walls)
            audit(x, nxt, gaps[-1], gn)
            controls.append(u); traj.append(nxt.copy()); gaps.append(gn)
            x = nxt
            if wp == len(waypoints)-1 and np.linalg.norm(x[:2]-goal[:2]) < .05 and abs(C.wrap(x[2]-goal[2])) < .05:
                reached = (step+1)*dt
                break
            continue
        if rows_fn is None:                  # tabulated field: one barrier
            assert np.max(np.abs(x[:2])) < domain-C.ACT_DELTA
            phi, grad = f.eval(np.r_[x[:2], x[2] % C.TWO_PI], grad=True)
            hv, rows = np.array([float(phi[0]-level)]), grad[:1]
        else:                                # surface cover: one barrier per active box corner
            rows, hv = rows_fn(x)
        if dyn == 'uni':                     # eq. (hdot-uni) for a static field: [g_t . e(theta), g_theta]
            rows = np.c_[rows[:, :2] @ np.r_[np.cos(x[2]), np.sin(x[2])], rows[:, 2]]
        near = hv < .05
        if near.any():
            min_row_near = min(min_row_near, float(np.linalg.norm(rows[near], axis=1).min()))
        u, slack = solve_ldp_qp(nominal, np.vstack([rows,Gbox]), np.r_[-gamma*hv,hbox], len(hv),
                                 controls[-1] if controls else None)
        timings.append(time.perf_counter()-t0)
        n_active.append(len(hv))
        slack_steps += int(slack > 0)
        barriers.append(float(hv.min()) if len(hv) else np.inf)
        residuals.append(float((rows@u+gamma*hv).min()) if len(hv) else np.inf)
        nominal_residuals.append(float((rows@nominal+gamma*hv).min()) if len(hv) else np.inf)
        nxt = x+dt*(np.r_[u[0]*np.cos(x[2]), u[0]*np.sin(x[2]), u[1]] if dyn == 'uni' else u)
        gn = gap_function(nxt, robot_gt, walls)
        audit(x, nxt, gaps[-1], gn)
        controls.append(u)
        traj.append(nxt.copy())
        gaps.append(gn)
        x = nxt
        if wp == len(waypoints)-1 and np.linalg.norm(x[:2]-goal[:2]) < (.30 if dyn == 'uni' else .05) \
                and (dyn == 'uni' or abs(C.wrap(x[2]-goal[2])) < .05):
            reached = (step+1)*dt
            break
    if summed_fn is not None:
        res = summed_fn(x)
        final_phi = res[4] if res is not None else np.inf
    elif rows_fn is None:
        final_phi = float(f.eval(np.r_[x[:2], x[2] % C.TWO_PI])[0])-level
    else:
        hv_end = rows_fn(x)[1]
        final_phi = float(hv_end.min()) if len(hv_end) else np.inf
    arrays = dict(poses=np.array(traj), controls=np.array(controls), gaps=np.array(gaps),
        barriers=np.r_[barriers,final_phi], timings=np.array(timings), waypoints=waypoints,
        nominal_residuals=np.array(nominal_residuals), n_active=np.array(n_active), n_boxes=np.array(n_boxes))
    result = dict(dt=dt, dynamics=dyn, recoveries=recoveries, min_row_norm_near_boundary=min_row_near,
        speed_bound=vmax, angular_speed_bound=wmax, gamma=gamma,
        reached_s=reached, steps=len(controls), min_physical_gap_m=float(np.min(gaps)),
        verified_interval_gap_lower_m=float(verified_lower), unresolved_intervals=unresolved,
        collision_intervals=collisions, min_h_m=float(min(min(barriers),final_phi)),
        min_constraint_residual=float(min(residuals)), slack_steps=slack_steps,
        barrier_intervention_steps=int(np.sum(np.array(nominal_residuals)<-1e-8)),
        barriers_max=int(max(n_active)) if n_active else 0,
        controller_median_ms=float(np.median(timings)*1000), controller_p95_ms=float(np.percentile(timings,95)*1000),
        path_length_m=float(np.linalg.norm(np.diff(arrays['poses'][:,:2],axis=0),axis=1).sum()),
        final_position_error_m=float(np.linalg.norm(x[:2]-waypoints[-1,:2])))
    print(json.dumps(result,indent=2), flush=True)
    (output_dir/'last_result.json').write_text(json.dumps(result,indent=2,default=float))
    # discrete-time dips of h are reported, not rejected; safety is judged on the physical geometry
    if reached is None or collisions or unresolved or slack_steps or result['verified_interval_gap_lower_m'] <= 0:
        np.savez_compressed(output_dir/'diagnostic.npz', **arrays)
        raise RuntimeError('Escape/clearance/CBF validation failed; saved diagnostic, no paper export.')
    return arrays, result


def draw_walls(ax, walls):
    # Draw fills first, then ONLY exposed cubic boundaries. Internal beam seams
    # are not boundaries of a cross/T and must not appear as outlines.
    for wall in walls:
        ax.add_patch(Polygon(wall.outline, closed=True, facecolor='#92979e', edgecolor='none', zorder=3))
    for i, wall in enumerate(walls):
        points = np.vstack([wall.dense,wall.dense[:1]])
        exposed = np.ones(len(points),dtype=bool)
        for j,other in enumerate(walls):
            if i == j or GROUP_IDS[i] != GROUP_IDS[j]:
                continue
            ab = other.b-other.a
            t = np.clip(((points-other.a)@ab)/(ab@ab),0,1)
            exposed &= np.linalg.norm(points-other.a-t[:,None]*ab,axis=1) >= other.radius-1e-8
        xy = np.ma.array(points, mask=np.broadcast_to(~exposed[:,None],points.shape))
        ax.plot(xy[:,0],xy[:,1],color='#101820',lw=1.6,zorder=4)


def plot(arrays, robot_gt, walls, robot, preview=False):
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'pdf.fonttype':42})
    fig, ax = plt.subplots(figsize=(7.4,6.9))
    draw_walls(ax,walls)
    poses = arrays['poses']
    ax.plot(poses[:,0], poses[:,1], color='#286dcc', lw=1.8, ls='--' if preview else '-', zorder=2)
    # Uniform arc-length spacing avoids piles of silhouettes at slow corners.
    arc = np.r_[0,np.linalg.norm(np.diff(poses[:,:2],axis=0),axis=1).cumsum()]
    ids = np.searchsorted(arc, np.arange(1.1, arc[-1]-.9, 1.65)) if not preview else []
    for idx in ids:
        points = robot_gt @ C.rot(poses[idx,2]).T + poses[idx,:2]
        ax.add_patch(Polygon(points, facecolor='#dce1e8', edgecolor='none', zorder=5))
        fit = robot.world(poses[idx]); fit = np.vstack([fit,fit[:1]])
        ax.plot(fit[:,0],fit[:,1],color='#124395',lw=1.5,zorder=6)
    for idx in (0,len(poses)-1):
        points = robot_gt @ C.rot(poses[idx,2]).T + poses[idx,:2]
        ax.add_patch(Polygon(points, facecolor='#759bd3', edgecolor='none', zorder=5))
        fit = robot.world(poses[idx]); fit = np.vstack([fit,fit[:1]])
        ax.plot(fit[:,0],fit[:,1],color='#082d71',lw=1.9,zorder=6)
    for distance in np.arange(3.,arc[-1],5.):
        idx = np.searchsorted(arc,distance)
        if idx+20 < len(poses):
            ax.annotate('', xy=poses[idx+20,:2], xytext=poses[idx,:2],
                arrowprops=dict(arrowstyle='-|>',color='#174f9d',lw=1.5), zorder=6)
    ax.text(-4.1,-4.59,'Start',ha='center',color='#163d79',weight='bold')
    ax.text(5.3,4.45,'Exit',ha='center',color='#163d79',weight='bold')
    ax.set(xlim=(-5.65,6.0),ylim=(-5.45,5.45),xlabel='$x$ [m]',ylabel='$y$ [m]',aspect='equal')
    ax.set_xticks(np.arange(-4,7,2)); ax.set_yticks(np.arange(-4,5,2))
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout(pad=.5)
    if preview:
        ax.set_title('Layout preview / prescribed waypoints (not a simulated trajectory)')
        fig.savefig(OUT/'layout_preview.png',dpi=170,bbox_inches='tight')
        plt.close(fig)
        return
    for path in (FIG/'maze_escape.png',OUT/'maze_escape.png'):
        fig.savefig(path,dpi=240,bbox_inches='tight',facecolor='white')
    fig.savefig(FIG/'maze_escape.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    fig, axes = plt.subplots(2,1,figsize=(7,4.6),sharex=True)
    t = np.arange(len(poses))*.01
    axes[0].plot(t, arrays['gaps'].min(axis=1)*1000,color='#286dcc')
    axes[0].set_ylabel('Physical gap [mm]')
    axes[1].plot(t,arrays['barriers']*1000,color='#286dcc')
    axes[1].axhline(0,color='0.4',ls='--',lw=.8)
    axes[1].set(ylabel='$h$ [mm]',xlabel='Time [s]')
    fig.tight_layout(); fig.savefig(OUT/'clearance.png',dpi=180); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--rebuild',action='store_true')
    ap.add_argument('--layout-only',action='store_true'); args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    data = np.load(ROOT/'cache'/'shapes5.npz')
    robot, robot_gt = C.SmoothShape(data['ctrl'][0]), data['GT'][0]
    walls = [Capsule(*w) for w in WALLS]
    separation = validate_geometry(walls)
    if args.layout_only:
        plot(dict(poses=WAYPOINTS),robot_gt,walls,robot,preview=True)
        return
    f, field_meta = make_field(robot,walls,args.rebuild)
    arrays, result = simulate(f,field_meta['level'],robot_gt,walls)
    result.update(wall_count=len(WALL_GROUPS), capsule_beam_count=len(walls), walls=WALL_GROUPS,
        min_wall_separation_m=separation,
        field=field_meta, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        navigation='Prescribed waypoint sequence; CBF is the safety filter, not a route planner.',
        geometry='Eight disconnected walls: unions of exact capsules, with cubic enclosing beam boundaries. Cached robot GT polygon.',
        verification='Float64 contact/band/regularity checks and interval Lipschitz clearance audit; no outward rounding.')
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(OUT/'trajectory.npz',**arrays)
    plot(arrays,robot_gt,walls,robot)
    tex = ('%% Generated by maze_escape.py; separate from the obstacle-count study.\n'
        '\\newcommand{\\mazeTime}{%.1f}\n'
        '\\newcommand{\\mazeGap}{%.1f}\n'
        '\\newcommand{\\mazeLower}{%.1f}\n'
        '\\newcommand{\\mazeLevel}{%.1f}\n'
        '\\newcommand{\\mazeMedian}{%.2f}\n'
        '\\newcommand{\\mazeInterventions}{%d}\n') % (
            result['reached_s'],result['min_physical_gap_m']*1000,
            result['verified_interval_gap_lower_m']*1000,field_meta['level']*1000,
            result['controller_median_ms'],result['barrier_intervention_steps'])
    (ROOT/'tro'/'maze_results.tex').write_text(tex)


if __name__ == '__main__':
    from reference_maze_escape import main as reference_main
    reference_main()
