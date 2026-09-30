"""
Star tube: the Franka Panda inserts its hand into a star-shaped tube (6.7 cm long): a 5-pointed star hole
surrounded by walls of constant 5-cm thickness (parallel edges). Ground truth: outer star minus hole star,
extruded; its exact SDF is analytic. All eight links are constrained
by the surface-cover barriers (paper, Sec. VI); the nominal controller is Cartesian damped-least-squares IK
that drives the hand straight through the center of the hole.

    python prototype_3d/dolphin3d.py [ours spheres points kkt]
"""
import json
import os
import sys
import time

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np

import franka3d as F
from proto3d import Spline3
from sdf_cbf_utils import solve_ldp_qp

HOOP_C = np.array([.64, 0., .55])         # center; the cross-section lies in the y-z plane, axis along x
# star tube with walls of constant thickness: the hole is a 5-pointed star, the outer boundary is the same star
# with every edge moved outward by WALL (parallel edges), extruded along x
R_HOUT, R_HIN, N_PTS, HALF_T, WALL = .24, .15, 5, .10 / 3, .05   # hole star, half depth (6.7-cm tube), wall
START_EE, GOAL_EE = np.array([.36, -.22, .40]), np.array([.70, 0., .55])   # goal 16 cm inside the tube
TCP = np.array([0., 0., .1034])
Q_START = np.array([0., -.35, 0., -2.55, 0., 2.2, .785])   # refined by IK to reach START_EE
STEPS, KP = 1200, 1.5


def star_polygon(r_out, r_in):
    a = np.pi / 2 + np.arange(2 * N_PTS) * np.pi / N_PTS
    r = np.where(np.arange(2 * N_PTS) % 2 == 0, r_out, r_in)
    return np.c_[r * np.cos(a), r * np.sin(a)]


def offset_polygon(P, d):
    """Polygon whose edges are those of the counterclockwise polygon P moved outward by d (parallel edges;
    vertices at the intersections of consecutive moved edges)."""
    Q = np.roll(P, -1, axis=0)
    e = Q - P
    n = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1, keepdims=True)     # outward normals (CCW)
    A, E = P + d * n, e                                                            # moved edge k: A_k + s E_k
    out = []
    for k in range(len(P)):
        j = k - 1                                                                  # edges j and k meet at P_k
        M = np.c_[E[j], -E[k]]
        s_, _ = np.linalg.solve(M, A[k] - A[j])
        out.append(A[j] + s_ * E[j])
    return np.array(out)


HOLE = star_polygon(R_HOUT, R_HIN)
OUTER = offset_polygon(HOLE, WALL)
R_OUT = float(np.linalg.norm(OUTER, axis=1).max())                 # bounding radius of the cross-section
EDGES = np.concatenate([np.stack([P, np.roll(P, -1, axis=0)], 1) for P in (OUTER, HOLE)])


def sdf2(Q):
    """Exact signed distance to the planar region outer star minus hole star (even-odd rule over all edges)."""
    A, B = EDGES[:, 0], EDGES[:, 1]
    d = np.full(len(Q), np.inf)
    inside = np.zeros(len(Q), bool)
    for s in range(0, len(A), 64):
        a, b = A[s:s + 64], B[s:s + 64]
        ab = b - a
        t = np.clip(np.einsum('pkd,kd->pk', Q[:, None] - a, ab) / np.einsum('kd,kd->k', ab, ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(Q[:, None] - a - t[..., None] * ab, axis=2).min(axis=1))
        cond = (a[None, :, 1] > Q[:, None, 1]) != (b[None, :, 1] > Q[:, None, 1])
        with np.errstate(divide='ignore', invalid='ignore'):
            xint = a[None, :, 0] + (Q[:, None, 1] - a[None, :, 1]) * ab[None, :, 0] / ab[None, :, 1]
        inside ^= (cond & (Q[:, None, 0] < xint)).sum(axis=1) % 2 == 1
    return np.where(inside, -d, d)


def sdf2_edges(Q, edges):
    """Exact signed distance to the planar region bounded by the closed polygons whose edges are `edges`
    (even-odd rule)."""
    A, B = edges[:, 0], edges[:, 1]
    d = np.full(len(Q), np.inf)
    inside = np.zeros(len(Q), bool)
    for s in range(0, len(A), 64):
        a, b = A[s:s + 64], B[s:s + 64]
        ab = b - a
        t = np.clip(np.einsum('pkd,kd->pk', Q[:, None] - a, ab) / np.einsum('kd,kd->k', ab, ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(Q[:, None] - a - t[..., None] * ab, axis=2).min(axis=1))
        cond = (a[None, :, 1] > Q[:, None, 1]) != (b[None, :, 1] > Q[:, None, 1])
        with np.errstate(divide='ignore', invalid='ignore'):
            xint = a[None, :, 0] + (Q[:, None, 1] - a[None, :, 1]) * ab[None, :, 0] / ab[None, :, 1]
        inside ^= (cond & (Q[:, None, 0] < xint)).sum(axis=1) % 2 == 1
    return np.where(inside, -d, d)


# rounded variant: every corner of both stars filleted with circular arcs of radius R_FIL, and the edges where
# the faces meet the walls rounded with R_EDGE (rounded extrusion). The cross-section boundary is a C^1 curve of
# straight segments and arcs, so its signed distance is analytic: the distance to the nearest piece, signed by
# the side of that piece (exact because the curve has a tangent everywhere).
R_FIL, R_EDGE = .015, .01


def _filleted(P, r):
    """Pieces of the CCW polygon P with every corner filleted with radius r: segments (a, b) and arcs
    (center, radius, start angle, signed sweep, convex)."""
    n = len(P)
    arcs, tang = [], []
    for k in range(n):
        u, v, w = P[k - 1], P[k], P[(k + 1) % n]
        d1, d2 = (u - v) / np.linalg.norm(u - v), (w - v) / np.linalg.norm(w - v)
        convex = (v - u)[0] * (w - v)[1] - (v - u)[1] * (w - v)[0] > 0
        half = np.arccos(np.clip(d1 @ d2, -1, 1)) / 2
        t = r / np.tan(half)
        c = v + (d1 + d2) / np.linalg.norm(d1 + d2) * r / np.sin(half)
        p1, p2 = v + d1 * t, v + d2 * t
        a1, a2 = np.arctan2(*(p1 - c)[::-1]), np.arctan2(*(p2 - c)[::-1])
        arcs.append((c, r, a1, (a2 - a1 + np.pi) % (2 * np.pi) - np.pi, convex))
        tang.append((p1, p2))
    segs = [(tang[k][1], tang[(k + 1) % n][0]) for k in range(n)]
    return segs, arcs


PIECES_R = [_filleted(P, R_FIL) for P in (OUTER, HOLE)]


def _sdf_curve(Q, segs, arcs):
    """Signed distance to the region enclosed by one C^1 curve (negative inside)."""
    best, sign = np.full(len(Q), np.inf), np.ones(len(Q))
    for a, b in segs:
        ab = b - a
        t = np.clip((Q - a) @ ab / (ab @ ab), 0, 1)
        d = np.linalg.norm(Q - a - t[:, None] * ab, axis=1)
        s = np.where(ab[0] * (Q[:, 1] - a[1]) - ab[1] * (Q[:, 0] - a[0]) > 0, -1., 1.)   # left of CCW: inside
        m = d < best
        best[m], sign[m] = d[m], s[m]
    for c, r, a0, sw, convex in arcs:
        v = Q - c
        rho = np.linalg.norm(v, axis=1)
        ang = (np.arctan2(v[:, 1], v[:, 0]) - a0) * np.sign(sw) % (2 * np.pi)
        on = ang <= abs(sw)
        ends = [c + r * np.array([np.cos(a0), np.sin(a0)]), c + r * np.array([np.cos(a0 + sw), np.sin(a0 + sw)])]
        dend = np.minimum(*[np.linalg.norm(Q - e, axis=1) for e in ends])
        d = np.where(on, np.abs(rho - r), dend)
        s = np.where((rho < r) == convex, -1., 1.)          # convex fillet: inside toward the center
        m = d < best
        best[m], sign[m] = d[m], s[m]
    return sign * best


def sdf2_rounded(Q):
    """Exact signed distance to the rounded cross-section (outer curve minus hole curve)."""
    so = _sdf_curve(Q, *PIECES_R[0])
    sh = _sdf_curve(Q, *PIECES_R[1])
    inside = (so < 0) & (sh > 0)
    return np.where(inside, -1., 1.) * np.minimum(np.abs(so), np.abs(sh))


def hoop_sdf_rounded(P):
    """3D SDF of the rounded star tube (world frame): the rounded cross-section, extruded with edges of radius
    R_EDGE (exact outside the body, 1-Lipschitz everywhere)."""
    q = np.atleast_2d(P) - HOOP_C
    w = np.c_[sdf2_rounded(q[:, 1:3]) + R_EDGE, np.abs(q[:, 0]) - HALF_T + R_EDGE]
    return np.minimum(w.max(axis=1), 0) + np.linalg.norm(np.maximum(w, 0), axis=1) - R_EDGE


def hoop_sdf(P):
    """Exact 3D SDF of the extruded region (world frame)."""
    q = np.atleast_2d(P) - HOOP_C
    out = np.empty(len(q))
    for s in range(0, len(q), 20000):
        qs = q[s:s + 20000]
        w = np.c_[sdf2(qs[:, 1:3]), np.abs(qs[:, 0]) - HALF_T]
        out[s:s + 20000] = np.minimum(w.max(axis=1), 0) + np.linalg.norm(np.maximum(w, 0), axis=1)
    return out


def build_hoop():
    cache = os.path.join(F.HERE, 'cache_startube.npz')
    lo, hi = HOOP_C - np.array([HALF_T, R_OUT, R_OUT]), HOOP_C + np.array([HALF_T, R_OUT, R_OUT])
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)['d'].item()
        f = F.spline_from(d['f'])
        l = d['l']
    else:
        t0 = time.perf_counter()
        f = Spline3(lo - .15, hi + .15, .015).fit(hoop_sdf, .0075)
        l = F.tight_level_implicit(f, F.hessian_majorant(f), hoop_sdf, lo, hi, h0=.006)
        np.savez_compressed(cache, d=np.array(dict(f=F.spline_state(f), l=l), dtype=object))
        print(f'hoop field + tight level {1e3 * l:.2f} mm in {time.perf_counter() - t0:.0f} s', flush=True)
    links, _, _ = F.build()
    reach = max(L['crad'].max() for L in links)
    Gn, Mn = F.neighborhood_bounds(f, reach)
    hoop = dict(f=f, M=F.hessian_majorant(f), Gn=Gn, Mn=Mn, l=l, blo=lo, bhi=hi, dlo=f.lo, dhi=f.lo + f.K * f.h,
                sdf=hoop_sdf)
    return links, {'hoop': hoop}


def real_gap(q, links, obst, per_link=False):
    """Certified lower bound on the distance between the link meshes and the hoop (exact SDF at the centroids of
    the subdivided triangles minus their radii)."""
    import trimesh
    T, _, _ = F.fk(q)
    o = obst['hoop']
    gaps = []
    for L in links:
        if 'gt8' not in L:                 # triangles subdivided to 8-mm edges (bound stays certified)
            Vs, Fs = trimesh.remesh.subdivide_to_size(L['V'], L['F'], max_edge=.008)
            Tr = Vs[Fs]
            L['gt8'] = Tr.mean(axis=1)
            L['gt8_r'] = np.linalg.norm(Tr - L['gt8'][:, None], axis=2).max(axis=1)
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        C = L['PC'] @ R.T + p              # whole link far: skip (spline lower bound over its box cover)
        inb = np.all((C > o['dlo']) & (C < o['dhi']), axis=1)
        if not inb.any() or o['f'].eval(C[inb]).min() - L['r'] > GAP_CAP + .02:
            gaps.append(GAP_CAP)
            continue
        W = L['gt8'] @ R.T + p
        m = np.flatnonzero(np.all((W > o['blo'] - .05) & (W < o['bhi'] + .05), axis=1))
        if len(m):                         # the spline field (fit error of a few mm) preselects the candidates
            m = m[o['f'].eval(W[m]) < GAP_CAP + .01]
        gaps.append(min(float((o['sdf'](W[m]) - L['gt8_r'][m]).min()), GAP_CAP) if len(m) else GAP_CAP)
    return gaps if per_link else min(gaps)


GAP_CAP = .03                                          # distances above 3 cm are reported as 3 cm


def ee_pose(q):
    T, Z, O = F.fk(q)
    R8 = T[8][:3, :3]
    return R8 @ TCP + T[8][:3, 3], R8, Z, O


def log_so3(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    th = np.arccos(c)
    if th < 1e-9:
        return np.zeros(3)
    return th / (2 * np.sin(th)) * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])


R_GOAL = np.array([[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]])     # gripper axis (z) along +x, fingers in y


def nominal(q, goal):
    p, R8, Z, O = ee_pose(q)
    J = np.vstack([np.cross(Z, p - O).T, Z.T])
    e = np.r_[np.clip(KP * (goal - p), -.4, .4), np.clip(1.5 * log_so3(R_GOAL @ R8.T), -1., 1.)]
    Jp = J.T @ np.linalg.inv(J @ J.T + 1e-3 * np.eye(6))
    return np.clip(Jp @ e + (np.eye(7) - Jp @ J) @ (.3 * (Q_START - q)), -F.QD_MAX, F.QD_MAX), p


def start_config():
    global Q_START
    q = Q_START.copy()
    for _ in range(4000):                               # IK to the start pose within the joint limits
        u, p = nominal(q, START_EE)
        q = np.clip(q + .01 * u, F.Q_MIN + .05, F.Q_MAX - .05)
    Q_START = q.copy()                                  # null-space posture of the run
    return q


def simulate(q0, links, obst, method, prim=None, steps=STEPS, record=True):
    import franka_baselines as Bl
    q = np.array(q0, float)
    log = dict(t=[], rows=[], h=[], gap=[], slack=0, reached=None, q=[], ee=[])
    u, state = None, {}
    lo = lambda q: np.clip(-2. * (q - F.Q_MIN), -F.QD_MAX, 0.)      # joint limits as velocity bounds
    hi = lambda q: np.clip(2. * (F.Q_MAX - q), 0., F.QD_MAX)
    for k in range(steps):
        u_nom, p = nominal(q, GOAL_EE)
        t0 = time.perf_counter()
        if method == 'ours':
            G, h = F.barrier_rows(q, links, obst)
        elif method == 'kkt':
            G, h = Bl.kkt_rows(q, links, obst, state)
        elif method == 'nominal':
            G, h = np.zeros((0, 7)), np.zeros(0)
        else:
            G, h = Bl.rows(q, links, obst, prim)
        Gall = np.vstack([G, np.eye(7), -np.eye(7)])
        hall = np.r_[-F.GAMMA * h, lo(q), -hi(q)]
        u, s = solve_ldp_qp(u_nom, Gall, hall, len(h), u)
        log['t'].append(time.perf_counter() - t0)
        log['slack'] += s > 0
        log['rows'].append(len(h)); log['h'].append(float(h.min()) if len(h) else np.inf)
        log['gap'].append(real_gap(q, links, obst))
        if record:
            log['q'].append(q.copy()); log['ee'].append(p.copy())
        q = q + F.DT * u
        if np.linalg.norm(GOAL_EE - p) < .01 and log['reached'] is None:
            log['reached'] = (k + 1) * F.DT
    return log


def main():
    methods = sys.argv[1:] or ['nominal', 'ours', 'spheres', 'points', 'kkt']
    links, obst = build_hoop()
    print(f"hoop level {1e3 * obst['hoop']['l']:.2f} mm", flush=True)
    q0 = start_config()
    p0 = ee_pose(q0)[0]
    print('start EE', p0.round(3), 'gap', 1e3 * real_gap(q0, links, obst), flush=True)
    import franka_baselines as Bl
    out = {}
    path = os.path.join(F.HERE, 'startube_results.json')
    if os.path.exists(path):
        out = json.load(open(path))
    for m in methods:
        prim = Bl.primitives(links, m) if m in ('spheres', 'points') else None
        L = simulate(q0, links, obst, m, prim)
        ee = np.array(L['ee'])
        r = dict(method=m, reached=L['reached'], min_gap_mm=1e3 * min(L['gap']), final_ee_err_mm=1e3 * float(
            np.linalg.norm(GOAL_EE - ee[-1])), inserted_mm=1e3 * float(ee[:, 0].max() - (HOOP_C[0] - HALF_T)),
                 slack_steps=int(L['slack']), min_h_mm=1e3 * min(L['h']),
                 t_med_ms=1e3 * float(np.median(L['t'])), t_max_ms=1e3 * float(np.max(L['t'])),
                 rows_max=int(max(L['rows'])))
        out[m] = r
        print(m, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        np.savez_compressed(os.path.join(F.HERE, f'startube_traj_{m}.npz'), q=np.array(L['q']), ee=ee,
                            gap=np.array(L['gap']), h=np.array(L['h']))
        json.dump(out, open(path, 'w'), indent=1)


if __name__ == '__main__':
    main()
