"""
Franka Panda (meshes and kinematics of github.com/idiap/RDF) among non-convex obstacles: one set of
barrier constraints for EVERY (link, obstacle) pair, no closest points, no optimization online.

Geometry. Link k's real surface is its visual mesh; each obstacle's real shape is a union of rounded
boxes. Every body has a smooth tensor cubic B-spline SDF phi, and a certified level l that encloses the
whole real surface, found by sampling: grid points within delta = s sqrt(3)/2 of the real surface
cover it, and for x within delta of a grid point g
    phi(x) <= phi(g) + |grad phi(g)| delta + M(g) delta^2 / 2 ,
with M a certified Hessian majorant (Bernstein coefficients of the second derivatives). So
    real surface  subset  {phi <= l}.
Offline, per link: the continuous level set {phi_A = l_A} is covered by small boxes whose Bernstein
coefficients of phi_A - l_A change sign (convex-hull property), halved to <= 1 cm.
Online, per (link, obstacle) pair and box (center c, half-diagonal r, corners v), one smooth barrier per
corner
    h_v(q) = phi_O(c) + grad phi_O(c)^T (v - c) - M_O(c) r^2 / 2 - l_O   <=  min_{x in box} phi_O(x) - l_O,
so h_v >= 0 for all v keeps {phi_A = l_A}, hence the real link, out of {phi_O < l_O}, which contains the
real obstacle.

Evaluation: certified lower bound on the real distance between every link mesh and every obstacle
(exact obstacle SDF at the centroids of the subdivided mesh triangles minus the triangle radii).

    python prototype_3d/franka3d.py [n_trials]
"""
import json
import os
import sys
import time

os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np
import trimesh
from scipy.cluster.vq import kmeans2
from scipy.ndimage import maximum_filter
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from proto3d import Spline3, sdf_boxes, surface_boxes   # noqa: E402
from sdf_cbf_utils import solve_ldp_qp                   # noqa: E402

MESH = os.path.join(HERE, 'RDF', 'panda_layer', 'meshes')
CACHE = os.path.join(HERE, 'cache_franka.npz')

# modified DH (A, alpha, D) of panda_layer.py; frame 8 (hand) is fixed to frame 7 with theta = -pi/4
DH = [(0., 0., .333), (0., -np.pi / 2, 0.), (0., np.pi / 2, .316), (.0825, np.pi / 2, 0.),
      (-.0825, -np.pi / 2, .384), (0., np.pi / 2, 0.), (.088, np.pi / 2, 0.), (0., 0., .107)]
Q_MIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
Q_MAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
DOF = 7
LINK_IDS = [1, 2, 3, 4, 5, 6, 7, 8]            # link0 is the fixed base

OBSTACLES = {'arch': [((.55, -.25, .38), (.04, .04, .38)),                # gate
                      ((.55, .25, .38), (.04, .04, .38)),
                      ((.55, 0., .76), (.04, .29, .04))],
             'hook': [((0., .62, .30), (.04, .04, .30)),                  # inverted J
                      ((.20, .62, .56), (.24, .04, .04)),
                      ((.40, .62, .42), (.035, .035, .16))],
             'rack': [((0., -.62, .20), (.30, .15, .02)),                 # plate + thin rod on a leg
                      ((0., -.62, .48), (.30, .012, .012)),
                      ((.28, -.62, .34), (.015, .015, .14))]}
GAMMA, DT, QD_MAX, KQ, ACT, BOX_SIDE, S_LEVEL = 5., .01, 1., 1.5, .03, .01, .004
CLUSTER = 48                                   # boxes per cluster for the coarse pruning (radius <= 3 cm)
CORNERS = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float) / 2


def dh(A, al, D, th):
    ct, st, ca, sa = np.cos(th), np.sin(th), np.cos(al), np.sin(al)
    return np.array([[ct, -st, 0., A], [st * ca, ct * ca, -sa, -sa * D],
                     [st * sa, ct * sa, ca, ca * D], [0., 0., 0., 1.]])


def fk(q):
    """Frames T[0..8] (T[0] = base) and joint axes / origins (world) of joints 1..7."""
    T = [np.eye(4)]
    for i, (A, al, D) in enumerate(DH):
        T.append(T[-1] @ dh(A, al, D, q[i] if i < DOF else -np.pi / 4))
    Z = np.array([T[i + 1][:3, 2] for i in range(DOF)])
    O = np.array([T[i + 1][:3, 3] for i in range(DOF)])
    return T, Z, O


# ---------------------------------------------------------------------------------------------
# certified quantities of a spline SDF
# ---------------------------------------------------------------------------------------------
def hessian_majorant(f):
    """C^2 spline M(c) >= certified Hessian bound on every cell within one cell of c's cell: control point
    m takes the max over cells m-4..m+1 per axis; cubic B-spline weights are a partition of unity."""
    W = np.pad(f.hessian_bound(cellwise=True), 4, mode='edge')
    for a in range(3):
        W = np.lib.stride_tricks.sliding_window_view(W, 6, axis=a).max(axis=-1)
    g = Spline3.__new__(Spline3)
    g.lo, g.h, g.K, g.n, g.W = f.lo, f.h, f.K, f.n, W
    return g


def neighborhood_bounds(f, reach):
    """Per-cell certified bounds of |grad phi| and of the Hessian norm over all cells within `reach` of
    the cell (for pruning a whole ball of radius <= reach at once)."""
    C = f.cell_coeffs()
    G = np.sqrt(sum((3 * np.abs(np.diff(C, axis=3 + a)).max(axis=(3, 4, 5)) / f.h) ** 2 for a in range(3)))
    k = int(np.ceil(reach / f.h)) + 1                   # cells a ball of radius reach can touch
    # the majorant M(c) at a box center reads cells within 4 more of c's cell
    return (maximum_filter(G, size=2 * k + 1), maximum_filter(f.hessian_bound(cellwise=True), size=2 * (k + 4) + 1))


def clusters(PC, r):
    n = max(1, len(PC) // CLUSTER)
    cc, lab = kmeans2(PC, n, seed=0, minit='++')
    groups, cen, rad = [], [], []
    for j in range(n):
        m = np.flatnonzero(lab == j)
        if len(m):
            groups.append(m); cen.append(cc[j]); rad.append(np.linalg.norm(PC[m] - cc[j], axis=1).max() + r)
    return groups, np.array(cen), np.array(rad)


def point_triangle_distance(P, A, B, C):
    """Exact distance from points P to triangles (A, B, C), row-wise (closest-point regions, Ericson)."""
    dot = lambda x, y: np.einsum('nd,nd->n', x, y)
    ab, ac, ap = B - A, C - A, P - A
    d1, d2 = dot(ab, ap), dot(ac, ap)
    bp = P - B
    d3, d4 = dot(ab, bp), dot(ac, bp)
    cp = P - C
    d5, d6 = dot(ab, cp), dot(ac, cp)
    vc, vb, va = d1 * d4 - d3 * d2, d5 * d2 - d1 * d6, d3 * d6 - d5 * d4
    with np.errstate(divide='ignore', invalid='ignore'):
        e_ab = A + (d1 / (d1 - d3))[:, None] * ab
        e_ac = A + (d2 / (d2 - d6))[:, None] * ac
        e_bc = B + ((d4 - d3) / ((d4 - d3) + (d5 - d6)))[:, None] * (C - B)
        den = 1. / (va + vb + vc)
        face = A + (vb * den)[:, None] * ab + (vc * den)[:, None] * ac
    conds = [(d1 <= 0) & (d2 <= 0), (d3 >= 0) & (d4 <= d3), (vc <= 0) & (d1 >= 0) & (d3 <= 0),
             (d6 >= 0) & (d5 <= d6), (vb <= 0) & (d2 >= 0) & (d6 <= 0),
             (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)]
    X = face.copy()
    for cnd, val in reversed(list(zip(conds, [A, B, e_ab, C, e_ac, e_bc]))):
        X = np.where(cnd[:, None], val, X)
    return np.linalg.norm(P - X, axis=1)


def mesh_distance(V, F, P, upto):
    """Exact distance from P to the triangle mesh (V, F); points farther than `upto` get a value > upto."""
    T = V[F]
    cen = T.mean(axis=1)
    Rt = np.linalg.norm(T - cen[:, None], axis=2).max(axis=1)
    tc = cKDTree(cen)
    dc = tc.query(P)[0]
    out = np.full(len(P), np.inf)
    cand = np.flatnonzero(dc - Rt.max() <= upto)             # exact distance >= dc - Rmax
    ub = cKDTree(V).query(P[cand])[0]                        # exact distance <= nearest vertex distance
    lists = tc.query_ball_point(P[cand], ub + Rt.max())      # the closest triangle is in this ball
    cnt = np.array([len(x) for x in lists])
    pi = np.repeat(np.arange(len(cand)), cnt)
    ti = np.concatenate([np.asarray(x, int) for x in lists])
    d = np.full(len(cand), np.inf)
    for k in range(0, len(pi), 2_000_000):
        sl = slice(k, k + 2_000_000)
        dd = point_triangle_distance(P[cand][pi[sl]], T[ti[sl], 0], T[ti[sl], 1], T[ti[sl], 2])
        np.minimum.at(d, pi[sl], dd)
    out[cand] = d
    return out


def enclosing_level(f, Mf, dist, lo, hi, eps=0., s=S_LEVEL):
    """Certified l with (real surface) subset {phi <= l}; dist(P) = exact distance to the real surface
    (eps = 0), or a value >= it minus eps."""
    delta = s * np.sqrt(3) / 2
    ax = [np.arange(lo[a] - s, hi[a] + s + 1e-9, s) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    G = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    band = G[dist(G) <= delta + eps]                       # contains the grid point nearest to any surface point
    v, g = f.eval(band, order=1)
    return float(np.max(v + np.linalg.norm(g, axis=1) * delta + .5 * Mf.eval(band) * delta ** 2))


def tight_level_mesh(f, Mf, V, Fc, eps=1e-4, max_edge=.004):
    """Tightest certified level for a triangle mesh: branch-and-bound over surface triangles.
    For a triangle with centroid c and circumradius-like rho, every point x of it satisfies
        phi(x) <= phi(c) + |grad phi(c)| rho + kappa(c) rho^2 / 2   (UB),
    and phi at vertices and centroids (points ON the surface) is attained (LB). Triangles with UB > LB + eps
    are split in four until none is left; then max_surface phi <= l = LB + eps <= max_surface phi + eps."""
    Vs, Fs = trimesh.remesh.subdivide_to_size(V, Fc, max_edge=max_edge)
    T = Vs[Fs]
    LB = float(f.eval(Vs).max())
    rounds = 0
    while len(T):
        c = T.mean(axis=1)
        rho = np.linalg.norm(T - c[:, None], axis=2).max(axis=1)
        v, g = f.eval(c, order=1)
        LB = max(LB, float(v.max()))
        UB = v + np.linalg.norm(g, axis=1) * rho + .5 * Mf.eval(c) * rho ** 2
        T = T[UB > LB + eps]
        if not len(T):
            break
        a, b, d = T[:, 0], T[:, 1], T[:, 2]
        ab, bd, da = (a + b) / 2, (b + d) / 2, (d + a) / 2
        LB = max(LB, float(f.eval(np.vstack([ab, bd, da])).max()))
        T = np.concatenate([np.stack(x, 1) for x in ((a, ab, da), (ab, b, bd), (da, bd, d), (ab, bd, da))])
        rounds += 1
    return LB + eps


def tight_level_implicit(f, Mf, sdf, lo, hi, eps=1e-4, h0=.008):
    """Tightest certified level for a body given by an exact SDF (1-Lipschitz): branch-and-bound over an
    octree of boxes that may contain surface points (|sdf(center)| <= half-diagonal). UB as for meshes with the
    half-diagonal; LB from exact surface points obtained by projecting outside centers along the gradient."""
    def project(P):
        d = sdf(P)
        P = P[d >= 0]
        if not len(P):
            return P
        e = 1e-6
        g = np.stack([(sdf(P + e * u) - sdf(P - e * u)) / (2 * e) for u in np.eye(3)], 1)
        Q = P - sdf(P)[:, None] * g / np.linalg.norm(g, axis=1, keepdims=True)
        return Q[np.abs(sdf(Q)) < 1e-9]

    ax = [np.arange(lo[a] - h0 / 2, hi[a] + h0, h0) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    C = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    h = h0
    C = C[np.abs(sdf(C)) <= h * np.sqrt(3) / 2]
    S = project(C)
    LB = float(f.eval(S).max())
    while len(C):
        r = h * np.sqrt(3) / 2
        v, g = f.eval(C, order=1)
        UB = v + np.linalg.norm(g, axis=1) * r + .5 * Mf.eval(C) * r ** 2
        C = C[UB > LB + eps]
        if not len(C):
            break
        h /= 2
        C = (C[:, None, :] + (CORNERS * h)[None]).reshape(-1, 3)
        C = C[np.abs(sdf(C)) <= h * np.sqrt(3) / 2]
        S = project(C)
        if len(S):
            LB = max(LB, float(f.eval(S).max()))
    return LB + eps


def spline_state(f):
    return dict(lo=f.lo, h=f.h, K=f.K, W=f.W)


def spline_from(d):
    f = Spline3.__new__(Spline3)
    f.lo, f.h, f.K, f.W = np.asarray(d['lo']), float(d['h']), np.asarray(d['K']), np.asarray(d['W'])
    f.n = f.K + 3
    return f


# ---------------------------------------------------------------------------------------------
# offline setup
# ---------------------------------------------------------------------------------------------
def build_link(i):
    vis = trimesh.load(os.path.join(MESH, 'visual', f'link{i}_vis.stl'))
    wt = trimesh.load(os.path.join(MESH, 'voxel_128', f'link{i}.stl'))          # watertight: sign only
    dense = trimesh.sample.sample_surface_even(vis, int(vis.area * 2e6), seed=0)[0]
    tree = cKDTree(np.vstack([dense, vis.vertices]))         # approximate distance: fitting only
    vox = wt.voxelized(.002).fill()
    lo, hi = vis.bounds

    def sdf(P):
        d = tree.query(P)[0]
        return np.where(vox.is_filled(P), -d, d)

    f = Spline3(lo - .06, hi + .06, .012).fit(sdf, .006)
    Mf = hessian_majorant(f)
    V, Fc = np.asarray(vis.vertices), np.asarray(vis.faces)
    Vs, Fs = trimesh.remesh.subdivide_to_size(V, Fc, max_edge=.004)          # same surface, small triangles
    l = tight_level_mesh(f, Mf, V, Fc)                  # tightest certified level (branch-and-bound on the mesh)
    blo, side = surface_boxes(f, BOX_SIDE, level=l)
    T = Vs[Fs]
    cen = T.mean(axis=1)
    return dict(f=spline_state(f), l=l, PC=blo + side / 2, side=side, gt=cen,
                gt_r=np.linalg.norm(T - cen[:, None], axis=2).max(axis=1), V=V, F=Fc)


def build():
    if os.path.exists(CACHE):
        d = np.load(CACHE, allow_pickle=True)['d'].item()
    else:
        t0 = time.perf_counter()
        d = dict(links=[build_link(i) for i in LINK_IDS], obst={})
        for name, B in OBSTACLES.items():
            lo = np.min([np.array(c) - b for c, b in B], axis=0)
            hi = np.max([np.array(c) + b for c, b in B], axis=0)
            f = Spline3(lo - .15, hi + .15, .02).fit(lambda P, b=B: sdf_boxes(P, b), .01)
            l = tight_level_implicit(f, hessian_majorant(f), lambda P, b=B: sdf_boxes(P, b), lo, hi)
            d['obst'][name] = dict(f=spline_state(f), l=l, blo=lo, bhi=hi)
        d['t_setup'] = time.perf_counter() - t0
        np.savez_compressed(CACHE, d=np.array(d, dtype=object))
    links = []
    for i, L in zip(LINK_IDS, d['links']):
        r = L['side'] * np.sqrt(3) / 2
        groups, cen, crad = clusters(L['PC'], r)
        links.append(dict(L, frame=i, n_joints=min(i, DOF), r=r, groups=groups, cc=cen, crad=crad,
                          rad=float(np.linalg.norm(L['gt'], axis=1).max()) + .02))
    obst = {}
    for name, o in d['obst'].items():
        f = spline_from(o['f'])
        Gn, Mn = neighborhood_bounds(f, max(L['crad'].max() for L in links))
        obst[name] = dict(o, f=f, M=hessian_majorant(f), Gn=Gn, Mn=Mn, boxes=OBSTACLES[name],
                          dlo=f.lo, dhi=f.lo + f.K * f.h)
    info = dict(t_setup=d['t_setup'],
                links=[dict(link=L['frame'], l_mm=1e3 * L['l'], n_boxes=len(L['PC']),
                            side_mm=1e3 * L['side']) for L in links],
                obstacles={k: dict(l_mm=1e3 * o['l']) for k, o in obst.items()})
    return links, obst, info


# ---------------------------------------------------------------------------------------------
# barrier rows and the real distance
# ---------------------------------------------------------------------------------------------
def barrier_rows(q, links, obst):
    T, Z, O = fk(q)
    ZxO = np.cross(Z, O)
    rows, hs = [], []
    for L in links:
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        for o in obst.values():
            if np.any(p + L['rad'] < o['dlo']) or np.any(p - L['rad'] > o['dhi']):
                continue                                         # link sphere outside the obstacle domain
            # coarse pruning: a cluster ball (center cc, radius R) inside the domain has, for every box corner,
            # h_v >= phi(cc) - G R - M r^2 / 2 - l_O  (certified neighborhood bounds G, M)
            ccw = L['cc'] @ R.T + p
            keep = ~np.all((ccw - L['crad'][:, None] > o['dlo']) & (ccw + L['crad'][:, None] < o['dhi']), axis=1)
            ins = np.flatnonzero(~keep)
            if len(ins):
                ci = o['f'].cell_of(ccw[ins])
                low = (o['f'].eval(ccw[ins]) - o['Gn'][ci[:, 0], ci[:, 1], ci[:, 2]] * L['crad'][ins]
                       - .5 * o['Mn'][ci[:, 0], ci[:, 1], ci[:, 2]] * L['r'] ** 2 - o['l'])
                keep[ins[low < ACT]] = True
            if not keep.any():
                continue
            sub = np.concatenate([L['groups'][j] for j in np.flatnonzero(keep)])
            Cs = L['PC'][sub] @ R.T + p
            idx = np.flatnonzero(np.all((Cs - L['r'] > o['dlo']) & (Cs + L['r'] < o['dhi']), axis=1))
            Cw = Cs
            if not len(idx):                                     # boxes outside the domain are >= 14 cm away
                continue
            val, g = o['f'].eval(Cw[idx], order=1)
            marg = .5 * o['M'].eval(Cw[idx]) * L['r'] ** 2
            near = idx[val - np.linalg.norm(g, axis=1) * L['r'] - marg - o['l'] < ACT]
            if not len(near):
                continue
            c = Cw[near]
            v2, g2, H2 = o['f'].eval(c, order=2)
            Mv, gM = o['M'].eval(c, order=1)
            d = (CORNERS * L['side']) @ R.T                           # corner offsets in the world (8,3)
            hh = v2[:, None] + g2 @ d.T - .5 * Mv[:, None] * L['r'] ** 2 - o['l']
            bi, vi = np.nonzero(hh < ACT)
            if not len(bi):
                continue
            dd = d[vi]
            a = g2[bi] + np.einsum('nij,nj->ni', H2[bi], dd) - .5 * L['r'] ** 2 * gM[bi]
            row = (np.cross(c[bi], a) + np.cross(dd, g2[bi])) @ Z.T - a @ ZxO.T
            row[:, L['n_joints']:] = 0.                                 # later joints do not move the link
            rows.append(row); hs.append(hh[bi, vi])
    if not rows:
        return np.zeros((0, DOF)), np.zeros(0)
    return np.vstack(rows), np.concatenate(hs)


def real_gap(q, links, obst, per_link=False):
    """Certified lower bound on the distance between the real link meshes and the real obstacles: the
    obstacle SDF is exact and 1-Lipschitz, so over a triangle with centroid m and radius rho it is at
    least sdf(m) - rho (triangles subdivided to edges <= 4 mm; no sampling)."""
    T, _, _ = fk(q)
    gaps = []
    for L in links:
        p = T[L['frame']][:3, 3]
        near = [o for o in obst.values()
                if np.all(p + L['rad'] > o['blo'] - .05) and np.all(p - L['rad'] < o['bhi'] + .05)]
        g = np.inf
        if near:
            W = L['gt'] @ T[L['frame']][:3, :3].T + p
        for o in near:
            m = np.all((W > o['blo'] - .05) & (W < o['bhi'] + .05), axis=1)
            if m.any():
                g = min(g, float((sdf_boxes(W[m], o['boxes']) - L['gt_r'][m]).min()))
        gaps.append(g)
    return gaps if per_link else min(gaps)


# ---------------------------------------------------------------------------------------------
# closed loop
# ---------------------------------------------------------------------------------------------
BOX_G = np.vstack([np.eye(DOF), -np.eye(DOF)])
BOX_H = -QD_MAX * np.ones(2 * DOF)


def simulate(q0, qg, links, obst, method, steps=1000, record=False):
    q = np.array(q0, float)
    log = dict(t=[], rows=[], h=[], gap=[], slack=0, reached=None, q=[], boxes=[])
    u = z = None
    for k in range(steps):
        u_nom = np.clip(KQ * (qg - q), -QD_MAX, QD_MAX)
        t0 = time.perf_counter()
        if method == 'ours':                  # summed-field barriers with Bernstein coefficient rows (summed3d.py)
            import summed as SM
            import summed3d as S3
            Ar, Tr, Cr, nb, hl, Er = S3.franka_rows(q, links, obst)
            u, s, z = SM.solve(u_nom, Ar, Tr, Cr, Er, BOX_G, BOX_H, z)
            log['slack'] += s > 0
            log['rows'].append(len(Cr)); log['h'].append(hl); log['boxes'].append(nb)
        elif method == 'corner':              # previous corner barriers
            G, h = barrier_rows(q, links, obst)
            u, s = solve_ldp_qp(u_nom, np.vstack([G, BOX_G]), np.r_[-GAMMA * h, BOX_H], len(h), u)
            log['slack'] += s > 0
            log['rows'].append(len(h)); log['h'].append(float(h.min()) if len(h) else np.inf)
        else:
            u = u_nom
        log['t'].append(time.perf_counter() - t0)
        log['gap'].append(real_gap(q, links, obst))
        if record:
            log['q'].append(q.copy())
        q = q + DT * u
        if np.linalg.norm(qg - q) < .05:
            log['reached'] = (k + 1) * DT
            break
    if record:
        log['q'].append(q.copy())
    return log


def sample_config(rng, links, obst, clear=.03):
    mid, half = (Q_MIN + Q_MAX) / 2, (Q_MAX - Q_MIN) / 2 * .8
    while True:
        q = rng.uniform(mid - half, mid + half)
        T, _, _ = fk(q)
        if T[8][2, 3] < .1:                                        # keep the hand above the floor level
            continue
        if real_gap(q, links, obst) > clear and barrier_rows(q, links, obst)[1].min(initial=1.) > 0:
            return q


def main():
    n_trials = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    links, obst, info = build()
    print('setup', json.dumps(info), flush=True)
    rng = np.random.default_rng(0)
    trials = []
    reuse = os.path.join(HERE, 'franka_trials.json')    # fixed trial set (start, goal, unfiltered min gap)
    if os.path.exists(reuse):
        trials = [(np.array(a), np.array(b), g) for a, b, g in json.load(open(reuse))][:n_trials]
    while len(trials) < n_trials:                       # keep pairs whose unfiltered motion collides
        q0, qg = sample_config(rng, links, obst), sample_config(rng, links, obst)
        nom = simulate(q0, qg, links, obst, 'nominal')
        if min(nom['gap']) < 0:
            trials.append((q0, qg, min(nom['gap'])))
    res = []
    import summed3d as S3
    S3.prep_all(links, obst.values())
    for i, (q0, qg, gnom) in enumerate(trials):
        L = simulate(q0, qg, links, obst, 'ours')
        r = dict(trial=i, q0=q0.tolist(), qg=qg.tolist(), nominal_min_gap_mm=1e3 * gnom,
                 min_gap_mm=1e3 * min(L['gap']), reached=L['reached'], slack_steps=int(L['slack']),
                 min_h_mm=1e3 * min(L['h']), t_med_ms=1e3 * float(np.median(L['t'])),
                 t_max_ms=1e3 * float(np.max(L['t'])), rows_max=int(max(L['rows'])),
                 boxes_max=int(max(L['boxes'])))
        res.append(r)
        print(i, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k not in ('q0', 'qg')},
              flush=True)
    gaps = np.array([r['min_gap_mm'] for r in res])
    summ = dict(n=len(res), collisions=int((gaps < 0).sum()), min_gap_mm=float(gaps.min()),
                median_min_gap_mm=float(np.median(gaps)), nominal_collisions=len(res),
                reached=int(sum(r['reached'] is not None for r in res)),
                slack_trials=int(sum(r['slack_steps'] > 0 for r in res)),
                t_med_ms=float(np.median([r['t_med_ms'] for r in res])),
                t_max_ms=float(max(r['t_max_ms'] for r in res)),
                rows_max=int(max(r['rows_max'] for r in res)), boxes_max=int(max(r['boxes_max'] for r in res)))
    print('summary', json.dumps(summ), flush=True)
    json.dump(dict(setup=info, summary=summ, trials=res), open(os.path.join(HERE, 'franka_results.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
