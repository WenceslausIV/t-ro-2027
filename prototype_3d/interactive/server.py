"""
Interactive safety-filter demo: Franka Panda with surface-cover barriers (paper, Sec. VI).

    python prototype_3d/interactive/server.py        then open  http://localhost:8765

Browser (index.html, three.js): drag the end-effector goal and the objects with a 3D gizmo.
  Move          arm 1 follows the goal (damped-least-squares IK) through the safety filter
  Add object    a non-convex object (ground truth + translucent certified SDF surface) in front of the arm
  Dual arm      second Panda facing arm 1 (its links are SDF bodies for arm 1)
  Play arm 2    arm 2 moves between random configurations, unfiltered (the "adversary")
  Hide GT/SDF   show only the certified SDF surfaces / only the ground-truth shapes
Arm 1 is filtered against every object and every link of arm 2; the known motion of arm 2 and of the
dragged objects enters the CBF condition as a drift term.
"""
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import dual3d as D                                    # noqa: E402
import franka3d as F                                  # noqa: E402
from franka_gif import decimate, level_mesh           # noqa: E402
from proto3d import Spline3, sdf_boxes, rotm          # noqa: E402
from sdf_cbf_utils import solve_ldp_qp                # noqa: E402
from skimage import measure                           # noqa: E402

PORT, DT, GAMMA, ACT, QD_MAX = 8765, .01, 5., .05, 1.
Q_HOME = np.array([0., -.5, 0., -2.2, 0., 1.8, .785])
Q2_HOME = np.array([0., -1.2, 0., -2.8, 0., 1.7, .785])     # arm 2 starts folded, clear of arm 1
TCP = np.array([0., 0., .1034])
OBJECT_LIBRARY = {                                    # local frames, rounded-box unions (center, half extents)
    'hook': [((-.2, 0., -.13), (.04, .04, .30)), ((0., 0., .13), (.24, .04, .04)), ((.2, 0., -.01), (.035, .035, .16))],
    'U-bracket': [((0., 0., -.12), (.16, .14, .02)), ((-.14, 0., .03), (.02, .14, .15)), ((.14, 0., .03), (.02, .14, .15))],
    'rack': [((0., 0., -.14), (.25, .12, .02)), ((0., 0., .14), (.25, .012, .012)), ((.235, 0., 0.), (.015, .015, .14))],
    'L-beam': [((0., 0., 0.), (.03, .2, .03)), ((0., .17, .14), (.03, .03, .17))],
    'arch': [((0., -.2, -.05), (.035, .035, .25)), ((0., .2, -.05), (.035, .035, .25)), ((0., 0., .23), (.035, .235, .035))],
}
CACHE = os.path.join(HERE, 'cache_objects.npz')
SETUPS = os.path.join(HERE, 'setups')


def shape_spec(name):
    """Exact SDF (local frame), bounding box, and spline cell size of a library object."""
    if name in ('star tube', 'star tube (rounded)'):  # star prism with a star hole, axis along local x
        import dolphin3d as Dp
        ext = np.array([Dp.HALF_T, Dp.R_OUT, Dp.R_OUT])
        fun = Dp.hoop_sdf if name == 'star tube' else Dp.hoop_sdf_rounded
        return dict(sdf=lambda P: fun(np.atleast_2d(P) + Dp.HOOP_C), lo=-ext, hi=ext, h=.015)
    B = OBJECT_LIBRARY[name]
    return dict(sdf=lambda P, b=B: sdf_boxes(P, b), lo=np.min([np.array(c) - b for c, b in B], axis=0),
                hi=np.max([np.array(c) + b for c, b in B], axis=0), h=.02)


SHAPES = list(OBJECT_LIBRARY) + ['star tube', 'star tube (rounded)']


# ---------------------------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------------------------
def build_objects():
    d = np.load(CACHE, allow_pickle=True)['d'].item() if os.path.exists(CACHE) else {}
    missing = [n for n in SHAPES if n not in d]
    for name in missing:
        sp = shape_spec(name)
        f = Spline3(sp['lo'] - .15, sp['hi'] + .15, sp['h']).fit(sp['sdf'], sp['h'] / 2)
        l = F.tight_level_implicit(f, F.hessian_majorant(f), sp['sdf'], sp['lo'], sp['hi'])
        d[name] = dict(f=F.spline_state(f), l=l, blo=sp['lo'], bhi=sp['hi'])
        print(f'object {name}: level {1e3 * l:.1f} mm', flush=True)
    if missing:
        np.savez_compressed(CACHE, d=np.array(d, dtype=object))
    objs = {}
    for name in SHAPES:
        o = d[name]
        f = F.spline_from(o['f'])
        Gn, Mn = F.neighborhood_bounds(f, .05)
        objs[name] = dict(o, fs=f, Ms=F.hessian_majorant(f), Gn=Gn, Mn=Mn, sdf=shape_spec(name)['sdf'],
                          dlo=f.lo, dhi=f.lo + f.K * f.h)
    return objs


def mesh_json(V, Fc):
    """Binary (base64 float32 / uint32) mesh for the browser."""
    import base64
    return dict(v64=base64.b64encode(np.ascontiguousarray(V, np.float32).tobytes()).decode(),
                f64=base64.b64encode(np.ascontiguousarray(Fc, np.uint32).tobytes()).decode())


def level_surface(fun, grad, level, lo, hi, res=.0025, iters=3):
    """Display mesh of {fun = level}: marching cubes on a res grid, then Newton projection of every vertex onto
    the level set, so the drawn surface is the certified one up to the chord error of res-sized triangles."""
    lo, hi = np.asarray(lo) - 3 * res, np.asarray(hi) + 3 * res
    ax = [np.arange(lo[a], hi[a] + 1e-9, res) for a in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing='ij')
    V = fun(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)
    v, fc, _, _ = measure.marching_cubes(V, level, spacing=(res,) * 3)
    v = v + lo
    for _ in range(iters):
        r, g = fun(v) - level, grad(v)
        v = v - (r / np.maximum((g * g).sum(1), 1e-12))[:, None] * g
    return v, fc


def display_res(lo, hi):
    """Grid of the display surfaces: vertices are projected exactly onto the surface, so 5-mm triangles only
    add a chord error below 0.1 mm while keeping the browser's triangle count low."""
    return max(.005, float(np.max(np.asarray(hi) - np.asarray(lo))) / 150)


def spline_surface(f, level, lo, hi):
    return level_surface(f.eval, lambda P: f.eval(P, order=1)[1], level, lo, hi, res=display_res(lo, hi))


def exact_surface(o):
    """Ground truth of an object: its exact zero level set (projected with the exact SDF)."""
    fun = o['sdf']
    e = 1e-6
    grad = lambda P: np.stack([(fun(P + e * u) - fun(P - e * u)) / (2 * e) for u in np.eye(3)], 1)
    return level_surface(fun, grad, 0., o['blo'], o['bhi'], res=display_res(o['blo'], o['bhi']))


def geometry_payload(links, objs):
    import trimesh
    base = trimesh.load(os.path.join(F.MESH, 'visual', 'link0_vis.stl'))
    out = dict(base=mesh_json(np.asarray(base.vertices), np.asarray(base.faces)),
               links=[dict(frame=L['frame'], gt=mesh_json(L['V'], L['F']),       # the real (certified) mesh
                           sdf=mesh_json(*spline_surface(L['fs'], L['l'], L['V'].min(0) - L['l'],
                                                         L['V'].max(0) + L['l'])),
                           level_mm=1e3 * L['l']) for L in links],
               objects={name: dict(gt=mesh_json(*exact_surface(o)),
                                   sdf=mesh_json(*spline_surface(o['fs'], o['l'], o['blo'] - .01, o['bhi'] + .01)),
                                   level_mm=1e3 * o['l'])
                        for name, o in objs.items()},
               base2=D.BASE2.tolist())
    return out


# ---------------------------------------------------------------------------------------------
# barriers of arm-1 links against a body B (object or arm-2 link) with pose (RB, pB)
# ---------------------------------------------------------------------------------------------
def body_rows(A, TA, Z1, Z1xO1, B, RB, pB):
    """Corner barriers of link A (frame TA) against B. Returns rows for qd1 (n,7), h (n,), and the vectors
    (aw, w) with which B's own motion enters: dh/dt += -(aw.v_B + w_B.(w - p_B x aw)) for a twist (v_B, w_B)."""
    RA, pA = TA[:3, :3], TA[:3, 3]
    if np.linalg.norm(pA - pB) > A['rad'] + B['reach'] + ACT:
        return None
    RbA, tb = RB.T @ RA, RB.T @ (pA - pB)
    cc = A['cc'] @ RbA.T + tb
    keep = ~np.all((cc - A['crad'][:, None] > B['dlo']) & (cc + A['crad'][:, None] < B['dhi']), axis=1)
    ins = np.flatnonzero(~keep)
    if len(ins):
        ci = B['fs'].cell_of(cc[ins])
        low = (B['fs'].eval(cc[ins]) - B['Gn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['crad'][ins]
               - .5 * B['Mn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['r'] ** 2 - B['l'])
        keep[ins[low < ACT]] = True
    if not keep.any():
        return None
    sub = np.concatenate([A['groups'][j] for j in np.flatnonzero(keep)])
    C = A['PC'][sub] @ RbA.T + tb
    idx = np.flatnonzero(np.all((C - A['r'] > B['dlo']) & (C + A['r'] < B['dhi']), axis=1))
    if not len(idx):
        return None
    val, g = B['fs'].eval(C[idx], order=1)
    near = idx[val - np.linalg.norm(g, axis=1) * A['r'] - .5 * B['Ms'].eval(C[idx]) * A['r'] ** 2 - B['l'] < ACT]
    if not len(near):
        return None
    c = C[near]
    v2, g2, H2 = B['fs'].eval(c, order=2)
    Mv, gM = B['Ms'].eval(c, order=1)
    d = (F.CORNERS * A['side']) @ RbA.T
    hh = v2[:, None] + g2 @ d.T - .5 * Mv[:, None] * A['r'] ** 2 - B['l']
    bi, vi = np.nonzero(hh < ACT)
    if not len(bi):
        return None
    dd = d[vi]
    a = g2[bi] + np.einsum('nij,nj->ni', H2[bi], dd) - .5 * A['r'] ** 2 * gM[bi]
    aw, gw = a @ RB.T, g2[bi] @ RB.T
    cw, dw = c[bi] @ RB.T + pB, dd @ RB.T
    w = np.cross(cw, aw) + np.cross(dw, gw)
    r1 = w @ Z1.T - aw @ Z1xO1.T
    r1[:, A['n_joints']:] = 0.
    return r1, hh[bi, vi], aw, w


# ---------------------------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------------------------
def quat_to_R(q):
    x, y, z, w = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def _log(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    th = np.arccos(c)
    if th < 1e-9:
        return np.zeros(3)
    return th / (2 * np.sin(th)) * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])


class Sim:
    def __init__(self):
        t0 = time.perf_counter()
        self.links, _ = D.build()
        for L in self.links:
            L['reach'] = L['rad']
            import trimesh
            Vs, Fs = trimesh.remesh.subdivide_to_size(L['V'], L['F'], max_edge=.008)   # display distances
            Tr = Vs[Fs]
            L['gt8'] = Tr.mean(axis=1)
            L['gt8_r'] = np.linalg.norm(Tr - L['gt8'][:, None], axis=2).max(axis=1)
        self.objlib = build_objects()
        for o in self.objlib.values():
            o['reach'] = float(np.linalg.norm(np.r_[o['blo'], o['bhi']].reshape(2, 3), axis=1).max()) + .05
        print(f'setup {time.perf_counter() - t0:.1f} s', flush=True)
        self.geometry = json.dumps(geometry_payload(self.links, self.objlib))
        self.lock = threading.Lock()
        self.reset()
        self.snapshot = json.dumps(self.state())

    def reset(self):
        with getattr(self, 'lock', threading.Lock()):
            self.q1 = Q_HOME.copy()
            self.q2 = Q2_HOME.copy()
            self.qd2 = np.zeros(7)
            self.goal_q2 = Q2_HOME.copy()
            self.moving = False
            self.dual = False
            self.play2 = False
            self.objects = []                            # dict(id, type, p, R, p_prev, R_prev)
            T, _, _ = F.fk(self.q1)
            self.goal_p, self.goal_R = T[8][:3, :3] @ TCP + T[8][:3, 3], T[8][:3, :3].copy()
            self.metrics = {}
            self.u_prev = None
            self.t_last_goal2 = 0.
            self.n_obj = 0
            self.n_box = 0
            self.last_saved = None

    # --- commands from the browser
    def command(self, cmd):
        with self.lock:
            if cmd == 'move':
                self.moving = not self.moving
            elif cmd in ('add_object', 'add_star'):
                names = list(OBJECT_LIBRARY)
                name = 'star tube' if cmd == 'add_star' else names[self.n_box % len(names)]
                self.n_box += cmd == 'add_object'
                self.n_obj += 1
                p = self.free_spot(name)
                self.objects.append(dict(id=self.n_obj, type=name, p=p, R=np.eye(3), p_prev=p.copy(), R_prev=np.eye(3)))
            elif cmd == 'save':
                self.last_saved = self.save_setup()
            elif cmd == 'remove_objects':
                self.objects = []
            elif cmd == 'dual':
                self.dual = not self.dual
                self.play2 = self.play2 and self.dual
                self.q2 = Q2_HOME.copy(); self.goal_q2 = Q2_HOME.copy(); self.qd2 = np.zeros(7)
            elif cmd == 'play2':
                self.play2 = self.dual and not self.play2
                self.qd2 = np.zeros(7)
            elif cmd == 'goal_to_ee':
                T, _, _ = F.fk(self.q1)
                self.goal_p, self.goal_R = T[8][:3, :3] @ TCP + T[8][:3, 3], T[8][:3, :3].copy()
            elif cmd == 'reset':
                pass
        if cmd == 'reset':
            self.reset()

    def save_setup(self):
        """Write the current scene for replication (called with the lock held): arm-1 configuration (the start of
        an experiment), end-effector goal pose, objects (type, position, rotation), and arm 2 if present."""
        os.makedirs(SETUPS, exist_ok=True)
        stamp = time.strftime('%Y%m%d_%H%M%S')
        T, _, _ = F.fk(self.q1)
        d = dict(saved=stamp, q1=self.q1.tolist(),
                 ee=dict(p=(T[8][:3, :3] @ TCP + T[8][:3, 3]).tolist(), R=T[8][:3, :3].tolist()),
                 goal=dict(p=self.goal_p.tolist(), R=self.goal_R.tolist()),
                 objects=[dict(type=o['type'], p=o['p'].tolist(), R=o['R'].tolist()) for o in self.objects],
                 dual=self.dual, q2=self.q2.tolist() if self.dual else None, base2=D.BASE2.tolist(),
                 tcp=TCP.tolist())
        path = os.path.join(SETUPS, f'setup_{stamp}.json')
        json.dump(d, open(path, 'w'), indent=1)
        print('saved', path, flush=True)
        return os.path.relpath(path, os.path.dirname(os.path.dirname(HERE)))

    def free_spot(self, name):
        """Candidate position around the arm where the new object activates no barrier of arm 1 and is farthest
        (bounding spheres) from the other objects and from arm 2 (called with the lock held)."""
        B = self.objlib[name]
        T1, Z1, O1 = F.fk(self.q1)
        Z1xO1 = np.cross(Z1, O1)
        cands = [np.array([x, y, z]) for x in (.55, .75, .35, .15) for y in (-.5, .5, 0., -.3, .3, -.7, .7)
                 for z in (.35, .6, .2)]
        spheres = [(o['p'], self.objlib[o['type']]['reach'] - .05) for o in self.objects]
        if self.dual:
            (_, _, _), (T2, _, _) = D.fk2(np.r_[self.q1, self.q2])
            spheres += [(T2[L['frame']][:3, 3], L['rad']) for L in self.links]
        best, best_score = None, -np.inf
        for p in cands:
            score = min((np.linalg.norm(p - c) - r - (B['reach'] - .05) for c, r in spheres), default=1.)
            if score <= best_score:
                continue
            if all(body_rows(A, T1[A['frame']], Z1, Z1xO1, B, np.eye(3), p) is None for A in self.links):
                best, best_score = p, score
            if best_score > .05:
                break
        return best if best is not None else cands[0]

    def set_input(self, inp):
        with self.lock:
            if 'goal' in inp:
                self.goal_p = np.array(inp['goal']['p'], float)
                self.goal_R = quat_to_R(inp['goal']['q'])
            for o in inp.get('objects', []):
                for ob in self.objects:
                    if ob['id'] == o['id']:
                        ob['p'] = np.array(o['p'], float)
                        ob['R'] = quat_to_R(o['q'])

    # --- one control step (dt: the wall-clock time of the previous step, 10-30 ms, so motion looks real-time)
    def step(self, dt=DT):
        with self.lock:
            q1, q2 = self.q1.copy(), self.q2.copy()
            goal_p, goal_R, moving = self.goal_p.copy(), self.goal_R.copy(), self.moving
            objects = [dict(o) for o in self.objects]
            for o in self.objects:
                o['p_prev'], o['R_prev'] = o['p'].copy(), o['R'].copy()
            dual, play2 = self.dual, self.play2
        t0 = time.perf_counter()
        # arm 2: random joint targets, unfiltered
        if dual and play2:
            now = time.time()
            if np.linalg.norm(self.goal_q2 - q2) < .1 or now - self.t_last_goal2 > 4.:
                mid, half = (F.Q_MIN + F.Q_MAX) / 2, (F.Q_MAX - F.Q_MIN) / 2 * .7
                self.goal_q2 = np.random.uniform(mid - half, mid + half)
                self.t_last_goal2 = now
            qd2 = np.clip(1.2 * (self.goal_q2 - q2), -.8, .8)
        else:
            qd2 = np.zeros(7)
        # arm 1 nominal: damped-least-squares IK toward the goal pose, null space toward home
        T1, Z1, O1 = F.fk(q1)
        p_ee = T1[8][:3, :3] @ TCP + T1[8][:3, 3]
        Jv = np.cross(Z1, p_ee - O1).T
        Jw = Z1.T
        J = np.vstack([Jv, Jw])
        if moving:
            e = np.r_[np.clip(2. * (goal_p - p_ee), -.5, .5), np.clip(1.5 * _log(goal_R @ T1[8][:3, :3].T), -1., 1.)]
            Jp = J.T @ np.linalg.inv(J @ J.T + 1e-3 * np.eye(6))
            u_nom = Jp @ e + (np.eye(7) - Jp @ J) @ (.5 * (Q_HOME - q1))
            u_nom = np.clip(u_nom, -QD_MAX, QD_MAX)
        else:
            u_nom = np.zeros(7)
        # barriers
        Z1xO1 = np.cross(Z1, O1)
        rows, rhs = [], []
        min_h = np.inf
        bodies = []
        for o in objects:
            B = self.objlib[o['type']]
            v_B = (o['p'] - o['p_prev']) / dt
            w_B = _log(o['R'] @ o['R_prev'].T) / dt
            v_B, w_B = np.clip(v_B, -2, 2), np.clip(w_B, -3, 3)
            bodies.append((B, o['R'], o['p'], ('twist', v_B, w_B)))
        if dual:
            (_, _, _), (T2, Z2, O2) = D.fk2(np.r_[q1, q2])
            for L2 in self.links:
                bodies.append((L2, T2[L2['frame']][:3, :3], T2[L2['frame']][:3, 3], ('arm2', L2, Z2, O2)))
        for A in self.links:
            for B, RB, pB, motion in bodies:
                res = body_rows(A, T1[A['frame']], Z1, Z1xO1, B, RB, pB)
                if res is None:
                    continue
                r1, h, aw, w = res
                if motion[0] == 'twist':
                    _, v_B, w_B = motion
                    e = -(aw @ v_B + (w - np.cross(pB, aw)) @ w_B)
                else:
                    _, L2, Z2, O2 = motion
                    r2 = -(w @ Z2.T - aw @ np.cross(Z2, O2).T)
                    r2[:, L2['n_joints']:] = 0.
                    e = r2 @ qd2
                rows.append(r1); rhs.append(-GAMMA * h - e)
                min_h = min(min_h, float(h.min()))
        lo = np.maximum(-QD_MAX, -2. * (q1 - F.Q_MIN))
        hi = np.minimum(QD_MAX, 2. * (F.Q_MAX - q1))
        G = np.vstack(rows + [np.eye(7), -np.eye(7)])
        hvec = np.r_[np.concatenate(rhs) if rhs else np.zeros(0), lo, -hi]
        n_bar = int(sum(len(r) for r in rows))
        u, slack = solve_ldp_qp(u_nom, G, hvec, n_bar, self.u_prev)
        self.u_prev = u
        t_ctrl = time.perf_counter() - t0
        with self.lock:
            self.q1 = q1 + dt * u
            self.q2 = q2 + dt * qd2
            self.metrics = dict(t_ctrl_ms=1e3 * t_ctrl, barriers=n_bar, min_h_mm=1e3 * min_h if n_bar else None,
                                slack=bool(slack > 0), ee_err_mm=1e3 * float(np.linalg.norm(goal_p - p_ee)),
                                filter_active=bool(n_bar and np.linalg.norm(u - u_nom) > 1e-3))

    def display_gaps(self):
        """Distances for display: certified lower bound to objects (exact object SDF at the centroids of 8-mm
        mesh triangles minus their radii; the spline field preselects the links and triangles that can be within
        5 cm), sampled (approximate) mesh distance to arm 2, and the arm-1 links that touch an object."""
        with self.lock:
            q1, q2, objects, dual = self.q1.copy(), self.q2.copy(), [dict(o) for o in self.objects], self.dual
        T1, _, _ = F.fk(q1)
        g_obj, hit1 = np.inf, set()
        for i, L in enumerate(self.links):
            TA = T1[L['frame']]
            for o in objects:
                B = self.objlib[o['type']]
                if np.linalg.norm(TA[:3, 3] - o['p']) > L['rad'] + B['reach']:
                    continue
                Rl, tl = o['R'].T @ TA[:3, :3], o['R'].T @ (TA[:3, 3] - o['p'])   # link -> object frame
                C = L['PC'] @ Rl.T + tl
                inb = np.all((C > B['dlo']) & (C < B['dhi']), axis=1)
                if not inb.any() or B['fs'].eval(C[inb]).min() - L['r'] > .06:
                    continue
                P = L['gt8'] @ Rl.T + tl
                m = np.flatnonzero(np.all((P > B['dlo']) & (P < B['dhi']), axis=1))
                if len(m):
                    m = m[B['fs'].eval(P[m]) < .06]
                if not len(m):
                    continue
                g = float((B['sdf'](P[m]) - L['gt8_r'][m]).min())
                g_obj = min(g_obj, g)
                if g < 0:
                    hit1.add(i)
        g_arm = np.inf
        if dual:
            g_arm = D.approx_gap(np.r_[q1, q2], self.links)
        return (g_obj if g_obj < .05 else np.inf), g_arm, sorted(hit1)

    def gap_loop(self):
        """Display distances in their own thread, so the control loop never waits for them."""
        while True:
            t0 = time.perf_counter()
            try:
                g_obj, g_arm, hit = self.display_gaps()
                self.gaps = dict(gap_obj_mm=None if not np.isfinite(g_obj) else 1e3 * g_obj,
                                 gap_arm2_mm=None if not np.isfinite(g_arm) else 1e3 * g_arm, hit_links=hit,
                                 sim_hz=getattr(self, '_hz', None))
            except Exception as ex:
                print('gap error:', repr(ex), flush=True)
            time.sleep(max(.05, .3 - (time.perf_counter() - t0)))

    def state(self):
        with self.lock:
            T1, _, _ = F.fk(self.q1)
            (_, _, _), (T2, _, _) = D.fk2(np.r_[self.q1, self.q2])
            return dict(arm1=[T1[L['frame']].tolist() for L in self.links],
                        arm2=[T2[L['frame']].tolist() for L in self.links] if self.dual else None,
                        objects=[dict(id=o['id'], type=o['type'], p=o['p'].tolist(), R=o['R'].tolist())
                                 for o in self.objects],
                        moving=self.moving, dual=self.dual, play2=self.play2, saved=self.last_saved, t=getattr(self, 't_sim', 0.),
                        goal=dict(p=self.goal_p.tolist(), R=self.goal_R.tolist()),
                        metrics=dict(self.metrics, **getattr(self, 'gaps', {})))

    def run(self):
        threading.Thread(target=self.gap_loop, daemon=True).start()
        dt, t_last = DT, time.perf_counter()
        while True:
            t0 = time.perf_counter()
            try:
                self.step(dt)
                self.t_sim = getattr(self, 't_sim', 0.) + dt  # simulation clock (browser interpolates on it)
                self.snapshot = json.dumps(self.state())  # served as is: requests do no work
            except Exception as ex:                     # keep the loop alive
                print('step error:', repr(ex), flush=True)
                time.sleep(.2)
            time.sleep(max(0., DT - (time.perf_counter() - t0)))
            now = time.perf_counter()
            dt = min(max(DT, now - t_last), .05)        # integrate the real elapsed time (motion stays real-time)
            t_last = now
            self._hz = .95 * getattr(self, '_hz', 1. / DT) + .05 / max(now - t0, 1e-4)


SIM = None


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, body, ctype='application/json'):
        b = body.encode() if isinstance(body, str) else body
        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(b)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(b)

    def do_GET(self):
        if self.path in ('/', '/index.html'):
            self._send(open(os.path.join(HERE, 'index.html'), 'rb').read(), 'text/html; charset=utf-8')
        elif self.path == '/geometry':
            self._send(SIM.geometry)
        elif self.path == '/state':
            self._send(SIM.snapshot)
        else:
            self.send_error(404)

    def do_POST(self):
        n = int(self.headers.get('Content-Length', 0))
        data = json.loads(self.rfile.read(n) or b'{}')
        if self.path == '/input':
            SIM.set_input(data)
            self._send(SIM.snapshot)
        elif self.path == '/cmd':
            SIM.command(data.get('cmd'))
            self._send(json.dumps(SIM.state()))
        else:
            self.send_error(404)


def main():
    global SIM
    SIM = Sim()
    threading.Thread(target=SIM.run, daemon=True).start()
    print(f'open http://localhost:{PORT}', flush=True)
    ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()


if __name__ == '__main__':
    main()
