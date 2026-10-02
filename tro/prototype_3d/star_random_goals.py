"""Star tube with random goals: ours vs. a full sphere-decomposition baseline.

The arm starts threaded deep in the rounded star tube (the paper's scene: interactive setup q1) and leaves to 30
random goal configurations. Goals are sampled as in the Franka benchmark (80% of the joint ranges, hand above
0.1 m, at least 3 cm mesh clearance from the tube) and kept only if the unfiltered motion collides (an actual mesh
point with negative tube SDF).

Methods (same nominal, dt, gamma, input bounds, horizon, activation threshold, QP = summed.solve with DAQP):
  ours     vertex certificate with the unit lift, box = native 6-mm SDF patch, no refinement; the tube is the
           fitted tube field phi_B with its certified level.
  spheres  links AND tube decomposed into certified enclosing spheres; one CBF row per sphere pair,
           h = |x_a - c_b| - r_a - r_b. Link spheres: k-means centers of the remeshed link vertices (radius = max
           assigned distance + delta). Tube spheres: voxels of side s whose cube may meet the tube solid
           (exact SDF at the center <= s*sqrt(3)/2), clustered by k-means; radius = max(|v - c| + s*sqrt(3)/2).
Audit at every state: certified mesh-distance lower bound to the analytic tube (startube_setup.gap) and a collision
witness (negative analytic SDF at a mesh centroid).

    python prototype_3d/star_random_goals.py goals
    python prototype_3d/star_random_goals.py run ours [--trials 0 30] [--no-audit] [--tag T]
    python prototype_3d/star_random_goals.py run spheres --link-spheres 64 --tube-spheres 128 [...]
Writes results/star_random_goals/goals.json and results/star_random_goals/<method>[_tag]/trial_XX.{json,npz}.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('SUMMED_REFINE_DEPTH', '0')
import numpy as np                                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'interactive'))
import franka3d as F                                              # noqa: E402
import summed as SM                                               # noqa: E402
import summed3d as S3                                             # noqa: E402
import certificate_upgrades as CU                                 # noqa: E402
import continuous_baselines as CBL                                # noqa: E402
import server as SV                                               # noqa: E402
from startube_setup import gap as tube_gap                        # noqa: E402

ROOT = Path(HERE).parent
OUT = ROOT / 'results' / 'star_random_goals'
SETUP = Path(HERE) / 'interactive' / 'setups' / 'setup_20260928_171837.json'
CACHE_6MM = ROOT / 'results' / 'fine_native_6mm_trial' / 'cache_franka_6mm.npz'
FIELDS = {'6mm': CACHE_6MM, '8mm': ROOT / 'results' / 'native_patch_sizes' / '8mm' / 'cache_franka_8mm.npz', '12mm': None}
STEPS, N_GOALS, CLEAR, SEED = 1500, 30, .03, 0
UMAX7 = np.full(F.DOF, F.QD_MAX)


def scene(field='6mm'):
    d = json.load(open(SETUP))
    o = d['objects'][0]
    links, _, _ = F.build(cache_path=FIELDS[field]) if FIELDS[field] else F.build()
    assert all(abs(L['side'] - int(field[:-2]) / 1000) < 1e-12 for L in links)
    B = SV.build_objects()['star tube (rounded)']
    return np.array(d['q1'], float), links, B, np.array(o['R'], float), np.array(o['p'], float)


def witness(q, links, B, RB, pB):
    """Minimum analytic tube SDF over the link mesh centroids (negative: an actual mesh point inside the tube)."""
    T, _, _ = F.fk(q)
    w = np.inf
    for L in links:
        P = (L['gt'] @ T[L['frame']][:3, :3].T + T[L['frame']][:3, 3] - pB) @ RB
        lo, hi = np.asarray(B['blo']) - .01, np.asarray(B['bhi']) + .01
        m = np.all((P > lo) & (P < hi), axis=1)
        if m.any():
            w = min(w, float(B['sdf'](P[m]).min()))
    return w


class Audit:
    """startube_setup.gap (certified whole-mesh lower bound, capped at 5 cm) with per-link reuse: a link's previous
    bound minus a rigid-motion bound on its displacement (|dt| + |dR|_F r_L, tube frame) is still a lower bound; the
    link is evaluated again only when that transferred bound falls below the cap, so every bound below the cap is a
    direct evaluation, identical to tube_gap. The collision witness is evaluated only when the bound is <= 0."""
    CAP = .05

    def __init__(self, links, B, RB, pB, reuse=None):
        self.links, self.B, self.RB, self.pB = links, B, RB, pB
        self.reuse = self.CAP if reuse is None else reuse     # recompute a link when its transferred bound < reuse
        self.rL = [float(np.linalg.norm(L['gt'], axis=1).max() + L['gt_r'].max()) for L in links]
        self.prev = [None] * len(links)

    def __call__(self, q):
        T, _, _ = F.fk(q)
        B, RB, pB, g = self.B, self.RB, self.pB, self.CAP
        for j, L in enumerate(self.links):
            Rl, tl = RB.T @ T[L['frame']][:3, :3], RB.T @ (T[L['frame']][:3, 3] - pB)
            if self.prev[j] is not None:
                g0, R0, t0 = self.prev[j]
                moved = g0 - (np.linalg.norm(tl - t0) + np.linalg.norm(Rl - R0) * self.rL[j])
                if moved >= self.reuse:                      # still a certified lower bound for this link
                    g = min(g, moved)
                    continue
            P = L['gt'] @ Rl.T + tl
            outside = np.linalg.norm(np.maximum(np.maximum(B['blo'] - P, P - B['bhi']), 0.), axis=1)
            m = outside - L['gt_r'] < self.CAP
            gl = min(self.CAP, float((B['sdf'](P[m]) - L['gt_r'][m]).min())) if m.any() else self.CAP
            self.prev[j] = (gl, Rl, tl)
            g = min(g, gl)
        w = witness(q, self.links, B, RB, pB) if g <= 0 else np.inf
        return g, w


def sample_goal(rng, links, B, RB, pB):
    mid, half = (F.Q_MIN + F.Q_MAX) / 2, (F.Q_MAX - F.Q_MIN) / 2 * .8
    while True:
        q = rng.uniform(mid - half, mid + half)
        T, _, _ = F.fk(q)
        if T[8][2, 3] >= .1 and tube_gap(q, links, B, RB, pB) >= CLEAR:
            return q


def nominal(q, qg):
    return np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)


_W = {}


def _collides(qg):
    """Worker: does the unfiltered joint-space nominal from the start hit the tube?"""
    if not _W:
        _W['s'] = scene()
    q0, links, B, RB, pB = _W['s']
    q, wmin = q0.copy(), np.inf
    for _ in range(STEPS):
        wmin = min(wmin, witness(q, links, B, RB, pB))
        if wmin < 0 or np.linalg.norm(qg - q) < .05:
            break
        q = q + F.DT * nominal(q, qg)
    return bool(wmin < 0)


def make_goals(n_cand=150, workers=10):
    """Candidates in a fixed seeded order; their collision checks run in parallel and are saved as they finish
    (restartable); the goals are the first N_GOALS colliding candidates in candidate order."""
    from multiprocessing import Pool
    OUT.mkdir(parents=True, exist_ok=True)
    q0, links, B, RB, pB = scene()
    cand_path = OUT / 'goal_candidates.json'
    if cand_path.exists():
        cand = json.loads(cand_path.read_text())
    else:
        rng = np.random.default_rng(SEED)
        cand = dict(q=[sample_goal(rng, links, B, RB, pB).tolist() for _ in range(n_cand)], collides={})
        cand_path.write_text(json.dumps(cand))
    todo = [i for i in range(len(cand['q'])) if str(i) not in cand['collides']]
    with Pool(workers) as pool:
        for i, c in zip(todo, pool.imap(_collides, [np.array(cand['q'][i]) for i in todo])):
            cand['collides'][str(i)] = c
            cand_path.write_text(json.dumps(cand))
            print(f'candidate {i}: collides={c} ({sum(cand["collides"].values())} colliding so far)', flush=True)
    order = [i for i in range(len(cand['q'])) if cand['collides'][str(i)]]
    assert len(order) >= N_GOALS, f'only {len(order)} colliding candidates; increase n_cand'
    used = order[:N_GOALS]
    (OUT / 'goals.json').write_text(json.dumps(dict(
        q_start=q0.tolist(), goals=[cand['q'][i] for i in used], candidate_index=used,
        sampled=used[-1] + 1, seed=SEED, clearance_m=CLEAR,
        rule='80% joint ranges, hand z >= 0.1 m, mesh clearance >= 3 cm, unfiltered joint-space nominal collides '
             '(negative analytic tube SDF at a mesh centroid); first 30 colliding candidates in seeded order'),
        indent=1))


def tube_spheres(B, s, k):
    """Certified enclosing spheres of the tube solid (tube frame)."""
    from scipy.cluster.vq import kmeans2
    lo, hi = np.asarray(B['blo']) - s, np.asarray(B['bhi']) + s
    ax = [np.arange(lo[i], hi[i] + s, s) for i in range(3)]
    V = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
    half = s * np.sqrt(3) / 2
    V = V[B['sdf'](V) <= half]                         # every cube that may meet the solid
    C = kmeans2(V, k, seed=0, minit='++')[0]
    near = np.argmin(np.linalg.norm(V[:, None] - C[None], axis=2), axis=1)
    keep = np.unique(near)
    rad = np.array([np.linalg.norm(V[near == j] - C[j], axis=1).max() for j in keep]) + half
    return C[keep], rad


def sphere_rows(q, links, prim, Cw, rw):
    """Sphere-pair CBF rows n^T J_a u >= -gamma h for every pair with h < activation."""
    T, Z, O = F.fk(q)
    A, Cc, hmin = [], [], np.inf
    for L, (Pl, rl) in zip(links, prim):
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        X = Pl @ R.T + p
        D = X[:, None] - Cw[None]
        dist = np.linalg.norm(D, axis=2)
        h = dist - rl[:, None] - rw[None]
        hmin = min(hmin, float(h.min()))
        ia, ib = np.nonzero(h < F.ACT)
        if not len(ia):
            continue
        n = D[ia, ib] / dist[ia, ib][:, None]
        nj = L['n_joints']
        J = np.cross(Z[None, :nj], X[ia][:, None] - O[None, :nj])            # (rows, nj, 3)
        row = np.zeros((len(ia), F.DOF))
        row[:, :nj] = np.einsum('rjk,rk->rj', J, n)
        A.append(row); Cc.append(F.GAMMA * h[ia, ib])
    if not A:
        return np.zeros((0, F.DOF)), np.zeros(0), hmin
    return np.vstack(A), np.concatenate(Cc), hmin


def run(method, trials, audit, tag, n_link, n_tube, voxel, variant='unit', field='6mm'):
    q0, links, B, RB, pB = scene(field)
    goals = json.loads((OUT / 'goals.json').read_text())
    assert np.allclose(goals['q_start'], q0)
    SM.QP_SOLVER = 'daqp'
    info = dict(method=method, dt_s=F.DT, gamma=F.GAMMA, eta_m=F.ACT, steps=STEPS, qp='daqp, slack fallback')
    if method == 'ours':
        CU.configure(variant)                           # unit: vertex certificate; joint: joint Bernstein; w = 1
        for A in links:
            S3.prep_link(A)
        S3.prep_field(B, max(A['side'] * np.sqrt(3) / 2 for A in links))
        info.update(certificate=('vertex' if variant == 'unit' else 'joint Bernstein') + ', w = 1',
                    patch_mm=int(field[:-2]), tube_level_mm=1e3 * float(B['l']))
    else:
        prim = CBL.primitives(links, 'spheres_kmeans', .01, n_link)
        Ct, rt = tube_spheres(B, voxel, n_tube)
        Cw, rw = Ct @ RB.T + pB, rt
        info.update(link_spheres=int(sum(len(x[0]) for x in prim)), tube_spheres=len(Ct),
                    tube_radius_mm=[1e3 * float(rt.min()), 1e3 * float(rt.max())],
                    link_radius_mm=[1e3 * float(min(x[1].min() for x in prim)),
                                    1e3 * float(max(x[1].max() for x in prim))], voxel_mm=1e3 * voxel)
    out = OUT / (method + tag)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'setup.json').write_text(json.dumps(info, indent=1))
    for i in range(*trials):
        path = out / f'trial_{i:02d}.json'
        if path.exists():
            continue
        qg = np.array(goals['goals'][i])
        q, z = q0.copy(), None
        log = dict(q=[q.copy()], t=[], rows=[], gap=[], wit=[], h=[], slack=0, reached=None)
        for k in range(STEPS):
            u_nom = nominal(q, qg)
            t0 = time.perf_counter()
            if method == 'ours':
                T1, Z1, O1 = F.fk(q)
                res = [S3.pair_rows(A, T1[A['frame']][:3, :3], T1[A['frame']][:3, 3], B, RB, pB,
                                    [(Z1, O1, A['n_joints'], 1., 0)], F.ACT, UMAX7) for A in links]
                Ar, Tr, Cr, nb, hl, Er = S3.stack(res, F.DOF)[:6]
                u, s, z = SM.solve(u_nom, Ar, Tr, Cr, Er, F.BOX_G, F.BOX_H, z)
            else:
                Ar, Cr, hl = sphere_rows(q, links, prim, Cw, rw)
                u, s, z = SM.solve(u_nom, Ar, np.zeros((len(Cr), 1)), Cr, np.zeros((1, F.DOF)), F.BOX_G, F.BOX_H, z)
            log['t'].append(time.perf_counter() - t0)
            log['slack'] += s > 0
            log['rows'].append(len(Cr)); log['h'].append(hl)
            if audit:
                log['gap'].append(tube_gap(q, links, B, RB, pB)); log['wit'].append(witness(q, links, B, RB, pB))
            q = q + F.DT * u
            log['q'].append(q.copy())
            if np.linalg.norm(qg - q) < .05:
                log['reached'] = (k + 1) * F.DT
                break
        if audit:
            log['gap'].append(tube_gap(q, links, B, RB, pB)); log['wit'].append(witness(q, links, B, RB, pB))
        t = 1e3 * np.asarray(log['t'])
        m = dict(steps=len(t), reached_s=log['reached'], slack_steps=int(log['slack']),
                 h0_mm=1e3 * float(log['h'][0]), rows_max=int(max(log['rows'])),
                 t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)), t_max_ms=float(t.max()))
        if audit:
            m.update(min_gap_bound_mm=1e3 * float(min(log['gap'])), collision_witness_states=int(
                sum(w < 0 for w in log['wit'])), min_witness_sdf_mm=1e3 * float(min(log['wit'])))
        path.write_text(json.dumps(dict(trial=i, method=method, metrics=m), indent=1))
        np.savez_compressed(out / f'trial_{i:02d}.npz', q=np.asarray(log['q']), t=t, rows=np.asarray(log['rows']),
                            gap=np.asarray(log['gap']), h=np.asarray(log['h']))
        print(method + tag, i, {k_: (round(v, 2) if isinstance(v, float) else v) for k_, v in m.items()}, flush=True)


def _audit_one(args):
    """Worker: certified mesh-distance lower bounds of one saved trajectory (every state, in order)."""
    npz, reuse = args
    if not _W:
        _W['s'] = scene()
    q0, links, B, RB, pB = _W['s']
    A = Audit(links, B, RB, pB, reuse)
    gw = np.array([A(q) for q in np.load(npz)['q']])
    return npz, gw[:, 0], gw[:, 1]


def audit_folder(folder, reuse, workers):
    from multiprocessing import Pool
    d = OUT / folder
    todo = [str(p) for p in sorted(d.glob('trial_*.npz')) if not (d / p.name.replace('trial_', 'audit_')).exists()]
    with Pool(workers) as pool:
        for npz, g, w in pool.imap_unordered(_audit_one, [(p, reuse) for p in todo]):
            p = Path(npz)
            np.savez_compressed(p.parent / p.name.replace('trial_', 'audit_'), gap=g, witness=w)
            js = p.with_suffix('.json')
            rec = json.loads(js.read_text())
            rec['metrics'].update(min_gap_bound_mm=1e3 * float(g.min()), collision_witness_states=int((w < 0).sum()),
                                  min_witness_sdf_mm=1e3 * float(w.min()) if np.isfinite(w.min()) else None,
                                  audit=f'every state incl. final; certified lower bound, direct below '
                                        f'{1e3 * reuse:.0f} mm (link-bound reuse above), capped at 50 mm')
            js.write_text(json.dumps(rec, indent=1))
            print(folder, p.name, f'min gap {1e3 * g.min():.2f} mm, witness states {(w < 0).sum()}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cmd', choices=('goals', 'run', 'audit'))
    p.add_argument('--folder', default='')
    p.add_argument('--reuse-mm', type=float, default=20.)
    p.add_argument('--workers', type=int, default=12)
    p.add_argument('method', nargs='?', choices=('ours', 'spheres'))
    p.add_argument('--trials', nargs=2, type=int, default=[0, N_GOALS])
    p.add_argument('--no-audit', action='store_true')
    p.add_argument('--tag', default='')
    p.add_argument('--variant', choices=('unit', 'joint'), default='unit')
    p.add_argument('--field', choices=tuple(FIELDS), default='6mm')
    p.add_argument('--link-spheres', type=int, default=64)
    p.add_argument('--tube-spheres', type=int, default=128)
    p.add_argument('--voxel-mm', type=float, default=10.)
    a = p.parse_args()
    if a.cmd == 'goals':
        make_goals()
    elif a.cmd == 'audit':
        audit_folder(a.folder, a.reuse_mm / 1000, a.workers)
    else:
        run(a.method, a.trials, not a.no_audit, a.tag, a.link_spheres, a.tube_spheres, a.voxel_mm / 1000,
            a.variant, a.field)


if __name__ == '__main__':
    main()
