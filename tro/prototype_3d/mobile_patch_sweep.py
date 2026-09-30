"""One-versus-one mobile robots: native patch size sweep and disk baselines on the same crossing scenarios.

    python prototype_3d/mobile_patch_sweep.py [--scenarios 20] [--methods patch10 patch20 patch40 patch80 disk1 disk40 disk20]
Two non-convex robots (the shapes of two_robot_cover_comparison: seed 3, shapes 0 and 3), SE(2) single integrators,
|v_x|, |v_y| <= 1 m/s, |omega| <= 2 rad/s, dt = 10 ms, 12-s horizon, gamma = 5, the same nominal controller.
Each scenario swaps the robots across a circle of radius 1.5 m on crossing lines, with random headings.

  patch<h>  link SDF refitted with h-mm cells (training spacing 5 mm); B's cover = its native cells (box = patch);
            joint Bernstein rows with one optimized multiplier per active box, activation 5 cm (Sec. native)
  disk1     each robot as its circumscribed disk (classic mobile-robot CBF)
  disk<c>   each robot as the union of disks circumscribing the c-mm grid cells that meet it (multi-disk CBF);
            rows for disk pairs within 30 cm
All methods use the same QP (summed.solve, DAQP constraint generation, exact check of every input).
Audit: exact polygon distance at every saved state (vertex-to-edge distances plus crossing tests), final included.
Rows are enforced at the 10-ms samples (no sampled-data tightening in the plane).
Results: results/mobile_patch_sweep/<method>.json and summary.json.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'planar'))
import numpy as np
import torch
from matplotlib.path import Path as MPath

import cspace_experiments as E
from cspace_cbf_5robots import random_shape
from dock_cover import body_field, cover_2d, tight_level_2d, make_summed_fn
import franka3d as F
import summed as SM

OUT = ROOT / 'results' / 'mobile_patch_sweep'
DT, HORIZON, GAMMA, RADIUS = .01, 12., 5., 1.5
J2 = np.array([[0., -1.], [1., 0.]])


def shapes():
    rng = np.random.default_rng(3)
    gt = [random_shape(rng) for _ in range(5)]
    return [gt[0], gt[3]]


def scenarios(n, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        th = rng.uniform(0, 2 * np.pi)
        tb = th + np.pi / 2 + rng.uniform(-.4, .4)
        pa, pb = RADIUS * np.array([np.cos(th), np.sin(th)]), RADIUS * np.array([np.cos(tb), np.sin(tb)])
        h = rng.uniform(-np.pi, np.pi, 4)
        out.append((np.array([[*pa, h[0]], [*pb, h[1]]]), np.array([[*-pa, h[2]], [*-pb, h[3]]])))
    return out


def world(GT, x):
    return [g @ E.rot(p[2]).T + p[:2] for g, p in zip(GT, x)]


def disk_cover(g, cell):
    """Disks circumscribing the grid cells that meet the polygon g (body frame); their union contains g."""
    if cell is None:
        return np.zeros((1, 2)), np.array([np.linalg.norm(g, axis=1).max()])
    lo, hi = g.min(axis=0) - cell, g.max(axis=0) + cell
    X, Y = np.meshgrid(np.arange(lo[0], hi[0], cell) + cell / 2, np.arange(lo[1], hi[1], cell) + cell / 2)
    c = np.c_[X.ravel(), Y.ravel()]
    rc = cell / np.sqrt(2)
    inside = MPath(g).contains_points(c)
    A, B = g, np.roll(g, -1, axis=0)                     # exact distance from the centers to the boundary
    AB = B - A
    t = np.clip(((c[:, None] - A[None]) * AB[None]).sum(-1) / (AB ** 2).sum(-1)[None], 0, 1)
    d = np.linalg.norm(c[:, None] - (A[None] + t[..., None] * AB[None]), axis=-1).min(axis=1)
    keep = inside | (d <= rc)
    return c[keep], np.full(int(keep.sum()), rc)


def disk_fn(cov, act=.3):
    """CBF rows of all disk pairs within act: h = |d| - ri - rj, d = pA + R_A ci - pB - R_B cj."""
    (ca, ra), (cb, rb) = cov

    def fn(xa, xb):
        wa, wb = ca @ E.rot(xa[2]).T, cb @ E.rot(xb[2]).T
        d = (xa[:2] + wa)[:, None] - (xb[:2] + wb)[None]
        nd = np.linalg.norm(d, axis=-1)
        h = nd - ra[:, None] - rb[None]
        i, j = np.nonzero(h < act)
        if not len(i):
            return None
        n = d[i, j] / nd[i, j, None]
        A = np.c_[n, (n * (wa[i] @ J2.T)).sum(1), -n, -(n * (wb[j] @ J2.T)).sum(1)]
        return A, np.zeros((len(i), 1)), GAMMA * h[i, j], len(i), float(h[i, j].min()), np.zeros((1, 6))
    return fn


def simulate(GT, fn, start, goal):
    x = start.copy()
    G, h = E.box_rows(2, 'si', 1., 2.)
    z, arrival, t_ms, rows, slack, states, umod = None, None, [], [], 0, [x.copy()], 0.
    for k in range(round(HORIZON / DT)):
        tic = time.perf_counter()
        d = goal[:, :2] - x[:, :2]
        v = 1.5 * d
        v /= np.maximum(1., np.linalg.norm(v, axis=1))[:, None]
        nom = np.c_[v, np.clip(2. * E.wrap(goal[:, 2] - x[:, 2]), -2., 2.)].ravel()
        res = fn(x[0], x[1])
        if res is None:
            A, T, C, twist, Wc = np.zeros((0, 6)), np.zeros((0, 0)), np.zeros(0), np.zeros((0, 6)), None
        else:
            A, T, C, _, _, twist = res[:6]
            Wc = res[6] if len(res) == 7 else None
        u, s, z = SM.solve(nom, A, T, C, twist, G, h, z, Wc)
        t_ms.append(1e3 * (time.perf_counter() - tic)); rows.append(len(C)); slack += s > 0
        umod += DT * np.linalg.norm(u - nom)
        x = x + DT * u.reshape(2, 3)
        states.append(x.copy())
        if np.all(np.linalg.norm(x[:, :2] - goal[:, :2], axis=1) < .05) and np.all(np.abs(E.wrap(x[:, 2] - goal[:, 2])) < .05):
            arrival = (k + 1) * DT
            break
    gaps = np.array([E.C.true_distance(*world(GT, s)) for s in states])
    t = np.asarray(t_ms)
    return dict(arrival_s=arrival, collision_states=int((gaps <= 0).sum()), min_gap_mm=1e3 * float(gaps.min()),
                slack_steps=int(slack), steps=len(t), max_rows=int(max(rows)), input_modification=float(umod),
                t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)), t_max_ms=float(t.max()),
                t_steps_ms=t.tolist())


def method_fn(name, GT, info):
    if name.startswith('patch'):
        h = int(name[5:]) / 1000
        SM.ROW_MODE, SM.MULT_MODE = 'bernstein', 'free'
        fields, majors, levels = [], [], []
        for g in GT:
            f = body_field(g, h, .1, training_spacing=.005)
            m = F.hessian_majorant(f)
            fields.append(f); majors.append(m); levels.append(tight_level_2d(f, m, g))
        centers, side = cover_2d(fields[1], levels[1])            # native cells of B: box = patch
        assert abs(side - h) < 1e-12
        body = SM.prepare_body(fields[1], levels[1], centers, side, 2)
        info.update(cell_mm=1e3 * h, levels_mm=[1e3 * l for l in levels], patches=len(centers))
        return make_summed_fn(fields[0], majors[0], levels[0], body)
    cell = None if name == 'disk1' else int(name[4:]) / 1000
    cov = [disk_cover(g, cell) for g in GT]
    info.update(disks=[len(c[0]) for c in cov], disk_radius_mm=[1e3 * float(c[1].max()) for c in cov])
    return disk_fn(cov)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenarios', type=int, default=20)
    p.add_argument('--methods', nargs='+', default=['patch10', 'patch20', 'patch40', 'patch80', 'disk1', 'disk40', 'disk20'])
    a = p.parse_args()
    torch.set_num_threads(1)
    SM.REFINE_DEPTH, SM.QP_SOLVER = 0, 'daqp'
    OUT.mkdir(parents=True, exist_ok=True)
    GT = shapes()
    sc = scenarios(a.scenarios)
    for name in a.methods:
        path = OUT / f'{name}.json'
        if path.exists():
            continue
        info = {}
        t0 = time.perf_counter()
        fn = method_fn(name, GT, info)
        info['setup_s'] = time.perf_counter() - t0
        runs = [dict(scenario=i, **simulate(GT, fn, s, g)) for i, (s, g) in enumerate(sc)]
        t = np.concatenate([r.pop('t_steps_ms') for r in runs])
        summ = dict(method=name, **info, scenarios=len(runs),
                    both_goals=sum(r['arrival_s'] is not None for r in runs),
                    mean_arrival_s=float(np.mean([r['arrival_s'] for r in runs if r['arrival_s'] is not None]))
                    if any(r['arrival_s'] is not None for r in runs) else None,
                    collision_runs=sum(r['collision_states'] > 0 for r in runs),
                    min_gap_mm=min(r['min_gap_mm'] for r in runs),
                    median_min_gap_mm=float(np.median([r['min_gap_mm'] for r in runs])),
                    slack_runs=sum(r['slack_steps'] > 0 for r in runs),
                    mean_input_modification=float(np.mean([r['input_modification'] for r in runs])),
                    t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)),
                    t_p99_ms=float(np.percentile(t, 99)), t_max_ms=float(t.max()),
                    over_10ms_percent=float(100 * np.mean(t > 10)), max_rows=max(r['max_rows'] for r in runs))
        path.write_text(json.dumps(dict(summary=summ, runs=runs), indent=1))
        print(json.dumps(summ), flush=True)
    rows = {}
    for q in sorted(OUT.glob('*.json')):
        if q.name != 'summary.json':
            rows[q.stem] = json.loads(q.read_text())['summary']
    (OUT / 'summary.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
