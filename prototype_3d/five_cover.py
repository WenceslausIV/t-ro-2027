"""
Five-robot experiments of Sec. VIII-D with the SURFACE-COVER construction (dock_cover.py) next to the
tabulated configuration-space fields, on identical instances, simulator, nominal controller, and QP.

One workspace field per body with the tightest level certified on its ground-truth polygon
(Proposition 1), a Bernstein cover of its certified curve, and for every pair (i, j) the corner barriers
of j's cover against i's field (one direction per pair suffices, Theorem 2).

    python prototype_3d/five_cover.py stats [side_mm]   # 20 random instances of cspace_experiments.exp_stats
    python prototype_3d/five_cover.py swap  [side_mm]   # perturbed antipodal swap of Fig. 5 (2x shapes)

Results: ../results/five_cover_stats.json, ../results/five_cover_swap.{json,npz}
"""
import itertools
import json
import os
import sys
import time

import numpy as np

import franka3d as F
from dock_cover import body_field, cover_2d, make_cover_fn, make_summed_fn, tight_level_2d
import summed as SM

sys.path.insert(0, os.path.dirname(F.HERE))
import cspace_experiments as E                                          # noqa: E402

RES = os.path.join(os.path.dirname(F.HERE), 'results')


def cover_setup(GT, side, h, margin=.10):
    """Fields, certified levels, and covers of all bodies; counted corner-barrier functions of all pairs."""
    t0 = time.perf_counter()
    bodies = []
    for G in GT:
        f = body_field(G, h, margin)
        M = F.hessian_majorant(f)
        l = tight_level_2d(f, M, G)
        PC, sz = cover_2d(f, l, side)
        bodies.append((f, M, l, PC, sz))
    counts = []

    def counted(fn):
        def g(xi, xj):
            R, hv = fn(xi, xj)
            counts.append(len(hv))
            return R, hv
        return g

    fns = {(i, j): counted(make_cover_fn(*bodies[i][:3], *bodies[j][3:]))
           for i, j in itertools.combinations(range(len(GT)), 2)}
    t_cover = time.perf_counter() - t0
    # summed-field barriers: Taylor data on every box, third-derivative majorants, shared gradient bounds
    t1 = time.perf_counter()
    sb = [SM.prepare_body(f, l, PC, sz, 2) for f, M, l, PC, sz in bodies]
    reach = max(float(F.clusters(b['PC'], b['r'])[2].max()) for b in sb)
    Gn = [F.neighborhood_bounds(b[0], reach) for b in bodies]
    sfns = {(i, j): make_summed_fn(*bodies[i][:3], sb[j], bounds=Gn[i])
            for i, j in itertools.combinations(range(len(GT)), 2)}
    info = dict(levels_mm=[1e3 * b[2] for b in bodies], n_boxes=[len(b[3]) for b in bodies],
                side_mm=1e3 * bodies[0][4], cell_mm=1e3 * h, t_setup=t_cover,
                t_setup_summed=t_cover + time.perf_counter() - t1)
    return fns, counts, bodies, info, sfns


def max_barriers(counts, n_pairs):
    c = np.array(counts)
    return int(c[:len(c) // n_pairs * n_pairs].reshape(-1, n_pairs).sum(axis=1).max()) if len(c) else 0


def summarize(rec, n):
    reached = [r['reached'] for r in rec if r['reached'] is not None]
    return dict(collisions=int(sum(r['min_gt'] <= 0 for r in rec)), min_gt=min(r['min_gt'] for r in rec),
                min_h=min(r['min_h'] for r in rec), success=len(reached),
                mean_time=float(np.mean(reached)) if reached else None,
                t_med=float(np.median([r['t_med'] for r in rec])), t_max=float(max(r['t_max'] for r in rec)),
                slack_steps=int(sum(r['slack_steps'] for r in rec)), tv_med=float(np.median([r['tv'] for r in rec])),
                barriers_max=max(r.get('barriers_max', 0) for r in rec),
                barriers_med=float(np.median([r.get('barriers_max', 0) for r in rec])), n=n)


def exp_stats(side, n_trials=20, steps=1500):
    GT, shapes, fields, _ = E.load_cache()
    rho = np.array([s.rho for s in shapes])
    rng = np.random.default_rng(11)                        # identical instances to cspace_experiments.exp_stats
    trials = [(E.sample_poses(rng, rho, 3.0, 0.5), E.sample_poses(rng, rho, 3.0, 0.5)) for _ in range(n_trials)]
    fns, counts, _, info, sfns = cover_setup(GT, side, h=.01)
    print('cover setup', info, flush=True)
    out = dict(cover_setup=info)
    methods = [m for m in ('cspace', 'summed', 'cover') if m in (sys.argv[3:] or ['cspace', 'summed'])]
    for method, flds in ((m, {'cspace': fields, 'cover': fns, 'summed': sfns}[m]) for m in methods):
        rec, pooled = [], []
        t0 = time.perf_counter()
        for k, (s, g) in enumerate(trials):
            counts.clear()
            L = E.simulate_team(GT, shapes, flds, s, g, method, steps=steps)
            r = dict(min_gt=L['min_gt'], min_h=L['min_h'], reached=L['reached'], path=L['path'].sum(),
                     t_med=np.median(L['t_ctrl']) * 1e3, t_max=L['t_ctrl'].max() * 1e3,
                     t_p95=np.percentile(L['t_ctrl'], 95) * 1e3,
                     slack_steps=L['slack_steps'], tv=L['tv'])
            if method == 'cover':
                r['barriers_max'] = max_barriers(counts, len(fns))
            if method == 'summed':
                r['barriers_max'] = int(max(L['n_rows'])); r['boxes_max'] = int(max(L['n_box']))
            rec.append(r)
            pooled.append(L['t_ctrl'])
            print(f"  [{method}] trial {k:2d}: min GT {1e3 * r['min_gt']:+7.1f} mm  min h {1e3 * r['min_h']:+6.1f} mm  "
                  f"reached {r['reached']}  t_med {r['t_med']:.2f} ms  bars {r.get('barriers_max', '-')}", flush=True)
        out[method] = dict(summarize(rec, n_trials), trials=rec, wall=time.perf_counter() - t0,
                           t_med_pooled=1e3 * float(np.median(np.concatenate(pooled))),
                           t_p95_pooled=1e3 * float(np.percentile(np.concatenate(pooled), 95)))
        print(method, {k: v for k, v in out[method].items() if k != 'trials'}, flush=True)
    both = [k for k in range(n_trials) if all(out[m]['trials'][k]['reached'] is not None for m in methods)]
    out['common_trials'] = both
    for m in methods:
        out[m]['common_mean_time'] = float(np.mean([out[m]['trials'][k]['reached'] for k in both])) if both else None
        out[m]['common_mean_path'] = float(np.mean([out[m]['trials'][k]['path'] for k in both])) if both else None
    E.dump(f'five_{"_".join(methods)}_stats_{1e3 * side:g}', out)


def exp_swap(side):
    from cspace_cbf_5robots import five_robot_setup
    cfg = json.load(open(os.path.join(RES, 'five_jitter_si.json')))
    GT, shapes, fields, starts, goals, _ = five_robot_setup(
        scale=cfg['scale'], radius=cfg['radius'], jitter_seed=cfg['chosen']['jitter_seed'],
        jitter=tuple(cfg['jitter']), goal_shift=cfg.get('goal_shift'), n_robots=cfg.get('n_robots', 5))
    fns, counts, bodies, info, sfns = cover_setup(GT, side, h=.01)
    print('cover setup', info, flush=True)
    out = dict(cover_setup=info)
    methods = [m for m in ('cspace', 'summed', 'cover') if m in (sys.argv[3:] or ['cspace', 'summed'])]
    for method, flds in ((m, {'cspace': fields, 'cover': fns, 'summed': sfns}[m]) for m in methods):
        counts.clear()
        L = E.simulate_team(GT, shapes, flds, starts, goals, method, steps=cfg['steps'], record=True)
        nb = (np.array(counts).reshape(-1, len(fns)).sum(axis=1) if method == 'cover'
              else np.array(L['n_rows']) if method == 'summed' else np.zeros(len(L['traj']), int))
        out[method] = dict(reached=L['reached'], min_gt_mm=1e3 * L['min_gt'], min_h_mm=1e3 * float(np.nanmin(L['h_t'])),
                           t_med_ms=1e3 * float(np.median(L['t_ctrl'])), t_max_ms=1e3 * float(L['t_ctrl'].max()),
                           t_p95_ms=1e3 * float(np.percentile(L['t_ctrl'], 95)),
                           tv=L['tv'], slack=L['slack_steps'], barriers_max=int(nb.max()))
        print(method, out[method], flush=True)
        np.savez_compressed(os.path.join(RES, f'five_cover_swap_{method}_{1e3 * side:g}.npz'), traj=np.array(L['traj']),
                            h=np.array(L['h_t']), d=np.array(L['d_t']), nb=nb, t=L['t_ctrl'], goals=goals,
                            nbox=np.array(L.get('n_box', [])),
                            **({f'PC_{i}': b[3] for i, b in enumerate(bodies)} if method != 'cspace' else {}))
    E.dump(f'five_{"_".join(methods)}_swap_{1e3 * side:g}', out)


if __name__ == '__main__':
    what = sys.argv[1]
    side = float(sys.argv[2]) / 1e3 if len(sys.argv) > 2 else .01
    {'stats': exp_stats, 'swap': exp_swap}[what](side)
