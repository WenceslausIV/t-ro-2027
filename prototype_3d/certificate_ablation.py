"""Reproducible certificate comparisons; never overwrite the paper's result files.

Examples:
  python prototype_3d/certificate_ablation.py --scenes dock star --variants vertex-one bernstein-one vertex-proj
  python prototype_3d/certificate_ablation.py --scenes franka --variants vertex-one bernstein-one
All variants use the same geometry and unit-lift screening, with distance-independent
refinement (near=inf) unless --near is supplied. Timings are single-thread local measurements.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time

for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import summed as S
import summed3d as S3
import franka3d as F

ROOT = Path(__file__).resolve().parent.parent


def digest(arrays):
    h = hashlib.sha256()
    for x in arrays:
        a = np.ascontiguousarray(x)
        h.update(str((a.shape, a.dtype)).encode()); h.update(a.tobytes())
    return h.hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def summary(L, dt, planar=False):
    ts = np.asarray(L['t_ctrl' if planar else 't']) * 1e3
    u, un = np.asarray(L['u_t' if planar else 'u']), np.asarray(L['u_nom_t' if planar else 'u_nom'])
    gaps = np.asarray(L['d_t' if planar else 'gap'])
    hs = np.asarray(L['h_t' if planar else 'h'])
    return dict(steps=len(ts), reached_s=L['reached'],
                min_gap_bound_mm=1e3 * (L['min_gt'] if planar else gaps.min()),
                distance_scope='polygon intervals and states' if planar else 'evaluated states; existing mesh lower-bound routine',
                min_h_mm=1e3 * (np.nanmin(hs) if np.isfinite(hs).any() else np.inf),
                slack_steps=L['slack_steps' if planar else 'slack'],
                uncertified_gap_states=int(np.sum(gaps < 0)),
                control_change_integral=float(dt * np.linalg.norm(u - un, axis=1).sum()),
                t_median_ms=float(np.median(ts)), t_p95_ms=float(np.percentile(ts, 95)),
                rows_max=int(max(L['n_rows' if planar else 'rows'], default=0)),
                boxes_max=int(max(L['n_box' if planar else 'boxes'], default=0)))


def save(path, L, meta, planar=False):
    result = dict(metadata=meta, metrics=summary(L, .01, planar))
    if planar:
        result['metrics']['final_pose'] = L['final'].tolist()
    arr = dict(q=np.asarray(L['traj' if planar else 'q']),
               u=np.asarray(L['u_t' if planar else 'u']),
               u_nom=np.asarray(L['u_nom_t' if planar else 'u_nom']),
               gap=np.asarray(L['d_t' if planar else 'gap']),
               time_s=np.asarray(L['t_ctrl' if planar else 't']))
    np.savez_compressed(path.with_suffix('.npz'), **arr)
    path.write_text(json.dumps(clean(result), indent=2, allow_nan=False), encoding='utf-8')
    print(path.name, json.dumps(clean(result['metrics'])), flush=True)


def settings(args, variant, geometry_hash, scene):
    S.ROW_MODE, S.MULT_MODE = variant.split('-')
    S.REFINE_NEAR = args.near
    S.REFINE_DEPTH = args.depth
    return dict(scene=scene, geometry_sha256=geometry_hash, row_mode=S.ROW_MODE, multiplier=S.MULT_MODE,
                refine_near_m='inf' if np.isinf(args.near) else args.near,
                refine_theta=S.REFINE_THETA, refine_depth=S.REFINE_DEPTH, dt_s=.01,
                python=platform.python_version(), platform=platform.platform(), threads=1,
                certificate_sha256=hashlib.sha256((ROOT / 'prototype_3d/summed.py').read_bytes()).hexdigest())


def dock(args):
    import dock_cover as D
    GS, GC = D.station_body(), D.craft_body()
    fs, fc = D.body_field(GS, .0125, .15), D.body_field(GC, .01, .06)
    ms, mc = F.hessian_majorant(fs), F.hessian_majorant(fc)
    ls, lc = D.tight_level_2d(fs, ms, GS), D.tight_level_2d(fc, mc, GC)
    pc, side = D.cover_2d(fc, lc, .005)
    body = S.prepare_body(fc, lc, pc, side, 2)
    fn = D.make_summed_fn(fs, ms, ls, body)
    sha = digest([fs.W, fc.W, [ls, lc], pc, GS, GC])
    shapes = [type('Shape', (), dict(rho=float(np.linalg.norm(g, axis=1).max())))() for g in (GS, GC)]
    for variant in args.variants:
        path = args.out / f'dock_{variant}.json'
        if path.exists():
            print('skip existing', path.name, flush=True); continue
        meta = settings(args, variant, sha, 'dock')
        print('running dock', variant, flush=True)
        L = D.simulate_team([GS, GC], shapes, {(0, 1): fn}, D.DOCK2_START, D.DOCK2_GOAL, 'summed',
                            steps=args.dock_steps, v_max=np.array([0., 1.]), w_max=np.array([0., 2.]), record=True)
        save(path, L, meta, planar=True)


def star(args):
    import startube_setup as ST
    import trimesh
    links, _ = ST.S.D.build(cover_side='cached')
    for A in links:
        vs, faces = trimesh.remesh.subdivide_to_size(A['V'], A['F'], max_edge=.008)
        triangles = vs[faces]
        A['gt8'] = triangles.mean(axis=1)
        A['gt8_r'] = np.linalg.norm(triangles - A['gt8'][:, None], axis=2).max(axis=1)
        S3.prep_link(A)
    B = ST.S.build_objects()['star tube (rounded)']
    S3.prep_field(B)
    setup = ROOT / 'prototype_3d/interactive/setups/setup_20260928_171837.json'
    _, q0, pB, RB = ST.load(str(setup))
    Th, _, _ = F.fk(ST.S.Q_HOME)
    gp, gr = Th[8][:3, :3] @ ST.S.TCP + Th[8][:3, 3], Th[8][:3, :3]
    arrays = [B['fs'].W, [B['l']], q0, pB, RB]
    for A in links:
        arrays += [A['sb']['f'].W, [A['l']], A['PC'], A['V'], A['F']]
    sha = digest(arrays)
    for variant in args.variants:
        path = args.out / f'star_{variant}.json'
        if path.exists():
            print('skip existing', path.name, flush=True); continue
        meta = settings(args, variant, sha, 'rounded star tube')
        meta['obstacle_level_mm'] = 1e3 * B['l']
        print('running rounded star', variant, flush=True)
        L = ST.simulate(q0, gp, gr, links, B, RB, pB, 'ours', steps=args.star_steps, record_controls=True)
        save(path, L, meta)


def franka(args):
    links, obst, _ = F.build(cover_side='cached')
    S3.prep_all(links, obst.values())
    old = json.loads((ROOT / 'prototype_3d/franka_results.json').read_text())['trials']
    if args.trials:
        ids = args.trials
    else:
        ids = ([r['trial'] for r in old if r['reached'] is not None][:3]
               + [r['trial'] for r in old if r['reached'] is None][:3])
    trials = json.loads((ROOT / 'prototype_3d/franka_trials.json').read_text())
    arrays = []
    for A in links:
        arrays += [A['sb']['f'].W, [A['l']], A['PC']]
    for B in obst.values():
        arrays += [B['f'].W, [B['l']]]
    sha = digest(arrays)
    print('Franka fixed trial IDs (zero-based):', ids, flush=True)
    for variant in args.variants:
        for idx in ids:
            path = args.out / f'franka_{idx:02d}_{variant}.json'
            if path.exists():
                print('skip existing', path.name, flush=True); continue
            meta = settings(args, variant, sha, 'franka')
            meta['trial_zero_based'] = idx
            q0, qg, _ = trials[idx]
            print('running franka', idx, variant, flush=True)
            L = F.simulate(q0, np.asarray(qg), links, obst, 'ours', steps=args.franka_steps, record=True)
            save(path, L, meta)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenes', nargs='+', choices=('dock', 'star', 'franka'), default=['dock', 'star'])
    p.add_argument('--variants', nargs='+', default=['vertex-one', 'bernstein-one', 'vertex-norm', 'bernstein-norm', 'vertex-proj', 'bernstein-proj'])
    p.add_argument('--near', type=float, default=float('inf'))
    p.add_argument('--depth', type=int, default=2, help='Optional refinement depth for this historical ablation')
    p.add_argument('--out', type=Path, default=ROOT / 'results/certificate_ablation/pilot')
    p.add_argument('--dock-steps', type=int, default=900)
    p.add_argument('--star-steps', type=int, default=600)
    p.add_argument('--franka-steps', type=int, default=1000)
    p.add_argument('--trials', nargs='+', type=int)
    args = p.parse_args()
    for v in args.variants:
        row, mult = v.split('-')
        if row not in S.ROW_MODES or mult not in S.MULT_MODES:
            p.error(f'unknown variant {v}')
    args.out.mkdir(parents=True, exist_ok=True)
    for scene in args.scenes:
        globals()[scene](args)


if __name__ == '__main__':
    main()
