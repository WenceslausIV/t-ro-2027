"""Certificate upgrades on the fixed Franka trials (native SDF patches unless stated).

    python prototype_3d/certificate_upgrades.py --variant free --trials 0 30
Variants (all keep the fields, levels, covers, nominal controller, and trials of the paper):
  unit      vertex rows with the unit lift (the native-patch default of the paper)
  free      joint Bernstein rows with one multiplier w in [0, W_MAX] per active box, a QP variable
  sampled   'free' plus the sampled-data tightening of Theorem sampled (held inputs, dt = 10 ms)
Infeasible steps use a slack (as in the paper) unless --fallback zero, which applies u = 0.
Every saved state, including the final one, is audited with the mesh lower-bound routine of franka3d.
Results: results/certificate_upgrades/<variant>[_zero]/franka_XX.json; historical files are untouched.
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
import numpy as np

import franka3d as F
import summed as S
import summed3d as S3
import sampled_data as SD

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'certificate_upgrades'


def simulate(q0, qg, links, obst, variant, fallback, steps=1000):
    q = np.array(q0, float)
    log = dict(t=[], rows=[], boxes=[], h=[], gap=[], slack_steps=0, zero_steps=0, reached=None, q=[q.copy()],
               w_min=[], w_max=[], caps=[])
    z = None
    state = SD.State() if variant == 'sampled' else None
    for k in range(steps):
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter()
        if variant == 'sampled':
            A, T, C, nb, hl, E, Wc = SD.franka_rows(q, links, obst, state)
        else:
            res = S3.franka_rows(q, links, obst)
            A, T, C, nb, hl, E = res[:6]
            Wc = res[6] if len(res) == 7 else None
        if fallback == 'zero':
            u, feasible, z = SD.solve_or_stop(u_nom, A, T, C, E, F.BOX_G, F.BOX_H, z, Wc)
            log['zero_steps'] += not feasible
        else:
            u, s, z = S.solve(u_nom, A, T, C, E, F.BOX_G, F.BOX_H, z, Wc)
            log['slack_steps'] += s > 0
        log['t'].append(time.perf_counter() - t0)
        if variant == 'sampled':
            state.update(u, E, z)
            log['caps'].append(state.last_caps())
        if Wc is not None and Wc.shape[1] and z is not None and len(z) > len(u_nom) + len(E):
            w = 1 + z[len(u_nom) + len(E):] / np.sqrt(S.EPS_W)
            log['w_min'].append(float(w.min())); log['w_max'].append(float(w.max()))
        log['rows'].append(len(C)); log['boxes'].append(nb); log['h'].append(hl)
        log['gap'].append(F.real_gap(q, links, obst))
        q = q + F.DT * u
        log['q'].append(q.copy())
        if np.linalg.norm(qg - q) < .05:
            log['reached'] = (k + 1) * F.DT
            break
    log['gap'].append(F.real_gap(q, links, obst))           # final state
    return log


def summarize(log):
    t = 1e3 * np.asarray(log['t'])
    rows = np.asarray(log['rows'])
    gaps = np.asarray([g for g in log['gap'] if np.isfinite(g)] or [np.inf])
    return dict(steps=len(t), reached_s=log['reached'], slack_steps=int(log['slack_steps']),
                zero_steps=int(log['zero_steps']), min_gap_bound_mm=1e3 * float(gaps.min()),
                nonpositive_gap_states=int((gaps <= 0).sum()),
                min_h_mm=1e3 * float(np.min(log['h'])) if np.isfinite(np.min(log['h'])) else None,
                t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)),
                active_t_median_ms=float(np.median(t[rows > 0])) if (rows > 0).any() else None,
                active_t_p95_ms=float(np.percentile(t[rows > 0], 95)) if (rows > 0).any() else None,
                rows_max=int(rows.max()), boxes_max=int(max(log['boxes'])),
                w_range=[min(log['w_min']), max(log['w_max'])] if log['w_min'] else None)


def configure(variant):
    S.REFINE_DEPTH = 0
    if variant == 'unit':
        S.ROW_MODE, S.MULT_MODE = 'vertex', 'one'
    else:
        S.ROW_MODE, S.MULT_MODE = 'bernstein', 'free'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant', choices=('unit', 'free', 'sampled'), required=True)
    p.add_argument('--fallback', choices=('slack', 'zero'), default='slack')
    p.add_argument('--trials', nargs=2, type=int, default=[0, 30])
    p.add_argument('--reverse', action='store_true', help='run the trial range in descending order')
    p.add_argument('--tag', default='', help='suffix of the output folder (default: the QP solver)')
    p.add_argument('--qp', choices=('nnls', 'clarabel'), default='nnls',
                   help='QP solver: dense least-distance NNLS, or sparse interior point with exact acceptance test')
    p.add_argument('--field', choices=('12mm', '6mm'), default='12mm',
                   help='link SDF cells; 6mm uses the refitted fields of results/fine_native_6mm_trial')
    a = p.parse_args()
    configure(a.variant)
    S.QP_SOLVER = a.qp
    out = OUT / (a.variant + ('_zero' if a.fallback == 'zero' else '') + ('_6mm' if a.field == '6mm' else '')
                 + ('_' + a.tag if a.tag else ('_clarabel' if a.qp == 'clarabel' else '')))
    out.mkdir(parents=True, exist_ok=True)
    if a.field == '6mm':
        links, obst, info = F.build(cache_path=ROOT / 'results' / 'fine_native_6mm_trial' / 'cache_franka_6mm.npz')
        assert all(abs(L['side'] - .006) < 1e-12 for L in links)
    else:
        links, obst, info = F.build()
    S3.prep_all(links, obst.values())
    if a.variant == 'sampled':
        SD.prepare(links, obst)
    trials = json.loads((Path(F.HERE) / 'franka_trials.json').read_text())
    for i in (range(a.trials[1] - 1, a.trials[0] - 1, -1) if a.reverse else range(*a.trials)):
        path = out / f'franka_{i:02d}.json'
        if path.exists():
            continue
        q0, qg, _ = trials[i]
        log = simulate(q0, np.asarray(qg), links, obst, a.variant, a.fallback)
        rec = dict(trial=i, variant=a.variant, fallback=a.fallback, metrics=summarize(log),
                   cover=f'native SDF patches ({a.field}), no refinement', dt_s=F.DT, gamma=F.GAMMA, eta_m=F.ACT,
                   sampled=SD.describe() if a.variant == 'sampled' else None, qp_solver=a.qp,
                   qp_stats=dict(S.QP_STATS),
                   note='Saved-state mesh lower bounds include the final state. One thread; this machine.')
        path.write_text(json.dumps(rec, indent=1))
        np.savez_compressed(out / f'franka_{i:02d}.npz', q=np.asarray(log['q']), t=np.asarray(log['t']),
                            rows=np.asarray(log['rows']), gap=np.asarray(log['gap']), h=np.asarray(log['h']))
        m = rec['metrics']
        print(i, a.variant, a.fallback, {k: m[k] for k in ('reached_s', 'slack_steps', 'zero_steps',
                                                          'min_gap_bound_mm', 't_median_ms', 'rows_max')},
              flush=True)


if __name__ == '__main__':
    main()
