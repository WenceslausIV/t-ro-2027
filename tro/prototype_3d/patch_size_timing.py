"""Filter timing of the sampled-data certificate for several native patch sizes, one process, nothing else running.

    python prototype_3d/patch_size_timing.py [--fields 6mm 12mm 18mm 24mm] [--trials 0 30]
Closed loop as in certificate_upgrades.py (--variant sampled --fallback zero --qp daqp), without the mesh audit,
so the trajectories are those of results/certificate_upgrades/sampled_zero_<field>_fixed. Box = native patch.
Writes results/patch_size_timing/<field>.npz (per-step times) and summary.json.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('SUMMED_REFINE_DEPTH', '0')
import argparse
import json
import platform
from pathlib import Path

import numpy as np

import certificate_upgrades as CU
import franka3d as F
import sampled_data as SD
import summed as S
import summed3d as S3

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'patch_size_timing'
CACHES = {'6mm': ROOT / 'results' / 'fine_native_6mm_trial' / 'cache_franka_6mm.npz', '12mm': None,
          '4mm': ROOT / 'results' / 'native_patch_sizes' / '4mm' / 'cache_franka_4mm.npz',
          '5mm': ROOT / 'results' / 'native_patch_sizes' / '5mm' / 'cache_franka_5mm.npz',
          '8mm': ROOT / 'results' / 'native_patch_sizes' / '8mm' / 'cache_franka_8mm.npz',
          '18mm': ROOT / 'results' / 'native_patch_sizes' / '18mm' / 'cache_franka_18mm.npz',
          '24mm': ROOT / 'results' / 'native_patch_sizes' / '24mm' / 'cache_franka_24mm.npz'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fields', nargs='+', default=list(CACHES))
    p.add_argument('--trials', nargs=2, type=int, default=[0, 30])
    p.add_argument('--variant', choices=('sampled', 'free', 'joint', 'unit'), default='sampled',
                   help='free: rows at the samples only (no sampled-data certificate), like the baselines')
    a = p.parse_args()
    fallback = 'zero' if a.variant == 'sampled' else 'slack'
    refine = int(os.environ.get('SUMMED_REFINE_DEPTH', '0'))
    suffix = ('' if a.variant == 'sampled' else '_' + a.variant) + (f'_refine{refine}' if refine else '')
    OUT.mkdir(parents=True, exist_ok=True)
    CU.configure(a.variant)
    S.QP_SOLVER = 'daqp'
    F.real_gap = lambda *args: 1.0                       # audit excluded (it is not part of the filter)
    trials = json.loads((Path(F.HERE) / 'franka_trials.json').read_text())
    summary = json.loads((OUT / 'summary.json').read_text()) if (OUT / 'summary.json').exists() else {}
    for field in a.fields:
        c = CACHES[field]
        links, obst, info = F.build(cache_path=c) if c else F.build()
        S3.prep_all(links, obst.values())
        if a.variant == 'sampled':
            SD.prepare(links, obst)
        t_all, per = [], {}
        for i in range(*a.trials):
            q0, qg, _ = trials[i]
            log = CU.simulate(q0, np.asarray(qg), links, obst, a.variant, fallback)
            t = 1e3 * np.asarray(log['t'])
            run = (f'{a.variant}{"_zero" if fallback == "zero" else ""}{"" if field == "12mm" else "_" + field}_fixed'
                   + (f'_refine{refine}' if refine else ''))
            ref = ROOT / 'results' / 'certificate_upgrades' / run / f'franka_{i:02d}.npz'
            dq = float(np.abs(np.load(ref)['q'] - np.asarray(log['q'])).max()) if ref.exists() and len(np.load(ref)['q']) == len(log['q']) else None
            per[i] = dict(steps=len(t), median_ms=float(np.median(t)), p95_ms=float(np.percentile(t, 95)),
                          max_ms=float(t.max()), zero_steps=int(log.get('zero_steps', 0)), max_dq_vs_run=dq)
            t_all.append(t)
            print(field, i, per[i], flush=True)
        t = np.concatenate(t_all)
        np.savez_compressed(OUT / f'{field}{suffix}.npz', t_ms=t, trial_lengths=[len(x) for x in t_all])
        summary[field + suffix] = dict(steps=len(t), median_ms=float(np.median(t)), p90_ms=float(np.percentile(t, 90)),
                              p95_ms=float(np.percentile(t, 95)), p99_ms=float(np.percentile(t, 99)),
                              max_ms=float(t.max()), mean_ms=float(t.mean()),
                              over_10ms_percent=float(100 * np.mean(t > 10)), trials=per,
                              variant=a.variant, fallback=fallback, refine_depth=refine,
                              machine=platform.processor() or platform.machine(),
                              note='one process, one thread, no other experiment running; audit excluded')
        (OUT / 'summary.json').write_text(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
