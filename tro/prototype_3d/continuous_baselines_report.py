"""Step 4 of results/continuous_baselines/PLAN.md: one table for the baselines and ours.

    python prototype_3d/continuous_baselines_report.py
Outcomes from results/continuous_baselines/<method>/ and results/certificate_upgrades/<ours>/; times from the
isolated runs (<method>_timing/ and results/patch_size_timing/summary.json) when present, else from the outcome runs.
Writes results/continuous_baselines/summary.json and README.md.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CB = ROOT / 'results' / 'continuous_baselines'
CU = ROOT / 'results' / 'certificate_upgrades'
PT = ROOT / 'results' / 'patch_size_timing' / 'summary.json'

METHODS = [  # (label, outcome folder, isolated-timing source, guarantee)
    ('capsule (1 per link)', CB / 'capsule', CB / 'capsule_timing', 'enclosing capsules, exact d_O'),
    ('spheres, RDF centers (54)', CB / 'spheres_enclosing', CB / 'spheres_enclosing_timing', 'enclosing spheres, exact d_O'),
    ('spheres, k-means 16/link', CB / 'spheres_kmeans16', CB / 'spheres_kmeans16_timing', 'enclosing spheres, exact d_O'),
    ('spheres, k-means 64/link', CB / 'spheres_kmeans64', CB / 'spheres_kmeans64_timing', 'enclosing spheres, exact d_O'),
    ('surface points + delta', CB / 'points_delta', CB / 'points_delta_timing', 'mesh within delta of points, exact d_O'),
    ('capsule, fitted field', CB / 'capsule_fitted', None, 'enclosing capsules, our fitted phi_O'),
    ('spheres RDF, fitted field', CB / 'spheres_enclosing_fitted', None, 'enclosing spheres, our fitted phi_O'),
    ('spheres k-means 16, fitted field', CB / 'spheres_kmeans16_fitted', None, 'enclosing spheres, our fitted phi_O'),
    ('spheres k-means 64, fitted field', CB / 'spheres_kmeans64_fitted', None, 'enclosing spheres, our fitted phi_O'),
    ('surface points + delta, fitted field', CB / 'points_delta_fitted', CB / 'points_delta_fitted_timing', 'mesh within delta of points, our fitted phi_O'),
    ('ours 12 mm, continuous-time rows', CU / 'free', None, 'surface cover, fitted fields'),
    ('ours 6 mm, continuous-time rows', CU / 'free_6mm_fixed', None, 'surface cover, fitted fields'),
    ('ours 12 mm, sampled-data', CU / 'sampled_zero_fixed', '12mm', 'surface cover + between samples'),
    ('ours 6 mm, sampled-data', CU / 'sampled_zero_6mm_fixed', '6mm', 'surface cover + between samples'),
]


def outcomes(d):
    recs = [json.loads(p.read_text()) for p in sorted(d.glob('franka_*.json'))]
    if not recs:
        return None
    m = [r['metrics'] for r in recs]
    gaps = np.array([x['min_gap_bound_mm'] for x in m])
    scale = 1e3 if d.parent == CU else 1.                # certificate_upgrades stores seconds, baselines ms
    t = np.concatenate([np.load(p)['t'] * scale for p in sorted(d.glob('franka_*.npz'))])
    return dict(trials=len(m), reached=sum(x['reached_s'] is not None for x in m),
                collision_trials=sum(x.get('nonpositive_gap_states', 0) > 0 for x in m),
                slack_or_zero_trials=sum((x.get('slack_steps', 0) + x.get('zero_steps', 0)) > 0 for x in m),
                min_gap_mm=float(gaps.min()), median_min_gap_mm=float(np.median(gaps)),
                rows_max=max(x['rows_max'] for x in m), t_median_ms=float(np.median(t)),
                t_p95_ms=float(np.percentile(t, 95)), timing='outcome runs (parallel)')


def isolated(src):
    if src is None:
        return None
    if isinstance(src, str):
        if not PT.exists():
            return None
        s = json.loads(PT.read_text()).get(src)
        return None if s is None else dict(t_median_ms=s['median_ms'], t_p95_ms=s['p95_ms'], t_max_ms=s['max_ms'])
    files = sorted(src.glob('franka_*.npz'))
    if len(files) < 30:
        return None
    t = np.concatenate([np.load(p)['t'] for p in files])
    return dict(t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)), t_max_ms=float(t.max()))


def main():
    rows = {}
    for label, d, tsrc, guarantee in METHODS:
        o = outcomes(d)
        if o is None:
            continue
        iso = isolated(tsrc)
        if iso:
            o.update(iso, timing='isolated (one process)')
        o['guarantee'] = guarantee
        setup = d / 'setup.json'
        if setup.exists():
            s = json.loads(setup.read_text())
            o['primitives'] = s['primitives']
            o['margin_mm_range'] = [round(min(s['margin_mm']), 1), round(max(s['margin_mm']), 1)]
        rows[label] = o
    (CB / 'summary.json').write_text(json.dumps(rows, indent=1))
    lines = ['| method | guarantee | margin [mm] | reached/30 | collision trials | min / median gap [mm] | '
             'time med / p95 [ms] | timing | max rows |', '|' + '---|' * 9]
    for k, o in rows.items():
        mg = '-'.join(str(x) for x in o.get('margin_mm_range', [])) or 'level 1-4'
        lines.append(f"| {k} | {o['guarantee']} | {mg} | {o['reached']}/{o['trials']} | {o['collision_trials']} | "
                     f"{o['min_gap_mm']:.1f} / {o['median_min_gap_mm']:.1f} | {o['t_median_ms']:.1f} / {o['t_p95_ms']:.1f} | "
                     f"{o['timing']} | {o['rows_max']} |")
    table = '\n'.join(lines)
    (CB / 'README.md').write_text(
        '# Continuous-boundary baselines vs. ours (Franka, 30 trials)\n\n'
        'See PLAN.md for the setup. Baselines use the exact obstacle SDF (1-Lipschitz); ours uses the fitted fields '
        'with certified levels. Baseline outcome runs audit every 10th state and the final state; ours every state. '
        'Times: filter only (rows + QP), Python, one thread.\n\n' + table + '\n')
    print(table)


if __name__ == '__main__':
    main()
