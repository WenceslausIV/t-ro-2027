"""Table for results/star_random_goals (python prototype_3d/star_random_goals_report.py).

Counts are out of all 30 goals. A trial whose start is outside the method's certified set (h0 < 0) counts as a
failure and is listed under invalid starts. Filter times come from the simulation runs, which are executed one
process at a time with nothing else running (no audit inside them). Writes README.md and summary.json.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / 'results' / 'star_random_goals'
BUDGET_MS = 50.
REUSE_MM = 20.   # audit evaluates directly below this distance; larger bounds are certified but not tight


def fmt(x):
    return 'n/a' if not np.isfinite(x) else (f'>= {REUSE_MM:.0f}' if x >= REUSE_MM - 1e-9 else f'{x:.1f}')
METHODS = [('ours: vertex certificate, box = 6-mm patch', 'ours'),
           ('spheres: 64/link, tube 256', 'spheres_L64_T256'),
           ('spheres: 64/link, tube 1024', 'spheres_L64_T1024'),
           ('spheres: 128/link, tube 1024', 'spheres_L128_T1024')]


def main():
    rows, lines = {}, ['| method | invalid starts | reached (valid) | collision trials | slack trials | '
                       'min / median gap [mm] | filter time med / p95 [ms] | 20-Hz budget | max rows |',
                       '|' + '---|' * 9]
    for label, f in METHODS:
        recs = [json.loads(p.read_text())['metrics'] for p in sorted((D / f).glob('trial_*.json'))]
        if not recs:
            continue
        valid = [m for m in recs if m['h0_mm'] >= 0]
        t = np.concatenate([np.load(p)['t'] for p in sorted((D / f).glob('trial_*.npz'))])
        audited = all('min_gap_bound_mm' in m for m in recs)
        gaps = np.array([m['min_gap_bound_mm'] for m in valid]) if audited and valid else np.array([np.nan])
        o = dict(trials=len(recs), invalid_starts=len(recs) - len(valid),
                 reached=sum(m['reached_s'] is not None for m in valid),
                 collision_trials=sum(m.get('collision_witness_states', 0) > 0 for m in valid) if audited else None,
                 slack_trials=sum(m['slack_steps'] > 0 for m in valid),
                 min_gap_mm=float(np.nanmin(gaps)), median_gap_mm=float(np.nanmedian(gaps)),
                 t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)),
                 rows_max=max(m['rows_max'] for m in recs))
        rows[label] = o
        lines.append(f"| {label} | {o['invalid_starts']} | {o['reached']}/30 ({o['reached']}/{len(valid)}) | "
                     f"{o['collision_trials'] if audited else 'pending'} | {o['slack_trials']} | "
                     f"{fmt(o['min_gap_mm'])} / {fmt(o['median_gap_mm'])} | {o['t_median_ms']:.1f} / {o['t_p95_ms']:.1f} | "
                     f"{'yes' if o['t_p95_ms'] <= BUDGET_MS else 'no'} | {o['rows_max']} |")
    (D / 'summary.json').write_text(json.dumps(rows, indent=1))
    text = ('# Star tube, 30 random goals: ours vs full sphere decomposition\n\n'
            'The arm starts threaded in the rounded star tube and leaves to 30 random goals whose unfiltered motion '
            'hits the tube (goals.json). Spheres: links AND tube decomposed into certified enclosing spheres. Same '
            'nominal, dt, gamma, input bounds, horizon (15 s), activation, QP (DAQP, slack fallback). Reached counts '
            'are out of 30; a start outside the certified set (h0 < 0) is a failure. Gaps: certified mesh-distance '
            'lower bounds at every state over valid trials; the audit is exact below 20 mm, so larger values are shown '
            'as >= 20. Times: isolated simulation runs, filter only.\n\n'
            + '\n'.join(lines) + '\n')
    (D / 'README.md').write_text(text)
    print(text)


if __name__ == '__main__':
    main()
