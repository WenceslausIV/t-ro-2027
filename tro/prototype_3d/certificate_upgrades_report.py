"""Aggregate results/certificate_upgrades/<variant>/franka_XX.json into summary.json and a Markdown table.

    python prototype_3d/certificate_upgrades_report.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'certificate_upgrades'
HIST = ROOT / 'results' / 'native_patch_evaluation'


def load(d):
    recs = [json.loads(p.read_text()) for p in sorted(d.glob('franka_*.json'))]
    recs = [r for r in recs if 'metrics' in r]
    for r in recs:                                   # pooled step timings from the saved arrays
        z = d / f"franka_{r['trial']:02d}.npz"
        if z.exists():
            r['t_steps_ms'] = (1e3 * np.load(z)['t']).tolist()
    return recs


def agg(recs):
    m = [r['metrics'] for r in recs]
    ts = [x['t_median_ms'] for x in m]
    pooled = np.concatenate([r['t_steps_ms'] for r in recs if 't_steps_ms' in r]) if any(
        't_steps_ms' in r for r in recs) else None
    return dict(pooled_median_ms=float(np.median(pooled)) if pooled is not None else None,
                pooled_p95_ms=float(np.percentile(pooled, 95)) if pooled is not None else None,
                pooled_over_10ms_percent=float(100 * np.mean(pooled > 10)) if pooled is not None else None,
                trials=len(m), reached=sum(x['reached_s'] is not None for x in m),
                reached_trials=[r['trial'] for r in recs if r['metrics']['reached_s'] is not None],
                slack_trials=sum(x.get('slack_steps', 0) > 0 for x in m),
                slack_steps=sum(x.get('slack_steps', 0) for x in m),
                zero_input_trials=sum(x.get('zero_steps', 0) > 0 for x in m),
                zero_input_steps=sum(x.get('zero_steps', 0) for x in m),
                min_gap_bound_mm=min(x['min_gap_bound_mm'] for x in m),
                median_trial_min_gap_mm=float(np.median([x['min_gap_bound_mm'] for x in m])),
                nonpositive_gap_states=sum(x['nonpositive_gap_states'] for x in m),
                median_of_trial_median_ms=float(np.median(ts)),
                max_trial_p95_ms=max(x['t_p95_ms'] for x in m),
                rows_max=max(x['rows_max'] for x in m), boxes_max=max(x['boxes_max'] for x in m),
                mean_arrival_s=float(np.mean([x['reached_s'] for x in m if x['reached_s'] is not None]))
                if any(x['reached_s'] is not None for x in m) else None)


def historical():
    recs = []
    for p in sorted(HIST.glob('franka_*.json')):
        r = json.loads(p.read_text())
        mt = r['metrics']
        recs.append(dict(trial=int(p.stem.split('_')[1]), metrics=dict(
            reached_s=mt['reached_s'], slack_steps=mt['slack_steps'], zero_steps=0,
            min_gap_bound_mm=min(mt['min_gap_bound_mm'], mt.get('final_gap_bound_mm') or np.inf),
            nonpositive_gap_states=mt.get('uncertified_gap_states', 0), t_median_ms=mt['t_median_ms'],
            t_p95_ms=mt['t_p95_ms'], rows_max=mt['rows_max'], boxes_max=mt['boxes_max'])))
    return recs


def main():
    rows = {'unit (historical, 12 mm)': agg(historical())}
    hs = json.loads((HIST / 'summary.json').read_text())      # pooled timings of the historical record
    rows['unit (historical, 12 mm)'].update(pooled_median_ms=hs['filter_median_ms'], pooled_p95_ms=hs['filter_p95_ms'])
    for d in sorted(p for p in OUT.iterdir() if p.is_dir()):
        recs = load(d)
        if recs:
            rows[d.name] = agg(recs)
    (OUT / 'summary.json').write_text(json.dumps(rows, indent=1))
    keys = ('trials', 'reached', 'slack_trials', 'slack_steps', 'zero_input_trials', 'zero_input_steps',
            'min_gap_bound_mm', 'median_trial_min_gap_mm', 'pooled_median_ms', 'pooled_p95_ms',
            'rows_max', 'mean_arrival_s')
    lines = ['| variant | ' + ' | '.join(keys) + ' |', '|---' * (len(keys) + 1) + '|']
    for name, a in rows.items():
        f = lambda v: f'{v:.2f}' if isinstance(v, float) else str(v)
        lines.append(f'| {name} | ' + ' | '.join(f(a[k]) for k in keys) + ' |')
    (OUT / 'summary.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
