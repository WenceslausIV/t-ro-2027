"""Fair comparison table for results/continuous_baselines/PLAN.md (FAIRNESS REQUIREMENTS, LOCAL CONTINUATION).

    python prototype_3d/continuous_baselines_report.py
Main table: every method on OUR fitted obstacle fields, rows enforced at the 10-ms samples, mesh audit at every
state, trials that start outside the method's certified set excluded (start_validity.json), filter times from the
isolated runs on this machine. Extra row: ours with the sampled-data certificate. Reference rows: baselines on the
exact obstacle SDF (different obstacle information; not part of the fair comparison).
Writes results/continuous_baselines/summary.json and README.md.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CB = ROOT / 'results' / 'continuous_baselines'
CU = ROOT / 'results' / 'certificate_upgrades'
PT = ROOT / 'results' / 'patch_size_timing' / 'summary.json'
VALID = CB / 'start_validity.json'
BUDGET_MS = 50.   # VLA-style 20-Hz loop: a setting qualifies if its isolated p95 filter time is within the budget

# (label, outcome folder, isolated timing: folder or patch_size_timing key, validity key, guarantee)
FAIR = [
    ('ours, 12-mm patches', CU / 'free_fixed', '12mm_free', 'ours_12mm', 'surface cover, rows at samples'),
    ('ours, 6-mm patches', CU / 'free_6mm_fixed', '6mm_free', 'ours_6mm', 'surface cover, rows at samples'),
    ('ours, 8-mm patches', CU / 'free_8mm_fixed', '8mm_free', 'ours_8mm', 'surface cover, rows at samples'),
    ('ours, 12 mm + one near-contact halving', CU / 'free_fixed_refine1', '12mm_free_refine1', 'ours_12mm',
     'surface cover, rows at samples'),
    ('ours, 6 mm + one near-contact halving', CU / 'free_6mm_fixed_refine1', '6mm_free_refine1', 'ours_6mm',
     'surface cover, rows at samples'),
    ('surface points + delta, 10-mm edges', CB / 'points_delta_fitted_audit1', CB / 'points_delta_fitted_timing_local',
     'points_delta_fitted', 'mesh within delta of points'),
    ('surface points + delta, 5-mm edges', CB / 'points_delta_fitted_e5', CB / 'points_delta_fitted_e5_timing_local',
     'points_delta_fitted_e5', 'mesh within delta of points'),
    ('surface points + delta, 15-mm edges', CB / 'points_delta_fitted_e15', CB / 'points_delta_fitted_e15_timing_local',
     'points_delta_fitted_e15', 'mesh within delta of points'),
    ('surface points + delta, 20-mm edges', CB / 'points_delta_fitted_e20', CB / 'points_delta_fitted_e20_timing_local',
     'points_delta_fitted_e20', 'mesh within delta of points'),
    ('capsule, 1 per link', CB / 'capsule_fitted_audit1', CB / 'capsule_fitted_timing_local', 'capsule_fitted',
     'enclosing capsules'),
    ('spheres, RDF centers (54)', CB / 'spheres_enclosing_fitted_audit1', CB / 'spheres_enclosing_fitted_timing_local',
     'spheres_enclosing_fitted', 'enclosing spheres'),
    ('spheres, k-means 16/link', CB / 'spheres_kmeans16_fitted_audit1', CB / 'spheres_kmeans16_fitted_timing_local',
     'spheres_kmeans16_fitted', 'enclosing spheres'),
    ('spheres, k-means 64/link', CB / 'spheres_kmeans64_fitted_audit1', CB / 'spheres_kmeans64_fitted_timing_local',
     'spheres_kmeans64_fitted', 'enclosing spheres'),
    ('spheres, k-means 128/link', CB / 'spheres_kmeans128_fitted_audit1', CB / 'spheres_kmeans128_fitted_timing_local',
     'spheres_kmeans128_fitted', 'enclosing spheres'),
]
EXTRA = [
    ('ours, 12 mm, plus between-sample guarantee', CU / 'sampled_zero_fixed', '12mm', 'ours_12mm',
     'surface cover + sampled-data'),
    ('ours, 6 mm, plus between-sample guarantee', CU / 'sampled_zero_6mm_fixed', '6mm', 'ours_6mm',
     'surface cover + sampled-data'),
]
REFERENCE = [  # exact obstacle SDF: more obstacle information than ours; audit every 10th state (cloud runs)
    ('[exact SDF] surface points + delta', CB / 'points_delta', CB / 'points_delta_timing_local', 'points_delta',
     'mesh within delta of points'),
    ('[exact SDF] capsule, 1 per link', CB / 'capsule', CB / 'capsule_timing_local', 'capsule', 'enclosing capsules'),
    ('[exact SDF] spheres, RDF centers', CB / 'spheres_enclosing', CB / 'spheres_enclosing_timing_local',
     'spheres_enclosing', 'enclosing spheres'),
    ('[exact SDF] spheres, k-means 16', CB / 'spheres_kmeans16', CB / 'spheres_kmeans16_timing_local',
     'spheres_kmeans16', 'enclosing spheres'),
    ('[exact SDF] spheres, k-means 64', CB / 'spheres_kmeans64', CB / 'spheres_kmeans64_timing_local',
     'spheres_kmeans64', 'enclosing spheres'),
]


def outcomes(d, invalid):
    recs = [json.loads(p.read_text()) for p in sorted(d.glob('franka_*.json'))]
    if not recs:
        return None
    done = len(recs)
    recs = [r for r in recs if r['trial'] not in invalid]
    m = [r['metrics'] for r in recs]
    gaps = np.array([x['min_gap_bound_mm'] for x in m])
    return dict(trials_done=done, valid_trials=len(m),
                reached=sum(x['reached_s'] is not None for x in m),
                collision_trials=sum(x.get('nonpositive_gap_states', 0) > 0 for x in m),
                slack_or_zero_trials=sum((x.get('slack_steps', 0) + x.get('zero_steps', 0)) > 0 for x in m),
                min_gap_mm=float(gaps.min()), median_min_gap_mm=float(np.median(gaps)),
                rows_max=max(x['rows_max'] for x in m))


def isolated(src):
    if isinstance(src, str):
        s = json.loads(PT.read_text()).get(src) if PT.exists() else None
        return None if s is None else dict(t_median_ms=s['median_ms'], t_p95_ms=s['p95_ms'], t_max_ms=s['max_ms'])
    files = sorted(src.glob('franka_*.npz'))
    if len(files) < 30:
        return None
    t = np.concatenate([np.load(p)['t'] for p in files])
    return dict(t_median_ms=float(np.median(t)), t_p95_ms=float(np.percentile(t, 95)), t_max_ms=float(t.max()))


def table(rows):
    lines = ['| method | guarantee | primitives | valid starts | reached | collision trials | slack trials | '
             'min / median gap [mm] | filter time med / p95 [ms] | 20-Hz budget | max rows |', '|' + '---|' * 11]
    for k, o in rows.items():
        t = (f"{o['t_median_ms']:.1f} / {o['t_p95_ms']:.1f}" if 't_median_ms' in o else 'pending')
        ok = '-' if 't_p95_ms' not in o else ('yes' if o['t_p95_ms'] <= BUDGET_MS else 'no')
        done = '' if o['trials_done'] == 30 else f" ({o['trials_done']} of 30 done)"
        lines.append(f"| {k} | {o['guarantee']} | {o.get('primitives', '-')} | {o['valid_trials']}{done} | "
                     f"{o['reached']}/{o['valid_trials']} | {o['collision_trials']} | {o['slack_or_zero_trials']} | "
                     f"{o['min_gap_mm']:.1f} / {o['median_min_gap_mm']:.1f} | {t} | {ok} | {o['rows_max']} |")
    return '\n'.join(lines)


def collect(spec, validity):
    rows = {}
    for label, d, tsrc, vkey, guarantee in spec:
        invalid = set(validity.get(vkey, {}).get('outside_trials', []))
        o = outcomes(d, invalid)
        if o is None:
            continue
        iso = isolated(tsrc)
        if iso:
            o.update(iso)
        o['guarantee'] = guarantee
        o['outcome_folder'] = str(d.relative_to(ROOT))
        setup = d / 'setup.json'
        if setup.exists():
            o['primitives'] = json.loads(setup.read_text())['primitives']
        rows[label] = o
    return rows


def main():
    validity = json.loads(VALID.read_text()) if VALID.exists() else {}
    fair, extra, ref = collect(FAIR, validity), collect(EXTRA, validity), collect(REFERENCE, validity)
    within = {k: o for k, o in fair.items() if o.get('t_p95_ms', 1e9) <= BUDGET_MS}
    over = {k: o for k, o in fair.items() if k not in within}
    (CB / 'summary.json').write_text(json.dumps(dict(fair=fair, extra=extra, reference=ref), indent=1))
    text = (
        '# Continuous-boundary baselines vs. ours (Franka, 30 fixed trials)\n\n'
        'Setup and fairness requirements: PLAN.md. Generated by `python prototype_3d/continuous_baselines_report.py`.\n\n'
        '## Fair comparison\n\n'
        'All methods use OUR fitted obstacle fields (level and gradient bounds), enforce their rows at the 10-ms samples, '
        'use the same trials, nominal controller, dt, gamma, input bounds, horizon, and QP (DAQP, exact check, slack '
        'fallback). Every state is audited (certified mesh-distance lower bound). Trials that start outside a '
        "method's certified set (h < 0 at q0; start_validity.json) are excluded and counted under valid starts. "
        'Filter times: rows + QP, one process, one thread, nothing else running, AMD Ryzen 9 3900X. '
        'Reached/collision/slack counts are over the valid trials. Real-time criterion (VLA-style 20-Hz loop): '
        f'a setting qualifies if its isolated p95 filter time is at most {BUDGET_MS:.0f} ms.\n\n'
        '### Settings within the 20-Hz budget (the comparison)\n\n' + table(within) + '\n\n'
        '### Settings over the budget or not yet timed\n\n' + table(over) + '\n\n'
        '## Extra row: ours with the sampled-data certificate (also certifies the motion between samples)\n\n'
        + table(extra) + '\n\n'
        '## Reference only: baselines on the exact obstacle SDF\n\n'
        'These use more obstacle information than ours (the analytic 1-Lipschitz distance of the boxes, no fitting '
        'level), and were audited at every 10th state (cloud runs). They are not part of the fair comparison.\n\n'
        + table(ref) + '\n')
    (CB / 'README.md').write_text(text)
    print(text)


if __name__ == '__main__':
    main()
