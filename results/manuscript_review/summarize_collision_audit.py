"""Combine witness replays and adaptive separation audits; preserve original results."""
import hashlib
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def main():
    witness_path = OUT / 'check_collision_witnesses.json'
    refined_path = OUT / 'refined_collision_audit.json'
    witness = json.loads(witness_path.read_text())
    refined = json.loads(refined_path.read_text())
    baseline_path = ROOT / 'prototype_3d/franka_baselines.json'
    ours_path = ROOT / 'prototype_3d/franka_results.json'
    baseline = json.loads(baseline_path.read_text())
    ours = json.loads(ours_path.read_text())
    summary = dict(
        scope='Historical evaluated saved states only; no between-step or final unlogged-state claim.',
        collision_definition='Negative attained primitive-SDF minimum at a mesh point; large negative lower bounds imply such a centroid witness using the maximum triangle radius.',
        separation_definition='Positive whole-triangle bounds, with adaptive subdivision when the coarse bound is inconclusive.',
        arithmetic='Double precision, without outward rounding.',
        nominal_collision_trials=witness['nominal_collisions_implied_by_centroid_witness'],
        methods={},
        sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in (witness_path, refined_path, baseline_path, ours_path)})
    for method, trials in witness['methods'].items():
        groups = dict(collision=[], separated=[], unresolved=[])
        for trial in trials:
            status = trial['status']
            if status.startswith('collision'):
                outcome = 'collision'
            elif status == 'positive_historical_saved_state_bound':
                assert trial['historical_lower_bound_mm'] > 0
                outcome = 'separated'
            else:
                audit = refined[f'{method}_trial_{trial["trial"]:02d}_witness']
                assert abs(trial['bound_difference_mm']) < 1e-9
                outcome = {'separated_at_saved_states': 'separated',
                           'collision': 'collision'}.get(audit['status'], 'unresolved')
            groups[outcome].append(trial['trial'])
        summary['methods'][method] = dict(counts={k: len(v) for k, v in groups.items()}, trials=groups)
        assert len(trials) == 30 and not groups['unresolved']
    for method in ('spheres', 'points'):
        trials = baseline[method]['trials']
        assert len(trials) == 30 and all(t['min_gap_mm'] > 0 for t in trials)
        summary['methods'][method + '_margin_1cm'] = dict(
            counts=dict(collision=0, separated=30, unresolved=0),
            minimum_saved_state_lower_bound_mm=min(t['min_gap_mm'] for t in trials))
    assert len(ours['trials']) == 30 and all(t['min_gap_mm'] > 0 for t in ours['trials'])
    summary['methods']['ours_subdivided'] = dict(
        counts=dict(collision=0, separated=30, unresolved=0),
        minimum_saved_state_lower_bound_mm=min(t['min_gap_mm'] for t in ours['trials']))
    (OUT / 'franka_collision_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v['counts'] for k, v in summary['methods'].items()}, indent=2))


if __name__ == '__main__':
    main()
