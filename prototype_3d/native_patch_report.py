"""Summarize the saved native-patch evaluations without rerunning a controller."""
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/native_patch_evaluation'


def main():
    paths = [OUT / f'franka_{i:02d}.json' for i in range(30)]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        raise SystemExit('Awaiting: ' + ', '.join(missing))
    records = [json.loads(p.read_text()) for p in paths]
    metrics = [r['metrics'] for r in records]
    assert len({r['metadata']['geometry_sha256'] for r in records}) == 1
    assert all(r['metadata']['refinement_depth'] == 0 for r in records)
    assert all(abs(l['side_mm'] - 12.) < 1e-9
               for r in records for l in r['metadata']['setup']['links'])

    # real_gap prunes centroids outside the obstacle AABB expanded by 50 mm.
    # Those triangles have distance >= 50 mm - their radius. Keep a 45-mm
    # cap in the reported bounds, including states for which all pairs prune.
    cache = np.load(ROOT / 'prototype_3d/cache_franka.npz', allow_pickle=True)['d'].item()
    max_radius = max(L['gt_r'].max() for L in cache['links'])
    assert max_radius < .005
    cap_mm = 45.
    minima = [min(cap_mm, m['min_gap_bound_mm'] if m['min_gap_bound_mm'] is not None else cap_mm,
                  m['final_gap_bound_mm'] if m['final_gap_bound_mm'] is not None else cap_mm)
              for m in metrics]
    times = []
    for p, m in zip(paths, metrics):
        z = np.load(p.with_suffix('.npz'))
        assert len(z['q']) == len(z['u']) + 1
        assert len(z['u']) == m['steps'] == len(z['time_s'])
        assert np.isfinite(z['q']).all() and np.isfinite(z['u']).all()
        times.extend(1000 * z['time_s'])
    summary = dict(trials=30, reached=sum(m['reached_s'] is not None for m in metrics),
                   trials_needing_slack=sum(m['slack_steps'] > 0 for m in metrics),
                   slack_steps=sum(m['slack_steps'] for m in metrics),
                   trials_with_negative_distance_bound=sum(g < 0 for g in minima),
                   min_gap_bound_mm=min(minima), median_min_gap_bound_mm=float(np.median(minima)),
                   filter_median_ms=float(np.median(times)), filter_p95_ms=float(np.percentile(times,95)),
                   max_rows=max(m['rows_max'] for m in metrics),
                   max_active_boxes=max(m['boxes_max'] for m in metrics),
                   cover_boxes=sum(l['n_boxes'] for l in records[0]['metadata']['setup']['links']),
                   distance_cap_mm=cap_mm, maximum_triangle_radius_mm=1000*max_radius,
                   note='Saved-state mesh bounds, including final states; no between-step guarantee. '
                        'One thread per process; trials ran in concurrent batches with other background '
                        'activity. Timings are exploratory, not isolated speed comparisons.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    with (OUT/'trials.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['trial','reached_s','mesh_gap_lower_bound_mm','slack_steps','max_rows',
                         'filter_median_ms','filter_p95_ms'])
        for i,(m,g) in enumerate(zip(metrics,minima)):
            writer.writerow([i,m['reached_s'],g,m['slack_steps'],m['rows_max'],
                             m['t_median_ms'],m['t_p95_ms']])
    dock=json.loads((OUT/'dock.json').read_text())['metrics']
    s=summary
    text=f'''# Native SDF patches as the default cover

The fitted fields, certified levels, nominal controllers and fixed trials are unchanged.
Cover side equals the source field cell size; the vertex certificate has unit lift and
no online refinement. Previous subdivided results and caches are preserved.

| Metric | Native patches | Historical subdivided implementation |
|---|---:|---:|
| Dock cover side / boxes | 10 mm / 772 | 5 mm / 1548 |
| Dock nose depth | {dock['tip_depth_m']:.3f} m | 0.698 m |
| Dock minimum polygon gap | {dock['min_gap_bound_mm']:.1f} mm | 1.7 mm |
| Dock maximum rows | {dock['rows_max']} | 1604 |
| Dock slack steps | {dock['slack_steps']} | 0 |
| Franka initial cover side / boxes | 12 mm / {s['cover_boxes']} | 6 mm / 26636 |
| Franka goals reached within 10 s | {s['reached']}/30 | 17/30 |
| Franka trials requiring slack | {s['trials_needing_slack']}/30 | 0/30 |
| Franka minimum saved-state mesh distance bound | {s['min_gap_bound_mm']:.1f} mm | 1.1 mm |
| Franka maximum rows | {s['max_rows']} | 34992 |

Native docking seats the nose (required depth 0.668 m); its nominal target deliberately lies
beyond feasible seating. Franka has {s['trials_with_negative_distance_bound']} trials with a negative
saved-state distance bound, but {s['trials_needing_slack']} trials use slack ({s['slack_steps']} steps).
Thus geometric audits and satisfaction of the hard CBF constraints must be distinguished.
The native setting is simpler and has fewer rows, with substantially more conservative motion.
It does not establish equivalent performance to subdivision.

Native filter median/p95: docking {dock['t_median_ms']:.2f}/{dock['t_p95_ms']:.2f} ms;
Franka {s['filter_median_ms']:.2f}/{s['filter_p95_ms']:.2f} ms, pooled across executed steps.
Each process uses one CPU thread. Franka trials were evaluated in concurrent batches with
other background activity; these timings are not isolated or directly comparable to the
historical paper timings. Collision audits and rendering are excluded from filter time.

Franka distance audits cover the subdivided link triangles at saved states, including the
final state, using the analytic obstacle distance and triangle radii. Bounds are capped at
45 mm to account for geometrically pruned centroids (50-mm AABB padding, maximum triangle
radius {s['maximum_triangle_radius_mm']:.3f} mm). They do not certify between-step motion.
All calculations use double precision without outward rounding.

Reproduce: `python prototype_3d/native_patch_evaluation.py`

Summarize and validate saved arrays: `python prototype_3d/native_patch_report.py`

The separate two-robot animation and performance comparison are in
`../two_robot_cover_comparison/`: identical 40-mm SDFs with native 40-mm covers versus
5-mm cover subdivision, no online refinement in either run.

For an explicit 6-mm Franka cover, use `F.build(cover_side=.006)`; optional online refinement
is enabled with `SUMMED_REFINE_DEPTH=2`. Historical certificate ablations explicitly retain
their cached covers and depth=2. Legacy Franka/dual figure and GIF scripts also select the
historical 6-mm covers and two refinement levels.
'''
    (OUT/'README.md').write_text(text,encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
