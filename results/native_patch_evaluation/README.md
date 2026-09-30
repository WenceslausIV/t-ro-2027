# Native SDF patches as the default cover

The fitted fields, certified levels, nominal controllers and fixed trials are unchanged.
Cover side equals the source field cell size; the vertex certificate has unit lift and
no online refinement. Previous subdivided results and caches are preserved.

| Metric | Native patches | Historical subdivided implementation |
|---|---:|---:|
| Dock cover side / boxes | 10 mm / 772 | 5 mm / 1548 |
| Dock nose depth | 0.691 m | 0.698 m |
| Dock minimum polygon gap | 9.1 mm | 1.7 mm |
| Dock maximum rows | 768 | 1604 |
| Dock slack steps | 0 | 0 |
| Franka initial cover side / boxes | 12 mm / 6847 | 6 mm / 26636 |
| Franka goals reached within 10 s | 6/30 | 17/30 |
| Franka trials requiring slack | 2/30 | 0/30 |
| Franka minimum saved-state mesh distance bound | 15.3 mm | 1.1 mm |
| Franka maximum rows | 2064 | 34992 |

Native docking seats the nose (required depth 0.668 m); its nominal target deliberately lies
beyond feasible seating. Franka has 0 trials with a negative
saved-state distance bound, but 2 trials use slack (1497 steps).
Thus geometric audits and satisfaction of the hard CBF constraints must be distinguished.
The native setting is simpler and has fewer rows, with substantially more conservative motion.
It does not establish equivalent performance to subdivision.

Native filter median/p95: docking 3.08/3.67 ms;
Franka 10.54/16.68 ms, pooled across executed steps.
Each process uses one CPU thread. Franka trials were evaluated in concurrent batches with
other background activity; these timings are not isolated or directly comparable to the
historical paper timings. Collision audits and rendering are excluded from filter time.

Franka distance audits cover the subdivided link triangles at saved states, including the
final state, using the analytic obstacle distance and triangle radii. Bounds are capped at
45 mm to account for geometrically pruned centroids (50-mm AABB padding, maximum triangle
radius 2.666 mm). They do not certify between-step motion.
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
