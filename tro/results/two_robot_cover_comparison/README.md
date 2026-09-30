# Two-robot fixed-cover comparison

Same ground-truth shapes (seed 3, shape indices 0 and 3), 40 mm spline patches,
5 mm training spacing, certified levels, starts/goals, input bounds, gamma=5,
50 mm activation threshold, unit lifted vertex certificate, dt=10 ms, 12 s horizon.
Only the cover size changes: 40 mm native patches vs 5 mm subdivisions.
Automatic refinement is OFF. One direction per pair: robot B's surface against A's field.
The GIF shows all B cover boxes; only active boxes contribute QP rows.

| Cover | Both goals [s] | Min saved gap [mm] | Interpolated lower bound [mm] | Max rows | Filter median/p95 [ms] | Active filter median/p95 [ms] | Input modification J | Slack steps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Patch as cover: 40 mm | None | 21.403 | 14.984 | 32 | 1.808/2.558 | 1.854/2.566 | 15.746 | 0 |
| Subdivided cover: 5 mm | None | 4.516 | 0.072 | 628 | 3.306/4.679 | 3.389/4.696 | 15.195 | 0 |

Unfiltered baseline: 117 of 1201 saved states intersect.
Goal tolerance: position < 5 cm AND orientation error < 0.05 rad, both robots.
J = integral ||u-u_nom||_2 dt; mixes translational and angular coordinates identically in both runs.
Timings cover nominal command, constraint assembly, and QP, excluding collision audit and drawing.
Single CPU thread; unrelated background activity was not excluded. Active-step timings use each run's own active steps.
This is one predetermined crossing scene, not a success-rate benchmark or proof of universal improvement.
Polygon vertex-to-edge distances, containment, and edge-crossing tests check saved states. Interpolated lower bounds use
rigid-motion bounds and recursive bisection for the simulated linear pose interpolation in double precision;
they are not a sampled-data controller theorem or an outward-rounded numerical certificate.

Reproduce: `python prototype_3d/two_robot_cover_comparison.py`
Render saved results: `python prototype_3d/two_robot_cover_comparison.py --render-only`
