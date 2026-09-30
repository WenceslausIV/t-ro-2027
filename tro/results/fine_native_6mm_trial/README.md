# Refit 6-mm native SDF patches: one Franka trial

Trial 2 was selected before execution: the first trial that succeeded with the historical
subdivided controller and stalled with 12-mm native patches. All eight link SDFs were **refitted**
with 6-mm cells and 3-mm training spacing using the existing distance-target routine; levels
were recertified on the meshes. Cover boxes are the native 6-mm cells. No knot insertion,
additional cover subdivision or online refinement was used. Obstacle SDFs remain at 20 mm.

| Setting | Goal arrival | Slack steps |
|---|---:|---:|
| Historical 12-mm SDF / 12-mm native cover | None | 0 |
| Historical 12-mm SDF / 6-mm cover + optional refinement | 5.21 s | 0 |
| Refitted 6-mm SDF / 6-mm native cover | 5.13 s | 0 |

Native cover: 26182 boxes; maximum 191 active boxes,
1528 rows. Minimum saved-state mesh-distance lower bound: 29.727 mm.
The audit includes the final state, covers mesh triangles, caps geometrically pruned distances
at 45 mm, and does not certify between-step motion. Double precision, no outward rounding.

| Controller timing | ms |
|---|---:|
| Median | 5.880 |
| p95 | 17.168 |
| p99 | 18.508 |
| Maximum | 20.443 |
| Mean | 7.397 |

10-ms budget exceeded on 138/513 steps (26.90%).
Active-step median/p95: 13.414/18.066 ms.
This is a measured deadline check in Python, not a hard-real-time guarantee.
One CPU thread, one evaluation process, other background activity not excluded.
Filter timing includes nominal command, FK, pruning, spline evaluations, constraint assembly
and QP; collision audit and rendering run separately. Offline link fitting/certification:
497.78 s; additional barrier preparation: 15.47 s.

Reproduce: `python prototype_3d/fine_native_trial.py` (reuses this directory's fine-field cache).
The historical cache and manuscript are unchanged. One selected trial does not estimate a success rate.
