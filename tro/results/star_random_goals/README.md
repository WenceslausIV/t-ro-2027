# Star tube, 30 random goals: ours vs full sphere decomposition

The arm starts threaded in the rounded star tube and leaves to 30 random goals whose unfiltered motion hits the tube (goals.json). Spheres: links AND tube decomposed into certified enclosing spheres. Same nominal, dt, gamma, input bounds, horizon (15 s), activation, QP (DAQP, slack fallback). Reached counts are out of 30; a start outside the certified set (h0 < 0) is a failure. Gaps: certified mesh-distance lower bounds at every state over valid trials; the audit is exact below 20 mm, so larger values are shown as >= 20. Times: isolated simulation runs, filter only.

| method | invalid starts | reached (valid) | collision trials | slack trials | min / median gap [mm] | filter time med / p95 [ms] | 20-Hz budget | max rows |
|---|---|---|---|---|---|---|---|---|
| ours: vertex certificate, box = 6-mm patch | 0 | 14/30 (14/30) | 0 | 0 | 6.9 / 15.6 | 22.9 / 36.2 | yes | 9222 |
| spheres: 64/link, tube 256 | 30 | 0/30 (0/0) | 0 | 0 | n/a / n/a | 5.1 / 6.2 | yes | 312 |
| spheres: 64/link, tube 1024 | 0 | 8/30 (8/30) | 0 | 0 | >= 20 / >= 20 | 18.9 / 23.2 | yes | 793 |
| spheres: 128/link, tube 1024 | 0 | 10/30 (10/30) | 0 | 0 | >= 20 / >= 20 | 41.6 / 52.8 | no | 1025 |
