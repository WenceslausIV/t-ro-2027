# Continuous-boundary baselines vs. ours (Franka, 30 trials)

See PLAN.md for the setup. Baselines use the exact obstacle SDF (1-Lipschitz); ours uses the fitted fields with certified levels. Baseline outcome runs audit every 10th state and the final state; ours every state. Times: filter only (rows + QP), Python, one thread.

| method | guarantee | margin [mm] | reached/30 | collision trials | min / median gap [mm] | time med / p95 [ms] | timing | max rows |
|---|---|---|---|---|---|---|---|---|
| capsule (1 per link) | enclosing capsules, exact d_O | 60.4-113.9 | 16/30 | 0 | 22.2 / 40.5 | 3.3 / 5.0 | isolated (one process) | 58 |
| spheres, RDF centers (54) | enclosing spheres, exact d_O | 44.4-92.1 | 21/30 | 0 | 5.2 / 12.9 | 3.0 / 4.5 | isolated (one process) | 22 |
| spheres, k-means 16/link | enclosing spheres, exact d_O | 40.1-67.9 | 22/30 | 0 | 15.6 / 26.5 | 3.1 / 4.6 | isolated (one process) | 21 |
| spheres, k-means 64/link | enclosing spheres, exact d_O | 25.3-44.1 | 22/30 | 0 | 12.1 / 18.2 | 3.0 / 4.8 | isolated (one process) | 54 |
| surface points + delta | mesh within delta of points, exact d_O | 5.8-5.8 | 22/30 | 0 | 3.2 / 4.4 | 22.4 / 77.2 | isolated (one process) | 6874 |
| capsule, fitted field | enclosing capsules, our fitted phi_O | 60.4-113.9 | 10/30 | 0 | 32.4 / 47.3 | 3.0 / 6.4 | outcome runs (parallel) | 93 |
| spheres RDF, fitted field | enclosing spheres, our fitted phi_O | 44.4-92.1 | 19/30 | 0 | 32.4 / 47.3 | 3.0 / 4.8 | outcome runs (parallel) | 25 |
| spheres k-means 16, fitted field | enclosing spheres, our fitted phi_O | 40.1-67.9 | 17/30 | 0 | 32.4 / 47.3 | 3.2 / 4.9 | outcome runs (parallel) | 18 |
| spheres k-means 64, fitted field | enclosing spheres, our fitted phi_O | 25.3-44.1 | 21/30 | 0 | 32.4 / 44.4 | 3.4 / 5.5 | outcome runs (parallel) | 63 |
| surface points + delta, fitted field | mesh within delta of points, our fitted phi_O | 5.8-5.8 | 22/30 | 0 | 9.2 / 12.8 | 48.5 / 76.9 | isolated (one process) | 7398 |
| ours 12 mm, continuous-time rows | surface cover, fitted fields | level 1-4 | 8/30 | 0 | 12.5 / 21.0 | 11.5 / 42.0 | outcome runs (parallel) | 9045 |
| ours 12 mm, sampled-data | surface cover + between samples | level 1-4 | 7/30 | 0 | 12.2 / 20.9 | 14.6 / 33.7 | outcome runs (parallel) | 7611 |
| ours 6 mm, sampled-data | surface cover + between samples | level 1-4 | 10/30 | 0 | 3.8 / 10.4 | 22.0 / 72.5 | outcome runs (parallel) | 37986 |
