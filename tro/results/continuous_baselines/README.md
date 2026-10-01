# Continuous-boundary baselines vs. ours (Franka, 30 trials)

See PLAN.md for the setup. Baselines use the exact obstacle SDF (1-Lipschitz); ours uses the fitted fields with certified levels. Baseline outcome runs audit every 10th state and the final state; ours every state. Times: filter only (rows + QP), Python, one thread.

| method | guarantee | margin [mm] | reached/30 | collision trials | min / median gap [mm] | time med / p95 [ms] | timing | max rows |
|---|---|---|---|---|---|---|---|---|
| capsule (1 per link) | enclosing capsules, exact d_O | 60.4-113.9 | 16/30 | 0 | 22.2 / 40.5 | 3.5 / 5.7 | outcome runs (parallel) | 58 |
| spheres, RDF centers (54) | enclosing spheres, exact d_O | 44.4-92.1 | 21/30 | 0 | 5.2 / 12.9 | 3.2 / 5.2 | outcome runs (parallel) | 22 |
| spheres, k-means 16/link | enclosing spheres, exact d_O | 40.1-67.9 | 22/30 | 0 | 15.6 / 26.5 | 3.2 / 5.0 | outcome runs (parallel) | 21 |
| spheres, k-means 64/link | enclosing spheres, exact d_O | 25.3-44.1 | 22/30 | 0 | 12.1 / 18.2 | 3.2 / 5.4 | outcome runs (parallel) | 54 |
| ours 12 mm, continuous-time rows | surface cover, fitted fields | level 1-4 | 8/30 | 0 | 12.5 / 21.0 | 11.5 / 42.0 | outcome runs (parallel) | 9045 |
| ours 12 mm, sampled-data | surface cover + between samples | level 1-4 | 7/30 | 0 | 12.2 / 20.9 | 14.6 / 33.7 | outcome runs (parallel) | 7611 |
| ours 6 mm, sampled-data | surface cover + between samples | level 1-4 | 10/30 | 0 | 3.8 / 10.4 | 22.0 / 72.5 | outcome runs (parallel) | 37986 |
