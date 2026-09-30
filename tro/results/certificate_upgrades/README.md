# Certificate upgrades on the 30 fixed Franka trials

Runner: `python prototype_3d/certificate_upgrades.py --variant {unit|free|sampled} [--fallback zero] [--field 6mm]`.
Aggregate: `python prototype_3d/certificate_upgrades_report.py` (writes `summary.json`, `summary.md`).

All variants keep the obstacle fields and levels, the nominal controller (gamma = 5, dt = 10 ms,
|qdot| <= 1 rad/s), the activation threshold eta = 3 cm, and the 30 trials of `prototype_3d/franka_trials.json`.
Covers are native SDF patches (no subdivision, no online refinement).

| Directory | Link fields | Certificate | Infeasible step |
|---|---|---|---|
| `unit/` | historical 12 mm | vertex rows, unit lift (paper's native default) | slack |
| `free/` | historical 12 mm | joint Bernstein rows, one multiplier w in [0, 3] per active box (QP variable) | slack |
| `sampled_zero/` | historical 12 mm | `free` + sampled-data tightening (Theorem sampled) | u = 0 |
| `sampled_zero_6mm/` | refitted 6 mm (`results/fine_native_6mm_trial`) | same as `sampled_zero` | u = 0 |

`unit/` reproduces `results/native_patch_evaluation/` bit for bit in trajectories and distance bounds
(checked on trials 2 and 5 before the runs); only timings differ (different machine). It was stopped after
five trials to free a core; the historical 30-trial record is used for the unit lift.

Sampled data (`prototype_3d/sampled_data.py`): inputs held for dt; per-link lever bound D (every link point and
joint origin of its chain), |yddot| <= 1.5 D |u|_1^2; caps on the relative twist and on |u|_1, enforced in the QP
and grown by 1.5x per step from the previous input (floors 0.05 m/s, 0.1 rad/s, 0.3 rad/s; certified point speed
<= 3 m/s); Hessian and gradient bounds of the obstacle field on the swept ball of radius r + travel; activation
threshold max(eta, G travel), G = 2.68 the largest cellwise gradient bound; domain-collar bounds 106.6--127 mm
exceed G travel_max = 80 mm. With the zero-input fallback, Corollary sampled covers every executed step.

Numerical validation of the implementation (not a proof): at the feasibility boundary of the rows (inputs scaled
until some multiplier barely satisfies them), F + gamma w g_A stayed positive at 7.2e5 random box points for both
the free and the sampled rows; the sampled-data curvature bound exceeded the observed |e''| by at least 8.5x and
the travel bound the observed travel by at least 2.5x along the saved native trajectories.

Saved-state mesh lower bounds use `franka3d.real_gap` at every executed state and the final state.
Timings: one thread per process, several processes concurrently on an Intel Xeon @ 2.1 GHz, not isolated.
All arithmetic is double precision without outward rounding.

## Final status (2026-09-30)

| Variant | Trials | Reached | Slack / zero-input steps | Min. saved-state bound | Pooled median/p95 filter time | Max rows |
|---|---:|---:|---:|---:|---:|---:|
| 12 mm, unit lift, slack (historical `native_patch_evaluation`) | 30 | 6 | 1497 (2 trials) | 15.3 mm | 10.5/16.7 ms | 2064 |
| `free/` | 30 | 8 | 0 | 12.5 mm | 11.5/42.0 ms | 9045 |
| `sampled_zero/` | 30 | 7 | 0 | 12.3 mm | 10.8/38.9 ms | 16,857 |
| `sampled_zero_6mm/` | 28 | 13 | 0 | 4.3 mm | 29.4/733.8 ms | 70,425 |

`sampled_zero_6mm/` stops at 28 trials by decision: trials 16 and 19 were stopped after about 1.7 h and 0.8 h
(steps 870 and 408 of 1000) because profiling showed 92--100% of the time in the dense NNLS least-distance QP
with 1260--1490 active boxes (one multiplier variable each). No method reaches these two trials (historical
subdivided, unit lift, free, sampled 12 mm), so the reach counts of all other variants are unchanged on the
same 28 trials (subdivided: 17). No partial output of the stopped trials is stored.
On trials 0--4 rerun on this machine, pooled median/p95 filter time is 8.2/14.4 ms (unit lift), 13.0/42.5 ms (free),
and 12.5/42.7 ms (sampled, 12 mm). Trial 18, which needed slack at all 1000 steps with the unit lift, reaches in
4.99 s with optimized multipliers.
Planned before rerunning: a sparse QP solver, multipliers shared per cluster, and local instead of global
gradient bounds in the activation threshold.
