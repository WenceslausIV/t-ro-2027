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
