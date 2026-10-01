# Continuous-boundary baselines vs. ours (Franka, 30 fixed trials) — handoff plan

Goal (user request, 2026-10-01): compare our certified surface-cover filter with similar methods that also
guarantee safety on a continuous boundary, under the same real-time criterion (VLA-style 20-Hz loop,
50-ms budget; our filter itself runs at dt = 10 ms). Work one step at a time; commit after each step so another
agent can continue. Answer the user in Korean (CLAUDE.md).

Common setup for every method (same as results/certificate_upgrades):
- trials: prototype_3d/franka_trials.json (30), nominal u = clip(KQ (qg - q), +-1 rad/s), dt = 10 ms, gamma = 5,
  1000 steps, goal tolerance |q - qg| < 0.05, activation 3 cm.
- obstacles: our certified obstacle fields phi_O with levels l_O (franka3d.build); audit = franka3d.real_gap
  (certified mesh-distance lower bound) at every state.
- QP: summed.solve with SUMMED_QP=daqp (exact check of every applied input), slack fallback reported.

Baselines (all certified on a continuous set that encloses the link mesh; rows enforced at the 10-ms samples):
- B1 `points_delta` (3D Poisson-safety-function style, arXiv 2604.21189): vertices of the link mesh remeshed to
  max edge e; every mesh point is within delta = e/sqrt(3) of a vertex; row h = phi_O(x) - l_O - G_loc delta.
- B2 `capsule` (capsule/line-segment links, e.g. arXiv 2507.01705): one capsule per link enclosing all mesh
  vertices (segment on the principal axis, radius = max vertex distance); segment sampled at spacing s,
  h = phi_O(x_j) - l_O - G_loc (r + s/2).
- B3 `spheres_enclosing` (sphere decomposition, certified): RDF sphere centers, each remeshed vertex assigned to
  its nearest center, radius = max assigned distance + delta.
- B3b `spheres_kmeans16`, `spheres_kmeans64`: same enclosure with 16 / 64 k-means centers per link.
Implementation choice (favors the baselines): baselines use the EXACT obstacle SDF d_O (union of boxes,
1-Lipschitz, no fitting error and no level), h = d_O(x) - m. Ours uses the fitted fields with certified levels.
Ours (existing results): results/certificate_upgrades/{free, sampled_zero_fixed, sampled_zero_6mm_fixed}.

Steps (status):
0. [done] stop other experiments, write this plan.
1. [done] prototype_3d/continuous_baselines.py written; smoke test trial 2, 100 steps, all OK:
   points_delta 97,835 points, m = 5.8 mm, 15 ms median; capsule m = 60-114 mm, 2.8 ms;
   spheres_enclosing (54) r = 44-92 mm, 2.8 ms; spheres_kmeans16 r <= 40-68 mm, 2.7 ms.
2. [running (audit every 10th state; parallel outcomes, then isolated timing) via `bash tro/results/continuous_baselines/run_all.sh` (restartable; commits after each
   method)] run B1-B3 on 30 trials, one process, sequential (isolated timing):
   `python prototype_3d/continuous_baselines.py --methods capsule spheres_enclosing spheres_kmeans points_delta`
   then `python prototype_3d/continuous_baselines.py --methods spheres_kmeans --spheres 64`
   -> results/continuous_baselines/<method>/franka_XX.json; commit.
3. [queued in run_all.sh] isolated timing of ours: `python prototype_3d/patch_size_timing.py --fields 12mm 6mm`; commit.
4. [todo] report table (reached, collisions, min/median gap, stop gap, time median/p95, rows) in
   results/continuous_baselines/README.md; answer the user.

Interim finding (2026-10-01, outcome runs, parallel timing): with the EXACT obstacle SDF, the certified
enclosing capsule/sphere baselines reach 16-22/30 goals, collision-free, ~3 ms/step, against 7-10/30 for ours
(fitted obstacle fields). The obstacle information differs (exact 1-Lipschitz SDF vs fitted field with level and
gradient bound up to 2.68), so this is not apples-to-apples.
5. [todo] fair variant: baselines on OUR fitted obstacle field, h = phi_O(x) - l_O - G_loc m
   (`--obstacle fitted`, folders <method>_fitted). Run capsule, spheres_enclosing, spheres_kmeans 16/64.
6. [todo] report: `python prototype_3d/continuous_baselines_report.py` (step 4 script, written).
