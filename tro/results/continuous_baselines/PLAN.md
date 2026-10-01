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
  its nearest center, radius = max assigned distance + delta; h = phi_O(c) - l_O - G_loc r.
G_loc: franka3d.neighborhood_bounds(phi_O, reach) cellwise gradient bound.
Ours (existing results): results/certificate_upgrades/{free, sampled_zero_fixed, sampled_zero_6mm_fixed}.

Steps (status):
0. [done] stop other experiments, write this plan.
1. [todo] write prototype_3d/continuous_baselines.py, smoke test 1 trial / 100 steps, commit.
2. [todo] run B1-B3 on 30 trials, one process, sequential (isolated timing):
   `python prototype_3d/continuous_baselines.py --methods points_delta capsule spheres_enclosing`
   -> results/continuous_baselines/<method>/franka_XX.json; commit.
3. [todo] isolated timing of ours: `python prototype_3d/patch_size_timing.py --fields 12mm 6mm`; commit.
4. [todo] report table (reached, collisions, min/median gap, stop gap, time median/p95, rows) in
   results/continuous_baselines/README.md; answer the user.
