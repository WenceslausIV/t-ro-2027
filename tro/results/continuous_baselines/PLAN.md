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
5. [running via `bash tro/results/continuous_baselines/run_fitted.sh`] fair variant: baselines on OUR fitted obstacle field, h = phi_O(x) - l_O - G_loc m
   (`--obstacle fitted`, folders <method>_fitted). Run capsule, spheres_enclosing, spheres_kmeans 16/64.
6. [todo] report: `python prototype_3d/continuous_baselines_report.py` (step 4 script, written).
NOTE: run_all.sh step 2b (timing) started while run_fitted.sh was still running, so the *_timing results and
step 3 may be contaminated by concurrent load. After everything finishes: check `pgrep -f continuous_baselines`
is empty, delete results/continuous_baselines/*_timing and results/patch_size_timing, and rerun only the timing
part of run_all.sh (the loops after "step 2b") alone.

Status 2026-10-01 (later):
- Steps 2a and 5 done (all outcome runs: 5 exact-SDF baselines, 4 fitted-field baselines).
- Finding: fitted-field baselines never get closer than their start pose (min gap = initial gap), and several
  baselines start OUTSIDE their certified set (h < 0 at q0, so their guarantee does not apply in those trials):
  capsule exact 7/30, capsule fitted 12/30, spheres RDF fitted 8/30, spheres k-means64 fitted 3/30.
  Valid from the start: spheres exact (RDF, k-means16, k-means64) and points_delta exact.
- points_delta exact (PSF style) reaches 22/30 with min/median gap 3.2/4.4 mm: better than ours on both.
  Likely cause on our side: the velocity remainder rows penalize the full twist |V|, |Omega| of a link near an
  obstacle (direction-independent), so tangential sliding is slowed; the point/sphere CBFs only constrain the
  normal component.
- The old *_timing folders were deleted (overlapped with other runs).
7. [done] points_delta on our fitted fields: 22/30 reached, 0 collisions, min/median gap 9.2/12.8 mm, starts
   valid in 30/30 (h >= 0 at q0). Ours 6 mm sampled: 10/30, 3.8/10.4 mm; ours 12 mm sampled: 7/30, 12.2/20.9 mm.
8. [running in run_final.sh] isolated timing of all baselines and ours, nothing else running.
9. [todo] `python prototype_3d/continuous_baselines_report.py` (add points_delta_fitted row) and answer the user.

## HANDOFF (2026-10-01, end of this session)
State: all outcome runs done; isolated timing done for all baselines (*_timing/, 30/30 each).
Table: results/continuous_baselines/README.md and summary.json (`python prototype_3d/continuous_baselines_report.py`).
Not done / next steps for the next agent:
a. Isolated timing of OURS (nothing else running): `python prototype_3d/patch_size_timing.py --fields 12mm 6mm`
   (writes results/patch_size_timing/, picked up by the report). Our rows in the table still use parallel-load times.
b. Fairness gap 1 (most important): the compared ours rows use the sampled-data certificate (extra caps and
   tightening), baselines use rows at the samples only. Run ours 6 mm WITHOUT sampled data for a like-for-like row:
   `python prototype_3d/certificate_upgrades.py --variant free --field 6mm --qp daqp --tag fixed --trials 0 30`
   (12 mm free exists: results/certificate_upgrades/free, 8/30) and add both to METHODS in the report script.
c. Fairness gap 2: baselines audited every 10th state (+ final), ours every state; baseline min gaps may be
   slightly optimistic.
d. Fairness gap 3: baseline parameters untuned (points edge 10 mm, 1 capsule/link, 16/64 spheres).
e. Key finding to address in theory: fair comparison (same fitted obstacle fields, valid starts):
   points+delta 22/30 reached (gap 9.2/12.8 mm, 48/77 ms isolated) vs ours 6 mm 10/30 (3.8/10.4 mm) and
   12 mm 7/30. Hypothesis: our velocity-remainder rows penalize the whole link twist near obstacles
   (direction-independent), so tangential sliding slows; point/sphere CBFs constrain only the normal component.
   A direction-aware remainder bound is the candidate fix.
Other agent's work: origin/main has tro/paper/concerns.md (non-smooth feedback from box activation); merged here.
Also see results/certificate_upgrades/README.md (solver/cap bugs fixed and audit), results/mobile_patch_sweep/,
results/native_patch_sizes/ (18/24-mm refits), prototype_3d/patch_size_timing.py.

Handoff steps a and b launched 2026-10-01 via `bash tro/results/continuous_baselines/run_handoff.sh` (restartable;
results in results/patch_size_timing/ and results/certificate_upgrades/free_6mm_fixed/). Rerun the same command if it stopped.
