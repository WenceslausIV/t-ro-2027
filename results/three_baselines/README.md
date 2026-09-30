# Three-method shared-geometry experiment

Run from the repository root:

```powershell
python compare_three_baselines.py --check-only
python compare_three_baselines.py --trials 6 --seconds 12 --dt 0.01
```

The second command regenerates `comparison.json`, `comparison_trajectories.npz`,
`table.md`, `trajectories.png`, `tro/figs/three_baselines.png`, and
`tro/three_baseline_table.tex`. Timings are measured, not copied from the papers.
The JSON stores versions, protocol, seed, controller adaptations, and script hash.

## Scope

This is a **fixed-orientation, one-obstacle, common C-space representation
benchmark**, not a reproduction of either publication's experimental setup or a
test of rotating multi-robot navigation. All generated runs are retained. There
is no deadlock detector, recovery controller, success-rate ranking, timeout
ranking, or goal-time statistic. All runs use the full 12-second horizon.

The translating robot is the union of rectangles with half extents
`(.26,.09)` and `(.09,.26)` m. The static obstacle is the union of rectangles
`(.50,.14)` and `(.14,.50)` m. Both are nonconvex crosses. At fixed orientation,
their exact C-space obstacle is the union of four centered rectangles obtained
by adding every robot/obstacle pair's half extents. Therefore the physical
body-to-body set distance is available exactly when they are disjoint.

Every method receives the SAME smooth star-shaped envelope of this C-space
obstacle, enlarged to include a 20-mm Euclidean clearance. Each rectangle is
first dilated by a 20-mm square (a conservative disk dilation), then the maximum
radial function is interpolated by a periodic cubic spline. A Lipschitz bound
between verification nodes supplies the necessary outward radial shift; it is
not just a sampled containment assertion. All arithmetic is float64, without
outward rounding. The extra radial shift is reported in JSON. Common geometry
preparation is kept separate from method-specific preparation.

Starts are `(-2.8,y)`, goals `(2.8,y)`, with offsets centered at
`[-.65,-.50,-.35,.35,.50,.65]` m and seeded jitter of +/-15 mm. All methods use
the same nominal proportional velocity, a 0.6-m/s inscribed 12-sided speed
polygon, gain 5, 10-ms timestep and a shared workspace radius of 4.2 m. These
choices are fixed across methods, not tuned separately to favor a controller.

## Method definitions and departures from source papers

1. **Ours.** Fit a uniform bicubic B-spline distance surrogate to the common
   envelope. Bound its maximum on the *continuous* envelope boundary with
   radial-spline interval boxes and Bernstein subdivision (1-mm optimization
   tolerance). Check activation band and level regularity. The orientation
   dimension is constant and is algebraically collapsed for online evaluation.
   Its values/gradients are checked against the original tensor-product evaluator.
   This is a translation-only specialization, not a measurement of full SE(2)
   field performance. Contact certification is a continuous-boundary check;
   fitting uses approximate distance samples.

2. **Sampling-based (adapted), N=64/256/1024.** Based on Lutkus, Chong and
   Lindemann, <https://arxiv.org/abs/2504.09038>. Apply their sampled-distance
   construction to the common C-space envelope, with a point representing the
   finite robot's reference translation. The covering radius is
   `eta = pi/N * sup ||boundary'(theta)|| + 1e-6` m. Use individual barriers
   `||x-s||^2-eta^2`; their minimum is the sampled barrier. Enforce every point
   row that is not automatically redundant for all permitted velocities. At
   speed limit V and gain gamma, distances above
   `(V+sqrt(V^2+gamma^2*eta^2))/gamma` are redundant. A KD-tree finds the remaining
   rows. This sufficient all-sample formulation covers all active minimizers
   and is stronger than constraining only instantaneous minima. The Euclidean
   covering argument is used directly rather than treating squared distance as
   a metric. It is NOT the repository's single-closest-pair baseline.

3. **Ball world (adapted).** Based on Notomista and Saveriano,
   <https://arxiv.org/abs/2106.06330>. Use their deforming-ball strategy in the
   single-obstacle case with the explicit radial map
   `F(x)=q_c+rho*(x-c)/r(angle(x-c))`. This maps the star obstacle exterior to
   the ball exterior and has `det J=rho^2/r^2>0`. The mapped center AND radius
   are controlled through the state-avoidance QP, with nominal restoring rates;
   the deformation mechanism has NOT been removed. A minimum mapped radius of
   0.1 prevents degeneracy. Use
   `J*u + F_parameters*parameters_dot = q_dot` to account for changing mapping
   parameters. Scale physical and virtual rates together to honor the common
   speed polygon/workspace bound. This differs from a literal discretization
   of Algorithm 1 and does not reproduce its multi-obstacle blended map or its
   compact ball-world boundary. A two-obstacle blended-map prototype failed
   Jacobian checks and its runs were NOT used as evidence about the paper.

The parameter choices are recorded in code/JSON. These are transparent,
independent reimplementations/specializations, not author-provided code.

## Measurements and verification

- **Physical gap:** exact disjoint rectangle-union distance, using original
  finite bodies; not the controller's own barrier value. Negative indicators
  would signify overlap, not minimum translational penetration depth.
- **Intersample check:** for each linear translation segment, use the
  1-Lipschitz distance bound with exact endpoint queries. Subdivide if needed;
  retain an explicit unresolved count instead of silently declaring a small
  interval safe. Store both the minimum evaluated gap and a conservative
  interval-wide lower bound. Collision checks and logging are outside timers.
- **Input TV:** `sum(norm(v[k+1]-v[k]))`, in m/s, averaged across fixed-horizon
  runs. It is not normalized by path length and is not a goal-progress metric.
- **Closed-loop online time:** wall-clock milliseconds for complete controller
  computation, including QPs, KD queries, mapping/Jacobian and virtual updates.
  All online methods use CPU NumPy/SciPy. Offline field bounds use PyTorch/GPU.
- **Matched-pose online time:** identical states from every coarse-sampling run
  at a fixed 20-step stride, repeated three times. Reset mapped parameters
  outside the timer, then include the complete online update. This removes
  differences in sampled physical states due to different paths or stopping;
  it is a query-cost diagnostic, not a trajectory replay guarantee.
- Source methods have different barrier scalings. Equal gain does not imply
  identical physical approach behavior. Compare physical gaps, not raw h values.
- Preparation times and coefficient/coordinate storage are saved separately.
  KD-tree internal allocations and Python object overhead are not included in
  array-storage figures. Runtime comparisons are implementation-specific.

Checks cover mapping state/parameter derivatives against central differences,
multiple simultaneous sample constraints, detected between-endpoint crossing,
radial containment, dense contact-level sanity checks, and exact collapse of the
constant orientation axis. These do not constitute verified-rounding proofs.

This small benchmark cannot establish general dominance, angular handling,
multi-obstacle scaling, or deadlock performance. In particular, a method that
stops can have small total variation; do not interpret that statistic as a
navigation-success claim.
