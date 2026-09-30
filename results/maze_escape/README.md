# Irregular branched-wall maze escape (replacement for Fig. 8)

**Superseded:** the current Fig. 8 uses the supplied curved-maze reference.
See [the current experiment](../reference_maze/README.md). This folder retains
the earlier branched-wall results only; `python maze_escape.py` now runs the
reference-maze experiment.

Run `python maze_escape.py` from the repository root. Use `--rebuild` to
regenerate the cached field. A geometry signature prevents reuse of a field
with different walls, robot control points, domain, or grid settings.

The scene has eight **disconnected** walls: four outer wall pieces and four
non-convex branched walls. The internal walls are an elongated asymmetric cross,
a left-facing T, a downward T, and a smaller upward T. Their lengths, branch
positions, and thicknesses differ; beam thicknesses range from 0.40 to 0.60 m.
Each connected wall is a union of one or two overlapping capsule beams (12 beams
in total). No beams from different walls touch. The robot is cached random shape
0, as in the original static-obstacle study. It starts at the lower left, rounds
the cross and passes between the T-shaped walls, then leaves at the upper right.
One orientation-dependent C-space field and one scalar CBF constraint encode
the entire obstacle union, including all branches.

The field domain is `[-6.2,6.2]^2 x S1`, with `160 x 160 x 48` cells.
Fitting uses 96 angular slices, a 2-cm raster, and a `310 x 310` spatial grid.
Control uses a 10-ms step, speed limits 0.8 m/s and 1.4 rad/s, and gain 5.

The route is supplied as seven poses (including the initial pose), not discovered
by the CBF. Nominal translation tracks the active waypoint and nominal rotation
tracks its orientation. Thus this is a waypoint-guided navigation and collision
avoidance demonstration, not a global planner or a deadlock benchmark.

## Geometry and safety checks

- Physical walls are unions of exact capsules. Cubic Bezier quarter-circle segments
  enclose their circular ends; radial polynomial extrema and monotone angular
  traversal are checked in floating point before fitting.
- Contact certification includes all enclosing beam curves, including covered
  pieces inside a branched wall. This is conservative and still includes every
  exposed union boundary; internal overlaps do not create missing contact cases.
- Group connectivity and pairwise separation of distinct wall groups are checked.
- The cached robot ground-truth polygon and its fitted enclosing boundary are
  inherited from the existing five-shape experiment.
- The contact-level branch-and-bound, field-boundary band, and regularity
  checks must all pass before a field cache is accepted.
- Physical clearance is independently measured between the robot polygon and
  exact capsules, including boundary crossings and containment.
- Each integration interval is audited using the rigid-motion bound
  `||delta_position|| + robot_radius * |delta_angle|`. Endpoint distances give
  the interval lower bound `min(d0,d1) - motion/2`. Inconclusive intervals are
  recursively subdivided; unresolved intervals are explicitly counted.
- Paper results are exported only if the robot reaches the exit with no
  collisions, unresolved intervals, QP slack, or negative sampled barrier.
- These are float64 numerical checks, not outward-rounded interval arithmetic.

## Files

- `results.json`: measured results, field verification, geometry, and source hash.
- `trajectory.npz`: every pose/control, each wall clearance, barrier, and timing.
- `maze_escape.png`: paper trajectory figure; also in `tro/figs` as PNG and PDF.
- Dark navy curves show the actual fitted robot boundary over gray ground truth;
  dark wall curves show the exposed cubic enclosing boundary, with internal beam
  seams omitted. Robot outlines are placed every 1.65 m of arc length.
- `layout_preview.png`: geometry/waypoint preview only; generate with
  `python maze_escape.py --layout-only`. Its dashed route is not simulation data.
- `clearance.png`: physical clearance and barrier histories for auditing.
- `field.npz` and `field.json`: reusable fitted field and certification metadata.
- `tro/maze_results.tex`: generated numerical macros used by the manuscript.

The older `results/static.json` obstacle-count study is retained and is not
relabeled as maze data. Its timings and certified levels describe that separate
experiment. The maze is a new demonstration in the same manuscript subsection.
