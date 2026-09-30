# Reference-image curved maze (current Fig. 8)

Run `python planar/reference_maze_escape.py` (or `python planar/maze_escape.py`).
Use `--layout-only` to preview geometry and the planned route, or `--rebuild`
to regenerate the field. The paper files are `paper/figs/reference_maze_tight.pdf`,
`paper/figs/reference_maze_tight.png`, and the generated metrics `tro/maze_results.tex`.

## Geometry

The user's supplied maze is reconstructed from manually traced image-coordinate
landmarks in `maze_reference_geometry.py`. It preserves the two long opposing
boundaries, winding passage, curled blind branches, upper-left entrance, and
lower-right exit. Rocket, moon, text, and the screenshot border are omitted.
This is a vector reconstruction, not a claim of pixel-exact extraction.

Piecewise cubic Hermite curves through the landmarks define the actual wall
geometry; straight outer edges close the two wall regions. No geometric
guarantee is claimed for the original raster picture. The cached random robot 0
and its enclosing B-spline are both scaled by 0.55 (fitted radius about 0.272 m).
This is 1.447 times the previous robot's linear size (2.095 times its area).

Gray regions represent physical walls; dark outlines are separate periodic cubic
B-splines interpolating 1-mm outward-offset wall samples with rounded convex
joins. Nonuniform knots concentrate near corners (4405 and 4176 control points),
avoiding the previous uniform fit's 8/12-mm global offsets. Navy curves enclose
the robot body.
The paper figure shows workspace body boundaries, as in the docking example;
configuration-space contours are not overlaid.

Enclosure verification recursively subdivides each fitted Bezier span. An
axis-aligned box enclosing its control hull must remain outside the reference
polygon by more than the original wall's chord-error bound. A nonzero winding
number on the accepted ordered chords establishes enclosure of each connected
physical wall. The checks cover the continuous curves rather than only samples.
A continuous fitted-curve-to-physical-boundary distance upper bound is checked
to be below 3 mm (actual computed bounds: 1.464 and 1.466 mm). It combines
point-to-polygon distance upper bounds, a 0.25-mm derivative-based bound between
curve samples, and the reference chord error. These are float64 bounds, not
interval-arithmetic certificates. The fitted wall curves are used for both field training and contact certification;
physical collision checks use the original walls. No displayed offset includes
the robot radius.

Each wall curve is discretized for physical distance evaluation with the bound
`sup ||B''|| / (8*n^2)` per cubic span. The maximum chord deviation is at most
0.1 mm. Reported sampled gaps subtract this bound from polygon distances.
Exact disjoint polygon distances use segment endpoints, explicit crossing and
containment checks, and conservative KD-tree candidate pruning.

## Navigation and safety

- A* on a 3-cm occupancy grid inflated by the fitted robot radius plus 2 cm
  supplies the nominal route. The path stays inside the reference rectangle,
  so it cannot bypass the maze via the outside. Its geometry is resampled into
  pose waypoints; this planner is distinct from the CBF safety filter.
- One `160 x 160 x 48` field over `[-6,6]^2 x S1` represents both wall regions.
  Training uses 96 angular slices, a 2-cm raster, and a `300 x 300` spatial grid.
- Offline wall enclosure, contact, activation-band, and regularity checks must all succeed.
- Control: 10-ms step, 0.8 m/s and 1.4 rad/s input bounds, CBF gain 5.
- Every integration interval is checked with the rigid-motion Lipschitz bound,
  including the wall chord-error correction. Inconclusive intervals are
  subdivided; unresolved intervals are counted explicitly, not declared safe.
- Paper results are exported only after reaching the exit with no collisions,
  unresolved intervals, slack, or negative sampled barrier beyond tolerance.
- Calculations use float64 without outward rounding.

## Evidence

`results.json` records the measurements, field checks, and SHA256 hashes of the
three source modules. `trajectory.npz` stores every pose/control, wall gap lower
bound, barrier, timing, and the planned waypoints. `clearance.png` shows the
clearance and barrier histories. `layout_preview.png` is a planning preview,
not the simulated paper trajectory.

`wall_fits.npz` stores nonuniform B-spline control points, knots, Bezier spans,
and span lengths; `results.json` records
wall enclosure checks and offsets. `tight_wall_corner.png` enlarges the
upper-left physical corner and its smooth enclosing cubic B-spline.
The older `fitted_wall_corner.png`, `displayed_field_slice.npz`, and `field_corner.png` are superseded
diagnostics from the previous display, not used in the current paper figure.

The earlier branched/capsule layout and results in `results/maze_escape` are
superseded. The obstacle-count benchmark `results/static.json` remains separate
and is not relabeled as measurements of this maze.
