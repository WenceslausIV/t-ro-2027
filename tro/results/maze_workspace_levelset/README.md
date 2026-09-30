# Workspace maze enclosure comparison

Run `python maze_levelset_figure.py`.

Outputs: `paper/figs/maze_ground_truth_levelset.png` and `.pdf`.
The left panel is the physical maze in gray. The right panel displays the
sublevel set of a fitted bicubic workspace SDF and its enclosing level contour.
There is no robot footprint, C-space inflation, trajectory, or body-curve
substitution in this figure.

The field has 600 x 600 cells on [-6,6]^2 and uses signed-distance targets from
a 2.5-mm raster, sampled on a 1801 x 1801 fitting grid. A constant angular basis
allows reuse of the existing field evaluator while the polynomial is exactly
independent of orientation.

The smallest enclosing level for this fixed field lies in
[8.813355, 8.913355] mm. The upper endpoint is displayed. Boundary
branch-and-bound gives a lower bound attained on the true cubic wall boundary.
The computed upper level is additionally checked over the entire wall union:
each field box is either outside the walls or has all restricted Bernstein
coefficients below the level. Original cubic wall boundaries are discretized
with a 0.1-micrometer chord-error bound for this exclusion test. Thus interior
holes cannot be silently ignored when claiming enclosure.

These are float64 computations without outward rounding; the enclosure is not
a verified interval-arithmetic result. "Tightest" is qualified by the 0.1-mm
level bracket for this field, not an optimization over all possible fields.

`field.npz` contains the actual field coefficients. `displayed_field.npz`
contains the evaluated field values used for contouring. `results.json` records
the bracket, enclosure/regularity checks, and source hashes.

The rotating-robot escape results and tighter parametrized body enclosures in
`results/reference_maze` remain a separate experiment. Its confirmed 40.44-s
escape is not relabeled as a run using this two-dimensional field.
