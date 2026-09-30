# Rocket escape through the curved maze

Run `python reference_maze_escape.py` (or `python maze_escape.py`).
The paper exports are `tro/figs/rocket_maze.png` and `.pdf`.

The user-supplied rocket is manually reconstructed as a planar body in
`rocket_geometry.py`. Its silhouette includes the nose, fuselage, and two side
fins; exhaust is excluded. This is a geometric reconstruction, not exact raster
segmentation. Heading zero points the nose along +x. The physical polygon is
about 0.477 m long and 0.261 m wide. A 64-control-point periodic cubic B-spline
encloses it, with enclosure checked by subdivision and winding. Its outward
fitting offset is 3 mm. The window and color patches are clipped to the same
physical silhouette used for collision checking.

The walls use the original manually traced cubic ground truth and its separately
fitted, enclosing nonuniform cubic B-splines. Their maximum boundary-distance
upper bounds remain below 1.47 mm. The figure shows gray wall interiors and red
fitted wall outlines. The blue line is the simulated rocket reference-point
trajectory. No C-space contour or fixed-orientation safety boundary is displayed.

The controller's SE(2) field is retrained for the rocket and fitted wall union.
An offline A* route supplies nominal waypoints. Escape, physical collision
checks, interval rigid-motion clearance bounds, and nonnegative sampled barriers
are required before exporting paper results. Arithmetic is float64 without
outward rounding.

`robot_geometry.npz` stores the actual physical polygon and fitted controls.
`wall_fits.npz` stores the wall coefficients and knots. `field.npz` is the
controller field; `results.json` records measurements, certificates, and source
hashes; `trajectory.npz` stores poses, controls, physical gaps, and barriers.

Earlier star-shaped-robot results in `results/reference_maze` and the separate
workspace-only visualization in `results/maze_workspace_levelset` are retained
as previous experiments. Their measurements are not reused for this rocket.
