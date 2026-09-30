# Certificate comparison (2026-09-29)

The physical geometry is fixed. In particular, the tube is the **rounded** star tube;
the experiment does not replace it with the old sharp shape. Fields, levels and
initial covers are shared by all variants within a scene and hashed in each result.

## Commands

```powershell
python prototype_3d/certificate_check.py
python prototype_3d/certificate_ablation.py --scenes dock star
python prototype_3d/certificate_ablation.py --scenes franka --variants vertex-one bernstein-one
```

The runner forces OMP/MKL/OpenBLAS to one thread and runs variants sequentially.
Do not run other heavy work concurrently when measuring timings. It saves each
run independently and skips completed filenames on restart. Use a new `--out`
directory for a different horizon, refinement threshold or source version.
Original paper result files are never overwritten by this runner.

## Methods and invariants

- `vertex` uses 4/8 rows per 2D/3D cell; `bernstein` uses 9/27 joint coefficients.
- `one` fixes lambda=gamma.
- `norm` sets w=lambda/gamma to the ratio of regularized gradient norms.
- `proj` sets w=-a.b/(a.a+1e-12), permitting negative weights.
- Taylor value errors use **abs(w)** times the A-side third-derivative bound.
- Screening, activation and refinement use the original unit lift for every
  variant, so the cells at any fixed pose are independent of the multiplier.
- The comparison uses `near=inf`: active cells are refined only if their gradient
  remainder exceeds theta=0.1, with at most two halvings. This is not unconditional
  refinement of every active cell. The legacy default remains reproducible.

The multiplier is part of a pointwise certificate of F>=0 on g_A=0. It is not a new
barrier whose time derivative is used. It neither changes the safe set nor removes
the boxes. General second-order accuracy is not claimed.

## Interpretation

- Recorded `min_h_mm` is the **unit-lift screening lower bound**, including for
  weighted variants. It may be negative while a weighted certificate is valid;
  do not compare it to zero as a test of the weighted constraints.
- A negative distance lower bound means the distance check is inconclusive. The
  `uncertified_gap_states` count is not a count of actual collisions.
- The existing 3D distance routines are evaluated at discrete states and use
  screening. This pilot is not an audit of all geometry and time intervals.
- Time per filter step is wall-clock time; reached_s is simulated time with a
  10-ms integration step. A 77-ms filter is not demonstrated at 100 Hz.
- Input modification integrates the Euclidean norm in the coordinates of each
  scene. Planar translation and rotation have mixed units; compare variants
  within a scene, not raw integrals across different scenes.
- Docking's nominal target need not satisfy the certified clearance. Report
  achieved seating/clearance as well as the generic goal tolerance flag.

`baseline/` preserves the files before this work. `pilot/` contains JSON metadata,
metrics and NPZ trajectories/commands. Final observations are recorded in
`REPORT.md`; `DECISION.md` records the implementation decision and remaining scope.

## Independent rounded-tube audit

```powershell
python prototype_3d/certificate_audit_star.py prototype_3d/startube_setup_ours.npz
```

This checks all 1,501 saved poses and every modeled link triangle with analytic
distance bounds, geometric AABB pruning, and conservative rigid-motion transfer
of previous link bounds. No fitted-field cutoff is used. The certified lower
bound is 1.000127 mm with zero unresolved states. Transfers are accepted above
1 mm; the resulting bound is not an estimate of the actual minimum separation.
Use `--reuse-margin-mm inf` to evaluate every link directly at every saved pose.
This audit does not certify the motion between saved poses.
