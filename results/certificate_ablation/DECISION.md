# Implementation and manuscript decision — 2026-09-29

## Completed

- Rounded star-tube geometry retained, including 15-mm cross-section fillets and
  10-mm extrusion-edge rounding. The existing figure renderer and visual style
  are retained.
- Six certificate variants implemented: vertex/joint Bernstein, each with unit,
  regularized norm-ratio, or regularized projection multiplier.
- Negative multipliers use the absolute multiplier in the Taylor error bound.
  Screening and refinement still use the unit lift, independently of the chosen
  multiplier. No spline, level, or continuous safe set is changed.
- 24 pilot runs completed: six docking, six rounded-tube, and twelve Franka runs
  (two certificates on six selected trials). All runs required zero slack. The
  geometry hash is identical across variants of each scene.
- Validation passed: 29,376 numerical residual comparisons, an adversarial
  negative-multiplier regression, and exact equality of unit-weight rows with
  the saved original implementation in both 2D and 3D. Numerical comparisons
  test the implementation; the continuum argument is in the manuscript.
- An independent geometric audit of all 1,501 saved states of the original
  15-second rounded-tube run certifies a lower bound of 1.000127 mm, with no
  unresolved states. It covers all modeled link triangles, uses no fitted-field
  cutoff, and transfers valid link bounds using a rigid-motion displacement
  bound when possible (193 direct link checks, 11,815 transfers). This is a
  conservative lower bound, not an estimate of the true minimum distance.
  The result is in `../../prototype_3d/startube_setup_ours.audit.json`.

## Decision

Keep `vertex` with unit multiplier as the default. The joint Bernstein condition
is an implemented alternative and is presented before the vertex reduction in
the method section. It preserves value/velocity cancellation but costs more.
Fixed multipliers remain experimental options, not a claimed general gain.

On the six-second rounded-tube comparison, joint Bernstein with unit multiplier
reduces input modification by approximately 4% (1.948 to 1.870) and arrival time
from 3.99 to 3.93 s. Maximum rows increase from 26,440 to 89,424. Fixed multipliers
do not improve the vertex variant in this scene.

Both unit-weight methods reach in four of the six selected Franka trials. Joint
Bernstein reduces one arrival time from 7.38 to 6.63 s but increases filter time.
The most expensive trial still stalls: median time rises from 115 to 253 ms.
This does not justify replacing all existing paper results with a new default.

## Refinement and reproducibility

The new pilot consistently uses `SUMMED_REFINE_NEAR=inf`, theta=0.1, depth=2.
This removes the distance gate; it does not split every active box unconditionally.
Original planar, 30-trial Franka, and 20-trial dual-arm results were preserved.
The production legacy default remains near=0.01; the successful original rounded
tube uses near=inf. This difference is disclosed in the manuscript, alongside
the new comparison under a common rule. The original full experiment suite has
not been rerun under the common rule.

The 15-second tube run and six-second comparison have different distance minima
and timing distributions. Do not combine their statistics. Local pilot timings
use one BLAS thread, but an existing unrelated process was running; they are not
isolated timing results and do not replace the paper's isolated measurements.

## Scope of the claims

The paper now separates the continuous safe set from finite input certificates,
includes strict-feasibility convergence under arbitrary refinement, and states
that the finite-depth implementation need not produce nested input sets.
It does not claim general second-order conservatism or that multipliers remove
the cover boxes. The guarantee assumes continuous time and exact arithmetic.
Saved-state distance checks alone do not prove inter-step safety.

Remaining broader work includes an all-time dual-arm distance audit (the old
routine selects candidate states), isolated timing of new variants, and a full
experiment rerun if a uniform refinement policy is later promoted. The present
comparison provides no new 30-trial success-rate estimate or hardware result.

## Files

- `REPORT.md`: all 24 raw metric summaries.
- `pilot/*.json`: configuration, geometry/source hashes, and metrics.
- `pilot/*.npz`: trajectories and both applied/nominal inputs.
- `baseline/`: pre-edit sources and star figure.
- `../../tro/main.tex`, `../../tro/certificate_comparison.tex`: revised text.
- `../../tro/build_certificate/main.pdf`: compiled manuscript.

See `README.md` for commands and `prototype_3d/certificate_audit_star.py` for
the separate all-saved-state geometric star audit.
