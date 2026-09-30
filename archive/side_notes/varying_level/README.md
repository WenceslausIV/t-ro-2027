# Side note: spatially varying (per-patch) level λ(q)

Not in the paper. Idea: replace the single certified level l* by a smooth level λ(q) built from the
per-cell certified levels ℓ_c, to reduce conservatism where the field is accurate.

## Construction (sound, stays C²)
- λ lives in the field's own cubic B-spline basis (same knots), so ψ = φ − λ is again a B-spline
  (coefficients W − Λ) and h = ψ − l_ψ is C². No extra dλ/dt term: λ depends on q only.
- ℓ_c from `Certifier.certify_cells` (branch-and-bound with per-cell pruning), clipped from below at
  `FLOOR` (always safe, since raising ℓ_c keeps λ ≥ ℓ_c).
- Λ: spread ℓ_c to the controls that are interior to each cell, then one correction pass raising
  all 64 controls of every cell whose Bernstein minimum is still below ℓ_c (convex hull ⇒ λ ≥ ℓ_c
  on every cell).
- ψ is then certified with the paper's Algorithm 1 (plus band and regularity), so safety never
  depends on the construction being tight.

## Result on the long tapered docking (78×78×144 B-spline, 9.4×10⁵ coefficients)
| | nose depth (full 0.58 m) | final true gap | min h |
|---|---|---|---|
| single level l* = 12.5 mm | 0.509 m | 37.2 mm | 0 |
| per-patch λ(q), floor 0 | 0.519 m | 34.8 mm | 0 |

- Per-cell certificate: 262,740 cells meet M, 22 min on the RTX 4070; ℓ_c from −1.20 m to +12.5 mm.
- Without clipping, λ exploded (deep-interior cells at −1.2 m next to boundary cells at +10 mm) and
  the nose stopped at 0.139 m. Floors 0 / −10 / −30 mm: depths 0.519 / 0.502 / 0.480 m.
- Why the gain is small here: the cells along the docking path already have ℓ_c ≈ 4–8 mm
  (8–11 mm in their neighborhoods), and the fitted field is already −5 mm at full seating, i.e. the
  0.5-cm raster training targets, not the level, limit the insertion. The global max ℓ_c sits at an
  unrelated configuration (θ ≈ 271°).
- Likely more useful where the worst fit error is far from where the task needs to go, or together
  with exact training targets.

## Files
- `dock_varying.py` — builds/certifies λ, runs both closed loops, writes the GIF.
- `dock_varying.gif`, `dock_varying_final.png`, `result.json`.
- Caches: `cache/dock_craft_lc.npy` (per-cell levels), `cache/dock_craft_varying_floor+0mm.npz`.
