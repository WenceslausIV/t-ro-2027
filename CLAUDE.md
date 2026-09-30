# Project notes for Claude Code

- Answer the user in Korean.
- Paper: `tro/main.tex` (IEEE T-RO, double-anonymous review since Jan 2025, target <= 12 pages).
  Certified geometry-aware CBFs: first-contact lemma, certified B-spline SDF levels, surface-cover
  (lifted-barrier) Bernstein/vertex constraints.
- Hard constraints from the user:
  - The safe set must stay continuous (H = {phi_B(y(x)) >= l_B for all x on S_A}); never discretize
    it. Conservatism is allowed only on the inputs.
  - Guarantees hold on the continuous surface, never on samples.
  - Do not change figure visuals (color theme, style) unless asked.
  - Keep SDF fields; no scene-specific special treatment.
- Core code: `prototype_3d/summed.py`, `prototype_3d/summed3d.py`; experiments in `prototype_3d/`
  (Franka, dual arm, star tube) and the repo root (planar docking, five robots, maze).
- Review history: `review_claude_codex.md` (items C-01..C-35, with a status summary),
  `tro/claude_review_dialogue.md`, audits in `results/manuscript_review/`.
- `prototype_3d/RDF` (Franka meshes/kinematics, github.com/idiap/RDF) is a git submodule; after
  cloning run `git submodule update --init` before running the Franka/dual-arm/star scripts.
- Excluded from git: `*.gif`, `cache/dock_craft_T.npz` (> 100 MB), LaTeX build outputs.
