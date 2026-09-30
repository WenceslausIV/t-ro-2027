# Project notes for Claude Code

- Answer the user in Korean.
- The repository holds two unrelated projects: `tro/` (the T-RO paper and all its code and data) and `swarm/`
  (an earlier swarm / B-spline prototype).
- Paper: `tro/paper/main.tex` (IEEE T-RO, double-anonymous review since Jan 2025, target <= 12 pages).
  Certified geometry-aware CBFs: first-contact lemma, certified B-spline SDF levels, surface-cover
  (lifted-barrier) Bernstein/vertex constraints.
- Hard constraints from the user:
  - The safe set must stay continuous (H = {phi_B(y(x)) >= l_B for all x on S_A}); never discretize
    it. Conservatism is allowed only on the inputs.
  - Guarantees hold on the continuous surface, never on samples.
  - Do not change figure visuals (color theme, style) unless asked.
  - Keep SDF fields; no scene-specific special treatment.

## Layout of `tro/` (run all scripts from `tro/`)

- `paper/`: manuscript, bibliography, `figs/`, reference PDFs (`paper/refs/`), and the review record
  (`paper/review_claude_codex.md`, items C-01..C-35 with a status summary; `paper/claude_review_dialogue.md`).
- `planar/`: planar experiments and paper figures (docking, five robots, maze, Figs. 1-5),
  e.g. `python planar/make_paper_figs.py dock2`.
- `prototype_3d/`: the surface-cover barrier core (`summed.py`, `summed3d.py`), the Franka / dual-arm /
  star-tube experiments, and the planar surface-cover runners (`dock_cover.py`, `five_cover.py`).
  `prototype_3d/RDF` (Franka meshes/kinematics, github.com/idiap/RDF) is a git submodule; after cloning run
  `git submodule update --init`.
- `results/`: result JSON/NPZ and audit scripts (`results/manuscript_review/`).
- `cache/`: precomputed fields and shapes (npz) that the experiments load.
- `archive/`: old notes and outputs.
- Scripts locate `tro/` from their own file location, so moving one script out of its folder breaks its paths.
- Excluded from git: `*.gif`, `tro/cache/dock_craft_T.npz` (> 100 MB), LaTeX build outputs, logs.
