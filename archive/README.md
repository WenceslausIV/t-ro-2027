# Archive

These files are not used by the current manuscript (`tro/main.tex`) or by any script it depends on.
They were moved from the repository root on 2026-09-30.

- `scripts/`: older experiment and figure scripts. This includes the three-baseline comparison, the maze variants, the four-robot
  search, the level-set figure, and earlier copies of `make_paper_figs.py`. They import root modules (e.g. `sdf_cbf_utils`), so
  run them from the repository root with the root on the path:
  `PYTHONPATH=. python archive/scripts/<name>.py`.
  Some READMEs under `results/` (`three_baselines/`, `maze_workspace_levelset/`) still give the old root paths.
- `notes/`: planning and progress notes (`NOTES_paper_ideas.md`, `PLAN_cover_first.md`, `rework_progress.md`).
- `outputs/`: old GIFs, PNGs, and logs from the configuration-space experiments.
- `swarm/`, `side_notes/`: an earlier swarm/B-spline prototype and side notes.

Reference PDFs (the ICRA 2026 arXiv paper and the ICRA 2027 manuscript) are in `../refs/`.
