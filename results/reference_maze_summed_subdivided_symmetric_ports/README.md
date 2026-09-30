# Symmetric maze entrance and exit

Figure 5 was regenerated from a new simulation after moving the frame upward by 0.35 m. Its bounds are [-4.44, -3.34] to [4.44, 5.54] m. The upper and lower straight passage sections now both measure 1.14 m. The passage centerline and width profile are unchanged. Original scene results remain in their original directories.

Reproduction (PowerShell, from the repository root):

```powershell
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:SUMMED_REFINE_DEPTH='2'
python -u maze_cover.py --side 5 --symmetric-ports
```

The wall SDF and certified level were recomputed for the revised geometry. The run reached the exit in 32.42 s over 3242 steps with no slack, no collision intervals, and no unresolved intervals. The minimum saved-state gap is 4.69356 mm and the audited interval lower bound is 0.594399 mm. Filter median/p95 are 1.9058/3.7056 ms on one CPU thread; background activity was not excluded, so these are not isolated timings. Full measurements and geometry metadata are in results.json.

rocket_maze_cover.png is the source of tro/figs/rocket_maze.png. The image shows the newly simulated trajectory, not a transformed version of the previous trajectory.
