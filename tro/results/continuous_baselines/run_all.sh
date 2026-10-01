#!/bin/bash
# Steps 2-3 of PLAN.md. Restartable (finished trials are skipped); commits after every stage.
#   nohup bash tro/results/continuous_baselines/run_all.sh > /tmp/cb_run.log 2>&1 &
# 2a: outcomes, 4 methods in parallel, mesh audit every 10th state and the final state.
# 2b: timing, one process, no audit (same trajectories), into <method>_timing/.
# 3 : isolated timing of ours (patch_size_timing.py).
cd "$(dirname "$0")/../.."            # -> tro/
commit() {
  cd ..
  git add tro/results && git commit -q -m "$1

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TByxUTGsR8uPvFhDBetfgi" || true
  for d in 2 4 8 16; do git push -q origin claude/tro-2027-paper-review-wg3qzv && break; sleep $d; done
  cd tro
}
CB=prototype_3d/continuous_baselines.py
python $CB --methods capsule --audit-every 10 &
python $CB --methods spheres_enclosing --audit-every 10 &
python $CB --methods spheres_kmeans --audit-every 10 &
python $CB --methods spheres_kmeans --spheres 64 --audit-every 10 &
wait; commit "Continuous-boundary baselines: outcomes of capsule and sphere sets (step 2a)"
python $CB --methods points_delta --audit-every 10 --trials 0 15 &
python $CB --methods points_delta --audit-every 10 --trials 15 30 &
wait; commit "Continuous-boundary baselines: outcomes of points_delta (step 2a)"
for m in capsule spheres_enclosing spheres_kmeans points_delta; do
  python $CB --methods $m --audit-every 0 --tag _timing
done
python $CB --methods spheres_kmeans --spheres 64 --audit-every 0 --tag _timing
commit "Continuous-boundary baselines: isolated timing (step 2b)"
python prototype_3d/patch_size_timing.py --fields 12mm 6mm && commit "Isolated timing of our filter, 12 mm and 6 mm (step 3)"
echo ALL_DONE
