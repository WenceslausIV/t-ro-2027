#!/bin/bash
# Steps 7-8 of PLAN.md. Restartable; commits after every stage.
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
# 7: points_delta on our fitted obstacle fields (outcomes, audit every 10th state)
python $CB --methods points_delta --obstacle fitted --audit-every 10 --trials 0 15 &
python $CB --methods points_delta --obstacle fitted --audit-every 10 --trials 15 30 &
wait; commit "Continuous-boundary baselines: points_delta on our fitted obstacle fields (step 7)"
# 8: isolated timing, one process at a time, nothing else running
for m in capsule spheres_enclosing spheres_kmeans points_delta; do
  python $CB --methods $m --audit-every 0 --tag _timing
done
python $CB --methods spheres_kmeans --spheres 64 --audit-every 0 --tag _timing
python $CB --methods points_delta --obstacle fitted --audit-every 0 --tag _timing
python prototype_3d/patch_size_timing.py --fields 12mm 6mm
commit "Isolated timing of baselines and ours (step 8)"
echo FINAL_DONE
