#!/bin/bash
# Steps 2-3 of PLAN.md, one process at a time (isolated timing); commits after every stage.
# Run from anywhere:  nohup bash tro/results/continuous_baselines/run_all.sh > /tmp/cb_run.log 2>&1 &
# Restartable: finished trials are skipped.
cd "$(dirname "$0")/../.."            # -> tro/
commit() {
  cd ..
  git add tro/results && git commit -q -m "$1

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TByxUTGsR8uPvFhDBetfgi" || true
  for d in 2 4 8 16; do git push -q origin claude/tro-2027-paper-review-wg3qzv && break; sleep $d; done
  cd tro
}
for m in capsule spheres_enclosing spheres_kmeans; do
  python prototype_3d/continuous_baselines.py --methods $m && commit "Continuous-boundary baseline results: $m (step 2)"
done
python prototype_3d/continuous_baselines.py --methods spheres_kmeans --spheres 64 && commit "Continuous-boundary baseline results: spheres_kmeans64 (step 2)"
python prototype_3d/continuous_baselines.py --methods points_delta && commit "Continuous-boundary baseline results: points_delta (step 2)"
python prototype_3d/patch_size_timing.py --fields 12mm 6mm && commit "Isolated timing of our filter, 12 mm and 6 mm (step 3)"
echo ALL_DONE
