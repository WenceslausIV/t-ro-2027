#!/bin/bash
# Step 5 of PLAN.md: baselines on our fitted obstacle fields (same obstacle information as ours). Restartable.
cd "$(dirname "$0")/../.."            # -> tro/
CB=prototype_3d/continuous_baselines.py
python $CB --methods spheres_kmeans --spheres 64 --obstacle fitted --audit-every 10 &
python $CB --methods capsule --obstacle fitted --audit-every 10 &
wait
python $CB --methods spheres_kmeans --obstacle fitted --audit-every 10 &
python $CB --methods spheres_enclosing --obstacle fitted --audit-every 10 &
wait
cd .. && git add tro/results && git commit -q -m "Continuous-boundary baselines on our fitted obstacle fields (step 5)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TByxUTGsR8uPvFhDBetfgi"
for d in 2 4 8 16; do git push -q origin claude/tro-2027-paper-review-wg3qzv && break; sleep $d; done
echo FITTED_DONE
